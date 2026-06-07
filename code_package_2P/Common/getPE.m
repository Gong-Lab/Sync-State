function [pe_locs, peCount, t]=getPE(eventT, tRange_PE, win_size, win_step, smoothC, MinEventDistance, doPlot, shflMethod, shflT_baseline)

if nargin<5
    smoothC=true;
end

if nargin<6
    MinEventDistance=0.2;
end

if nargin<7
    doPlot=true;
end

if nargin<8
    shflMethod='cirShft';
end

if nargin<9
    shflT_baseline=tRange_PE;
end


winBins = round(win_size / win_step);
smoothKern = ones(1, winBins) / winBins;


[win_times, win_centers] = moving_windows(tRange_PE, win_size, win_step);
n_windows=length(win_centers);
t=win_centers;

n_cells=length(eventT);

peCount=zeros(1, n_windows);
for w=1:n_windows
    for i=1:n_cells
        peCount(w) = peCount(w) + sum(eventT{i} >= win_times(w,1) & eventT{i} < win_times(w,2));
    end
end
if smoothC
    peCount = conv(peCount, smoothKern, 'same');
end

%% shuffle
n_shfl = 10000;
i_b1=dsearchn(win_times(:,1), shflT_baseline(1));
i_b2=dsearchn(win_times(:,1), shflT_baseline(2));
winEdges_b = [win_times(i_b1:i_b2,1); win_times(i_b2,2)];  % ensure edges cover all windows
n_windows_b = length(winEdges_b)-1;

%Preallocate
peCount_shfl = zeros(n_shfl, n_windows_b);

%Parallel shuffle loop
switch shflMethod
    case 'cirShft'
        fprintf('Running %d shuffles using circular-shift null...\n', n_shfl);
    case 'perm'
        fprintf('Running %d shuffles using permutation null...\n', n_shfl);
end

parfor s = 1:n_shfl
    % Each worker runs an independent shuffle
    rng(s + 1000, 'twister'); % reproducibility
    
    % Shuffle each cell's event train
    pEvent_shfl = cell(1, n_cells);
    for i = 1:n_cells
        switch shflMethod
            case 'cirShft'
                pEvent_shfl{i} = shflEventTrain_cirShft(eventT{i}, shflT_baseline);
            case 'perm'
                pEvent_shfl{i} = shflEventTrain_perm(eventT{i}, shflT_baseline);
        end
    end

    % Count all shuffled events across cells
    allShufTimes = [pEvent_shfl{:}];
    % Use histogram to count within sliding windows
    counts = histcounts(allShufTimes, winEdges_b);
    
    if smoothC
        counts = conv(counts, smoothKern, 'same');
    end

    peCount_shfl(s, :) = counts;

    % if mod(s, max(1, round(n_shfl/10))) == 0
    %     fprintf('Completed %d/%d shuffles...\n', s, n_shfl);
    % end
end

fprintf('All %d shuffles finished.\n', n_shfl);

%%
% alpha=0.05;
% thr_shfl.ptwise=prctile(ccg_norm_shfl, (1-alpha)*100, 1);

target_globalAlpha=0.05;
alpha_adj=findAlphaBelowThreshold(peCount_shfl, target_globalAlpha, 1/n_shfl, 0.05, [], 'greater', 5);
% alpha_adj=0.05;
thr_=prctile(peCount_shfl, 100*(1 - alpha_adj), 1);

% alpha=0.01;
% thr_=prctile(peCount_shfl, 100*(1 - alpha), 1);

%%
% min_peCount=3;
min_peCount=max([2, median(thr_)+1]);
disp(['minCount:', num2str(min_peCount)]);
[pks,pe_locs]=findpeaks(peCount, win_centers, 'MinPeakHeight',min_peCount, 'MinPeakDistance',MinEventDistance);


%%
if doPlot
    figure;hold on
    plot(win_centers, peCount);
    plot([win_centers(1),win_centers(end)], [min_peCount,min_peCount], 'm--');
    plot(pe_locs, pks, 'r.');
    % plot(win_centers, peCount/n_cells);
    % plot([win_centers(1),win_centers(end)], [min_peCount/n_cells,min_peCount/n_cells], 'm--');
    % plot(pe_locs, pks/n_cells, 'r.');
    % ylim([0, 0.5]);
end