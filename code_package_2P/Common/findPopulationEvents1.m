function [sigWindows, thr, popCountSmooth, edges] = findPopulationEvents1(eventTimes, binSize, winSize, nShuffle, alpha, doPlot)
% findPopulationEvents identifies short time windows with significant 
% coincident events across multiple event trains using a shuffle test.
%
% INPUTS:
%   eventTimes : cell array {N x 1}, each cell = vector of event times (s)
%   binSize    : scalar, bin width (s)
%   winSize    : scalar, optional smoothing window (s)
%   nShuffle   : number of shuffle iterations
%   alpha      : significance level (e.g. 0.01)
%   doPlot     : logical, whether to plot (default = true)
%
% OUTPUTS:
%   sigWindows     : cell array of [start, end] times for significant bursts
%   thr            : threshold of population activity (from shuffle)
%   popCountSmooth : smoothed population activity time series
%   edges          : time bin edges corresponding to popCountSmooth
%
% Example:
%   eventTimes = {rand(1,30)*10, rand(1,35)*10, rand(1,40)*10};
%   [sigWindows, thr] = findPopulationEvents(eventTimes, 0.01, 0.05, 1000, 0.01, true);

%% Check inputs
if nargin < 6, doPlot = true; end
if nargin < 5, alpha = 0.01; end
if nargin < 4, nShuffle = 1000; end
if nargin < 3, winSize = 0; end

%% Remove empty cells and ensure row vectors
eventTimes = eventTimes(~cellfun(@isempty, eventTimes));
eventTimes = cellfun(@(x) sort(x(:)'), eventTimes, 'UniformOutput', false);
if isempty(eventTimes)
    error('All eventTrains are empty.');
end

%% Define time limits and binning
allTimes = [eventTimes{:}];
tmin = min(allTimes);
tmax = max(allTimes);
edges = tmin:binSize:tmax;
nbins = numel(edges) - 1;
counts = zeros(numel(eventTimes), nbins);

for i = 1:numel(eventTimes)
    counts(i,:) = histcounts(eventTimes{i}, edges);
end

%% Compute population activity
popCount = sum(counts, 1);

%% Optional smoothing
if winSize > binSize
    winBins = round(winSize / binSize);
    smoothKern = ones(1, winBins) / winBins;
    popCountSmooth = conv(popCount, smoothKern, 'same');
else
    popCountSmooth = popCount;
end

%% Shuffle test (ISI randomization per train)
shuffleMax = zeros(nShuffle, 1);
dur = tmax - tmin;

for s = 1:nShuffle
    shCounts = zeros(size(counts));

    for i = 1:numel(eventTimes)
        t = eventTimes{i};
        if numel(t) < 2
            % Too few events to shuffle
            shTimes = t;
        else
            isi = diff(t);
            permISI = isi(randperm(length(isi))); % random permutation
            shTimes = [t(1), t(1) + cumsum(permISI)];
            
            % Rescale or wrap within original duration
            if shTimes(end) > tmax
                shTimes = mod(shTimes - tmin, dur) + tmin;
            end
        end

        shCounts(i,:) = histcounts(shTimes, edges);
    end

    shPop = sum(shCounts, 1);
    if winSize > binSize
        shPop = conv(shPop, smoothKern, 'same');
    end
    shuffleMax(s) = max(shPop);
end

%% Determine significance threshold
thr = prctile(shuffleMax, 100*(1 - alpha));

%% Find significant windows
sigGroups = bwconncomp(popCountSmooth > thr);
sigWindows = cellfun(@(idx) [edges(idx(1)), edges(idx(end)+1)], ...
                     sigGroups.PixelIdxList, 'UniformOutput', false);

%% Plot if requested
if doPlot
    figure; hold on
    plot(edges(1:end-1), popCountSmooth, 'k', 'LineWidth', 1.5);
    yline(thr, 'r--', 'LineWidth', 1.2);
    xlabel('Time (s)'); ylabel('Population activity');
    title('Significant population events (ISI shuffle null)');
    for w = 1:numel(sigWindows)
        x = sigWindows{w};
        fill([x(1) x(2) x(2) x(1)], [0 0 max(popCountSmooth) max(popCountSmooth)], ...
             [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
    hold off
end

fprintf('Found %d significant population events.\n', numel(sigWindows));

end
