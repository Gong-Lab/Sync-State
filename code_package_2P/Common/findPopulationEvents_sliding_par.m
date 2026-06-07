function [sigWindows, thr, popCountSmooth, tCenters] = findPopulationEvents_sliding_par(eventTimes, winSize, stepSize, nShuffle, alpha, doPlot)
% findPopulationEvents_sliding_par
% Detects significant synchronous population events using sliding windows
% and an ISI-shuffled null distribution, parallelized with parfor.
%
% INPUTS:
%   eventTimes : cell array {N x 1}, each cell = vector of event times (s)
%   winSize    : scalar, sliding window width (s)
%   stepSize   : scalar, step between successive windows (s)
%   nShuffle   : number of shuffles for null distribution
%   alpha      : significance level (default = 0.01)
%   doPlot     : logical, whether to plot (default = true)
%
% OUTPUTS:
%   sigWindows     : cell array of [start, end] times of significant events
%   thr            : threshold value (percentile of shuffle max)
%   popCountSmooth : raw population count trace across sliding windows
%   tCenters       : vector of window center times
%
% Example:
%   eventTimes = {sort(rand(1,30)*10), sort(rand(1,35)*10), sort(rand(1,40)*10)};
%   [sigWindows, thr] = findPopulationEvents_sliding_par(eventTimes, 0.05, 0.005, 1000, 0.01, true);

%% Default arguments
if nargin < 6, doPlot = true; end
if nargin < 5, alpha = 0.01; end
if nargin < 4, nShuffle = 1000; end
if nargin < 3, stepSize = winSize / 10; end

%% Preprocess event trains
eventTimes = eventTimes(~cellfun(@isempty, eventTimes));
eventTimes = cellfun(@(x) sort(x(:)'), eventTimes, 'UniformOutput', false);
if isempty(eventTimes)
    error('All event trains are empty.');
end

%% Basic parameters
allTimes = [eventTimes{:}];
tmin = min(allTimes);
tmax = max(allTimes);
dur = tmax - tmin;

tCenters = tmin + winSize/2 : stepSize : tmax - winSize/2;
nWin = numel(tCenters);
nTrains = numel(eventTimes);

%% --- Population activity in sliding windows ---
popCount = zeros(1, nWin);
for w = 1:nWin
    t1 = tCenters(w) - winSize/2;
    t2 = tCenters(w) + winSize/2;
    for i = 1:nTrains
        popCount(w) = popCount(w) + sum(eventTimes{i} >= t1 & eventTimes{i} < t2);
    end
end
popCountSmooth = popCount;

%% --- Parallel shuffle (ISI randomization) ---
shuffleMax = zeros(nShuffle, 1);

fprintf('Running %d shuffles in parallel...\n', nShuffle);
parfor s = 1:nShuffle
    shPop = zeros(1, nWin);
    for i = 1:nTrains
        t = eventTimes{i};
        if numel(t) < 2
            shTimes = t;
        else
            isi = diff(t);
            permISI = isi(randperm(length(isi))); % random permutation of ISIs
            shTimes = [t(1), t(1) + cumsum(permISI)];
            % wrap if out of range
            if shTimes(end) > tmax
                shTimes = mod(shTimes - tmin, dur) + tmin;
            end
        end

        % Count shuffled events per sliding window
        for w = 1:nWin
            t1 = tCenters(w) - winSize/2;
            t2 = tCenters(w) + winSize/2;
            shPop(w) = shPop(w) + sum(shTimes >= t1 & shTimes < t2);
        end
    end
    shuffleMax(s) = max(shPop);
end

%% --- Compute significance threshold ---
thr = prctile(shuffleMax, 100*(1 - alpha));

%% --- Detect significant population activity ---
above = popCountSmooth > thr;
sigGroups = bwconncomp(above);
sigWindows = cellfun(@(idx) [tCenters(idx(1)) - winSize/2, tCenters(idx(end)) + winSize/2], ...
                     sigGroups.PixelIdxList, 'UniformOutput', false);

%% --- Plot results ---
if doPlot
    figure; hold on
    plot(tCenters, popCountSmooth, 'k', 'LineWidth', 1.5);
    yline(thr, 'r--', 'LineWidth', 1.3);
    xlabel('Time (s)'); ylabel('Population count');
    title('Significant population events (ISI-shuffle, sliding, parfor)');
    for w = 1:numel(sigWindows)
        x = sigWindows{w};
        fill([x(1) x(2) x(2) x(1)], [0 0 max(popCountSmooth) max(popCountSmooth)], ...
             [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
    hold off
end

fprintf('Found %d significant population events (p < %.3f).\n', numel(sigWindows), alpha);

end
