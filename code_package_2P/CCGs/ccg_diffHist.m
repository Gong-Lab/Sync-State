% filepath: c:\Users\YL\Code\vscode101\cross_correlogram_triangle_norm.m
function [lags, ccg_norm] = ccg_diffHist(spike_times1, spike_times2, bin_size, max_lag, t_range)
%CROSS_CORRELOGRAM_TRIANGLE_NORM Cross-correlogram with triangular and rate normalization.
%   [lags, ccg_norm] = cross_correlogram_triangle_norm(spike_times1, spike_times2, bin_size, max_lag, t_range)
%   spike_times1, spike_times2: vectors of spike times (in seconds)
%   bin_size: bin size for correlogram (in seconds)
%   max_lag: maximum lag to compute (in seconds)
%   t_range: [t_start t_end], time range to consider (in seconds)
%   lags: vector of lag times (in seconds)
%   ccg_norm: normalized cross-correlogram

if nargin < 5 || isempty(t_range)
    t_start = min([spike_times1(:); spike_times2(:)]);
    t_end = max([spike_times1(:); spike_times2(:)]);
else
    t_start = t_range(1);
    t_end = t_range(2);
end
duration = t_end - t_start;

% Restrict spikes to time range
spike_times1 = spike_times1(spike_times1 >= t_start & spike_times1 <= t_end);
spike_times2 = spike_times2(spike_times2 >= t_start & spike_times2 <= t_end);

% Compute all pairwise time differences
dt = [];
for i = 1:length(spike_times1)
    diffs = spike_times2 - spike_times1(i);
    diffs = diffs(abs(diffs) <= max_lag);
    dt = [dt; diffs(:)];
end

% Bin the time differences
edges = -max_lag:bin_size:max_lag;
ccg = histcounts(dt, edges);
lags = edges(1:end-1);
% lags = edges(1:end-1) + bin_size/2;

% Triangular normalization: for each lag, compute the overlap duration
overlap = duration - abs(lags);
overlap(overlap < 0) = 0; % No overlap outside recording

% Firing rates
rate1 = numel(spike_times1) / duration;
rate2 = numel(spike_times2) / duration;
geom_mean_rate = sqrt(rate1 * rate2);

% Normalize: correct for overlap and geometric mean rate
ccg_norm = ccg ./ (overlap * geom_mean_rate);

end