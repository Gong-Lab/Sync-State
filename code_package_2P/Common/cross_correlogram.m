% filepath: c:\Users\YL\Code\vscode101\cross_correlogram.m
function [lags, ccg] = cross_correlogram(spike_times1, spike_times2, bin_size, max_lag, t_range)
%CROSS_CORRELOGRAM Calculate cross-correlogram of two neurons' spike times.
%   [lags, ccg] = cross_correlogram(spike_times1, spike_times2, bin_size, max_lag, t_range)
%   spike_times1, spike_times2: vectors of spike times (in seconds)
%   bin_size: bin size for correlogram (in seconds)
%   max_lag: maximum lag to compute (in seconds)
%   t_range: [t_start t_end], time range to consider (in seconds)
%   lags: vector of lag times (in seconds)
%   ccg: cross-correlogram counts

% Restrict spikes to time range
if nargin < 5 || isempty(t_range)
    t_start = min([spike_times1(:); spike_times2(:)]);
    t_end = max([spike_times1(:); spike_times2(:)]);
else
    t_start = t_range(1);
    t_end = t_range(2);
end
spike_times1 = spike_times1(spike_times1 >= t_start & spike_times1 <= t_end);
spike_times2 = spike_times2(spike_times2 >= t_start & spike_times2 <= t_end);

% Compute all pairwise time differences
dt = [];
for i = 1:length(spike_times1)
    diffs = spike_times2 - spike_times1(i);
    diffs = diffs(abs(diffs) <= max_lag);
    dt = [dt; diffs(:)]; % Ensure column vector
end

% Bin the time differences
edges = -max_lag:bin_size:max_lag;
ccg = histcounts(dt, edges);

% Output lags as bin centers
lags = edges(1:end-1) + bin_size/2;

end