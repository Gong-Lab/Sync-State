% filepath: c:\Users\YL\Code\vscode101\jitter_spike_times.m
function jittered_spike_times = jitter_spike_times(spike_times, t_range, bin_size)
%JITTER_SPIKE_TIMES Randomly jitters spike times within bins.
%   spike_times: vector of spike times (in seconds)
%   t_range: [t_start t_end], time range for jittering (in seconds)
%   bin_size: size of jitter bin (in seconds)
%   jittered_spike_times: vector of jittered spike times

if nargin < 2 || isempty(t_range)
    t_range = [min(spike_times), max(spike_times)];
end

if nargin < 3 || isempty(bin_size)
    bin_size = 0.01; % default 10 ms
end

% Only jitter spikes within the specified time range
spike_times = spike_times(spike_times >= t_range(1) & spike_times <= t_range(2));

% Assign each spike to a bin
bin_edges = t_range(1):bin_size:t_range(2);
[~, ~, bin_idx] = histcounts(spike_times, bin_edges);

% Jitter each spike uniformly within its bin
jittered_spike_times = zeros(size(spike_times));
for i = 1:length(spike_times)
    if bin_idx(i) > 0 && bin_idx(i) < length(bin_edges)
        bin_start = bin_edges(bin_idx(i));
        bin_end = bin_edges(bin_idx(i)+1);
        jittered_spike_times(i) = bin_start + rand() * (bin_end - bin_start);
    else
        jittered_spike_times(i) = spike_times(i); % leave unchanged if outside bins
    end
end

end