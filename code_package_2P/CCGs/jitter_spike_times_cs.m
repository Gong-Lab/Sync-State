% filepath: c:\Users\YL\Code\vscode101\jitter_spike_times.m
function jittered_spike_times = jitter_spike_times_cs(spike_times, t_range, bin_size)
%JITTER_SPIKE_TIMES Circularly shifts spike times within each bin by a random time step.
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

jittered_spike_times = spike_times; % initialize

for b = 1:length(bin_edges)-1
    idx_in_bin = find(bin_idx == b);
    n_spikes = numel(idx_in_bin);
    if n_spikes > 1
        rel_times = spike_times(idx_in_bin) - bin_edges(b);
        t_shift=bin_size*rand;
        rel_times_shifted=rel_times+t_shift;
        rel_times_shifted(rel_times_shifted>bin_size)=rel_times_shifted(rel_times_shifted>bin_size)-bin_size;
        jittered_spike_times(idx_in_bin) = bin_edges(b) + rel_times_shifted;
    end
end

end