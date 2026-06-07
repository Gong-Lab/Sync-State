function burst_indices=getBurst(spike_times, isi_thr)

if nargin<2
    isi_thr=[0.08, 0.16];
end

isi = diff(spike_times);
in_burst = false;
burst_indices = {};
i = 1;
while i <= length(isi)
    if ~in_burst && isi(i) <= isi_thr(1)
        % Start of a burst
        burst_start_idx = i;
        in_burst = true;
    end
    if in_burst
        % Check for end of burst
        if isi(i) >= isi_thr(2)
            burst_end_idx = i + 1;
            burst_indices{end+1} = burst_start_idx:burst_end_idx;
            in_burst = false;
        elseif i == length(isi)
            % End of spike train, close burst
            burst_end_idx = length(spike_times);
            burst_indices{end+1} = burst_start_idx:burst_end_idx;
            in_burst = false;
        end
    end
    i = i + 1;
end