function [corrs, time_centers] = sliding_corr(fr1, fr2, time, win_size, step_size)
% Compute correlation between two firing rate vectors over overlapping windows
% Preallocate output arrays for efficiency and exclude NaN results

win_starts = time(1):step_size:(time(end)-win_size);
num_windows = length(win_starts);

corrs_tmp = zeros(1, num_windows);
time_centers_tmp = zeros(1, num_windows);
valid_count = 0;

for i = 1:num_windows
    win_start = win_starts(i);
    win_end = win_start + win_size;
    idx = find(time >= win_start & time < win_end);
    if length(idx) >= 2
        valid_count = valid_count + 1;
        corrs_tmp(valid_count) = corr(fr1(idx)', fr2(idx)');
        time_centers_tmp(valid_count) = mean(time(idx));
    end
end

% Trim unused preallocated space
corrs = corrs_tmp(1:valid_count);
time_centers = time_centers_tmp(1:valid_count);