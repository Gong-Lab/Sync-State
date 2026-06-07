function [sig_time_centers, corrs, time_centers, p_values] = sliding_corr_shuffle(fr1, fr2, time, win_size, step_size, n_shuffles, alpha)
% Compute correlation and test significance using shuffles
% Returns significant time centers, correlations, all time centers, and p-values

win_starts = time(1):step_size:(time(end)-win_size);
num_windows = length(win_starts);

corrs_tmp = zeros(1, num_windows);
time_centers_tmp = zeros(1, num_windows);
p_values_tmp = zeros(1, num_windows);
valid_count = 0;

for i = 1:num_windows
    win_start = win_starts(i);
    win_end = win_start + win_size;
    idx = find(time >= win_start & time < win_end);
    if length(idx) >= 2
        valid_count = valid_count + 1;
        obs_corr = corr(fr1(idx)', fr2(idx)');
        % Shuffling
        shuffle_corrs = zeros(1, n_shuffles);
        for s = 1:n_shuffles
            fr2_shuff = fr2(idx(randperm(length(idx))));
            shuffle_corrs(s) = corr(fr1(idx)', fr2_shuff');
        end
        % Two-tailed p-value
        p = mean(abs(shuffle_corrs) >= abs(obs_corr));
        corrs_tmp(valid_count) = obs_corr;
        time_centers_tmp(valid_count) = mean(time(idx));
        p_values_tmp(valid_count) = p;
    end
end

corrs = corrs_tmp(1:valid_count);
time_centers = time_centers_tmp(1:valid_count);
p_values = p_values_tmp(1:valid_count);

% Significant time centers
sig_time_centers = time_centers(p_values < alpha);