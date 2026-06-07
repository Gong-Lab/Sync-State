function [win_times, win_centers] = moving_windows(time_range, win_length, step_size)
    % moving_windows - Create moving windows within a given time range
    %
    % INPUTS:
    %   time_range  - [t_start, t_end] (s)
    %   win_length  - window length (s)
    %   step_size   - step size (s)
    %
    % OUTPUTS:
    %   win_times   - Nx2 array, [start_time end_time] per row
    %   win_centers - Nx1 array, center time of each window

    % Unpack
    t_start = time_range(1);
    t_end   = time_range(2);

    % Generate start times
    start_times = t_start:step_size:(t_end - win_length);

    % Generate end times
    end_times = start_times + win_length;

    % Combine
    win_times = [start_times(:), end_times(:)];

    % Compute centers
%     win_centers = mean(win_times, 2);
    win_centers = win_times(:, 2);
end
