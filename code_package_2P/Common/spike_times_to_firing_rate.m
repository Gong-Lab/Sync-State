% filepath: c:\Users\YL\Code\vscode101\spike_times_to_firing_rate.m
function [t_vec, inst_rate_mat] = spike_times_to_firing_rate(spike_times_cell, fs, sigma, t_range)
%SPIKE_TIMES_TO_FIRING_RATE Converts spike times to instantaneous firing rates for multiple neurons.
%   spike_times_cell: cell array, each cell contains spike times (in seconds) for one neuron
%   fs: sampling frequency for output rate (Hz)
%   sigma: standard deviation of Gaussian (in seconds)
%   t_range: [t_start t_end], time range for output (in seconds)
%   t_vec: time vector (in seconds)
%   inst_rate_mat: matrix of firing rates (neurons x time)

if nargin < 2 || isempty(fs)
    fs = 1000;
end
if nargin < 3 || isempty(sigma)
    sigma = 0.05;
end
if nargin < 4 || isempty(t_range)
    % Default: use min/max across all spike times
    all_spikes = cell2mat(spike_times_cell(:));
    t_range = [min(all_spikes)-3*sigma, max(all_spikes)+3*sigma];
end

t_vec = t_range(1):1/fs:t_range(2);
num_neurons = numel(spike_times_cell);
inst_rate_mat = zeros(num_neurons, numel(t_vec));

% Create Gaussian kernel
kernel_width = round(6*sigma*fs); % +/- 3 sigma
x = -kernel_width:kernel_width;
gauss_kernel = exp(-0.5*(x/(sigma*fs)).^2);
gauss_kernel = gauss_kernel / (sum(gauss_kernel) * (1/fs)); % Normalize for rate in Hz

% Preallocate sparse spike train matrix for efficiency
spike_train_mat = sparse(num_neurons, numel(t_vec));
for n = 1:num_neurons
    spike_times = spike_times_cell{n};
    if ~isempty(spike_times)
        spike_idx = round((spike_times - t_range(1)) * fs) + 1;
        spike_idx = spike_idx(spike_idx >= 1 & spike_idx <= numel(t_vec));
        spike_train_mat(n, spike_idx) = 1;
    end
end

% Use conv2 for efficient convolution across all neurons
inst_rate_mat = conv2(full(spike_train_mat), gauss_kernel, 'same');

end