% filepath: c:\Users\YL\Code\vscode101\plot_spike_trains.m
function plot_spike_trains(spike_times_cell, tlim, neuron_indices)
%PLOT_SPIKE_TRAINS Visualize spike trains for multiple neurons efficiently with compact spacing.
%   spike_times_cell: cell array, each cell contains spike times (in seconds) for one neuron
%   tlim: [t_start t_end], time window to plot (optional)
%   neuron_indices: vector of neuron indices to plot (optional)

if nargin < 2 || isempty(tlim)
    tlim = [];
end
if nargin < 3 || isempty(neuron_indices)
    neuron_indices = 1:numel(spike_times_cell);
end

num_neurons = length(neuron_indices);
spike_x = [];
spike_y = [];

spacing = 1; % Reduce spacing between neurons (default was 1, now 0.2)

for idx = 1:num_neurons
    n = neuron_indices(idx);
    st = spike_times_cell{n};
    if ~isempty(tlim)
        st = st(st >= tlim(1) & st <= tlim(2));
    end
    if ~isempty(st)
        spike_x = [spike_x; st(:)];
        spike_y = [spike_y; idx*spacing*ones(numel(st),1)];
    end
end

figure;
plot(spike_x, spike_y, 'k|', 'MarkerSize', 12);
xlabel('Time (s)');
ylabel('Neuron');
yticks(spacing*(1:num_neurons));
yticklabels(neuron_indices);
% ylim([0.5*spacing, (num_neurons+0.5)*spacing]);
title('Spike Trains for Multiple Neurons');
if ~isempty(tlim)
    xlim(tlim);
end
set(gcf, 'Color', 'w');
end