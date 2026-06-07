function [calcium, t] = spikesToCalcium_vfilter(spikeTimeArray, fs, tRange, tau)
% Efficiently convert cell array of spike times (multiple neurons) to calcium traces
% Inputs:
%   spikeTimeArray - cell array, each cell contains spike times (seconds) for one neuron
%   fs - sample rate (Hz)
%   tRange - [startTime endTime] in seconds
%   tau - decay time constant (seconds)
%
% Outputs:
%   calcium - matrix (neurons x time) of calcium traces
%   t - time vector (seconds)

if nargin < 4
    tau = 0.5; % default decay time constant (seconds)
end

t = tRange(1):1/fs:tRange(2);
numNeurons = numel(spikeTimeArray);
numTimePoints = numel(t);

% Exponential kernel (causal)
kernelLength = round(5 * tau * fs);
kernel = exp(-((0:kernelLength-1)/fs)/tau);

% Preallocate spike matrix for efficiency
spikeMat = zeros(numNeurons, numTimePoints);

for i = 1:numNeurons
    if ~isempty(spikeTimeArray{i})
        idx = round((spikeTimeArray{i} - tRange(1)) * fs) + 1;
        idx = idx(idx > 0 & idx <= numTimePoints);
        spikeMat(i, idx) = 1;
    end
end

% Convolve all neurons at once using filter (faster than conv for causal kernel)
calcium = filter(kernel, 1, spikeMat, [], 2);

end

% Example usage:
% fs = 100; tRange = [0 10]; tau = 0.7;
% spikeTimes = cell(3,1);
% for i = 1:3
%     spikeTimes{i} = sort(rand(20+i*5,1)*10);
% end
% [calcium, t] = spikesToCalcium(spikeTimes, fs, tRange, tau);
% plot(t, calcium'); xlabel('Time (s)'); ylabel('Calcium signal');