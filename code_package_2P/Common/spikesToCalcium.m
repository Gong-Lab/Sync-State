function [calcium, t] = spikesToCalcium(spikeTimeArray, fs, tRange, tau)
% Convert cell array of spike times (multiple neurons) to calcium traces
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
calcium = zeros(numNeurons, numTimePoints);

% Exponential kernel (causal)
kernelLength = round(5 * tau * fs);
kernel = exp(-((0:kernelLength-1)/fs)/tau);

for i = 1:numNeurons
    spikeVec = zeros(1, numTimePoints);
    if ~isempty(spikeTimeArray{i})
        idx = round((spikeTimeArray{i} - tRange(1)) * fs) + 1;
        idx = idx(idx > 0 & idx <= numTimePoints);
        spikeVec(idx) = 1;
    end
    cal_temp = conv(spikeVec, kernel, 'full');
    cal_temp(length(t)+1:end)=[];
    calcium(i,:)=cal_temp;
end

end

% Example usage:
% fs = 100; tRange = [0 10]; tau = 0.7;
% spikeTimes = cell(3,1);
% for i = 1:3
%     spikeTimes{i} = sort(rand(20+i*5,1)*10);
% end
% [calcium, t] = spikesToCalcium(spikeTimes, fs, tRange, tau);
% plot(t, calcium'); xlabel('Time (s)'); ylabel('Calcium