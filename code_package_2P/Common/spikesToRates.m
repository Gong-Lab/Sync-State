function [rates, t] = spikesToRates(spikeTimeArray, fs, kernelWidth, tRange, kernelType)
% Convert cell array of spike times to instantaneous firing rates using convolution
% Inputs:
%   spikeTimeArray - cell array, each cell contains spike times (in seconds) for one neuron
%   fs - sample rate in Hz
%   kernelWidth - kernel width in seconds
%   tRange - [startTime endTime] in seconds
%   kernelType - 'gaussian' or 'average' (default: 'gaussian')
%
% Outputs:
%   rates - matrix (neurons x time) of firing rates in Hz
%   t - time vector in seconds

if nargin < 5
    kernelType = 'gaussian';
end

t = tRange(1):1/fs:tRange(2);
numNeurons = numel(spikeTimeArray);
numTimePoints = numel(t);

% Create kernel
kernelSize = 2 * round(kernelWidth * fs/2) + 1;
switch lower(kernelType)
    case 'gaussian'
        sigma = kernelSize / 6;
        kernel = fspecial('gaussian', [1 kernelSize], sigma);
        kernel = kernel / (sum(kernel)/fs); % normalize for Hz
    case 'average'
        kernel = ones(1, kernelSize) / kernelSize * fs;
    otherwise
        error('kernelType must be ''gaussian'' or ''average''');
end

% % Preallocate spike matrix for efficiency
% spikeMat = zeros(numNeurons, numTimePoints);
% 
% for i = 1:numNeurons
%     if ~isempty(spikeTimeArray{i})
%         idx = round((spikeTimeArray{i} - tRange(1)) * fs) + 1;
%         idx = idx(idx > 0 & idx <= numTimePoints);
%         spikeMat(i, idx) = 1;
%     end
% end
% 
% rates = filter(kernel, 1, spikeMat, [], 2);

rates=zeros(numNeurons,numTimePoints);
for i = 1:numNeurons
    spikes = zeros(1, numTimePoints);
    if ~isempty(spikeTimeArray{i})
        idx = round((spikeTimeArray{i} - tRange(1)) * fs) + 1;
        idx = idx(idx > 0 & idx <= numTimePoints);
        spikes(idx) = 1;
    end
    rates(i,:) = conv(spikes, kernel, 'same');
end

end

% Example usage:
% fs = 1000;
% kernelWidth = 0.1;
% tRange = [0 10];
% spikeTimes = cell(3,1);
% for i = 1:3
%     spikeTimes{i} = sort(rand(50+i*10,1)*10);
% end
% [rates, t] = spikesToRates(spikeTimes, fs, kernelWidth, tRange, 'gaussian');