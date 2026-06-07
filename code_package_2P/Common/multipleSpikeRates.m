


function [rates, t] = multipleSpikeRates(spikeTimeArray, fs, kernelWidth, tRange, varargin)
    % Parse inputs
    p = inputParser;
    addRequired(p, 'spikeTimeArray', @iscell);
    addRequired(p, 'fs', @isnumeric);
    addRequired(p, 'kernelWidth', @isnumeric);
    addRequired(p, 'tRange', @isnumeric);
    addParameter(p, 'plotFlag', false, @islogical);
    addParameter(p, 'kernelType', 'gaussian', ...
        @(x) ischar(x) && any(strcmp(x, {'gaussian', 'average'})));
    
    parse(p, spikeTimeArray, fs, kernelWidth, tRange, varargin{:});
    plotFlag = p.Results.plotFlag;
    kernelType = p.Results.kernelType;
    
    
    % Pre-allocate arrays
    t = tRange(1):1/fs:tRange(2);
    numNeurons = length(spikeTimeArray);
    numTimePoints = length(t);
    rates = zeros(numNeurons, numTimePoints);
    
    % Create Gaussian kernel once
    kernelSize = 2 * round(kernelWidth * fs/2) + 1;
    sigma = kernelSize / 6;
    kernel = fspecial('gaussian', [1 kernelSize], sigma);
    kernel = kernel / sum(kernel) * fs;
    
    % Vectorized spike binning
    spikeMat = false(numNeurons, numTimePoints);
    for i = 1:numNeurons
        if ~isempty(spikeTimeArray{i})
            spike_idx = round((spikeTimeArray{i} - tRange(1)) * fs) + 1;
            valid_idx = spike_idx > 0 & spike_idx <= numTimePoints;
            spikeMat(i, spike_idx(valid_idx)) = true;
        end
    end
    
    % Efficient convolution using FFT
    kernelPadded = zeros(1, numTimePoints);
    kernelPadded(1:kernelSize) = kernel;
    kernelFFT = fft(kernelPadded);
    rates = ifft(fft(spikeMat, [], 2) .* kernelFFT, [], 2, 'symmetric');
    
    % Optional plotting
    if plotFlag
        figure('Position', [100 100 800 600]);
        
        % Raster plot
        subplot(2,1,1);
        hold on;
        [neuronIdx, timeIdx] = find(spikeMat);
        scatter(t(timeIdx), neuronIdx, 2, 'k', 'filled');
        ylim([0.5 numNeurons+0.5]);
        ylabel('Neuron');
        title('Multi-neuron Spike Trains and Firing Rates');
        hold off;
        
        % Firing rate plot
        subplot(2,1,2);
        colors = jet(numNeurons);
        plot(t, rates', 'LineWidth', 1.5);
        colororder(colors);
        xlabel('Time (s)');
        ylabel('Firing Rate (Hz)');
        legend(arrayfun(@(x) sprintf('Neuron %d', x), 1:numNeurons, 'UniformOutput', false), ...
               'Location', 'eastoutside');
        grid on;
    end
end