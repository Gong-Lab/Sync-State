function [C, tvec, fvec] = spike_coherence_multitaper(spikeTimes1, spikeTimes2, fs, timeRange, freqRange, varargin)
% SPIKE_COHERENCE_MULTITAPER_MATLAB
% Computes multitaper coherence spectrogram between two spike trains
% using pure MATLAB functions (no Chronux).
%
% Inputs:
%   spikeTimes1, spikeTimes2: spike time vectors (s) or cell arrays (trials)
%   fs       : sampling freq (Hz) for binning (e.g. 1000 Hz)
%   timeRange: [start end] (s)
%   freqRange: [fmin fmax] (Hz)
%
% Optional name-value pairs:
%   'Window'  : [winLength stepSize] in seconds (default [0.5 0.05])
%   'NW'      : time-half bandwidth product (default 3)
%   'K'       : number of tapers (default 5)
%   'Plot'    : true/false (default true)
%
% Outputs:
%   C    : coherence spectrogram (freq x time)
%   tvec : center times of each window (s)
%   fvec : frequencies (Hz)

    %% Parse inputs
    p = inputParser;
    addParameter(p,'Window',[0.5 0.05]);
    addParameter(p,'NW',3);
    addParameter(p,'K',5);
    addParameter(p,'Plot',true);
    parse(p,varargin{:});
    opts = p.Results;

    %% Crop spike trains to timeRange
    t_start = timeRange(1);
    t_end = timeRange(2);
    spk1 = crop_spikes(spikeTimes1, t_start, t_end);
    spk2 = crop_spikes(spikeTimes2, t_start, t_end);

    %% If cell input, concatenate trials (simple average approach)
    if iscell(spk1)
        spk1 = cell2mat(cellfun(@(x) x(:)', spk1, 'UniformOutput', false));
    end
    if iscell(spk2)
        spk2 = cell2mat(cellfun(@(x) x(:)', spk2, 'UniformOutput', false));
    end

    %% Define parameters
    winLength = opts.Window(1);
    stepSize  = opts.Window(2);
    NW = opts.NW;
    K = opts.K;

    % Time bins for the whole time range (at sampling rate fs)
    edges = t_start:1/fs:t_end;
    nBins = length(edges)-1;

    % Compute number of windows
    winSamples = round(winLength * fs);
    stepSamples = round(stepSize * fs);
    nWins = floor((nBins - winSamples)/stepSamples) + 1;

    % Precompute DPSS tapers
    [tapers, ~] = dpss(winSamples, NW, K);

    % Prepare output
    C_spec = zeros(length(1:winSamples/2+1), nWins);  % prealloc freq x time

    fvec_full = linspace(0, fs/2, winSamples/2+1);
    fInd = fvec_full >= freqRange(1) & fvec_full <= freqRange(2);
    fvec = fvec_full(fInd);

    %% Bin spikes into binary spike trains
    spikeTrain1 = histcounts(spk1, edges)*fs;
    spikeTrain2 = histcounts(spk2, edges)*fs;

    %% Smooth with Gaussian kernel
%     sigma=0.1;
%     g = fspecial('gaussian', [1 round(6*sigma*fs)], sigma*fs);
%     spikeTrain1 = conv(spikeTrain1, g, 'same');
%     spikeTrain2 = conv(spikeTrain2, g, 'same');

    %% Compute coherence window by window
    for w = 1:nWins
        idx = (1:winSamples) + (w-1)*stepSamples;

        x1 = spikeTrain1(idx);
        x2 = spikeTrain2(idx);

        % Detrend window data
        x1 = x1 - mean(x1);
        x2 = x2 - mean(x2);

        % Multitaper spectral estimation
        Sxx = zeros(length(fvec_full),1);
        Syy = zeros(length(fvec_full),1);
        Sxy = zeros(length(fvec_full),1);

        for k = 1:K
            x1t = x1(:) .* tapers(:,k);
            x2t = x2(:) .* tapers(:,k);

            X1 = fft(x1t);
            X2 = fft(x2t);

            % Keep positive freqs
            X1 = X1(1:length(fvec_full));
            X2 = X2(1:length(fvec_full));

            Sxx = Sxx + (X1 .* conj(X1));
            Syy = Syy + (X2 .* conj(X2));
            Sxy = Sxy + (X1 .* conj(X2));
        end

        % Average over tapers
        Sxx = Sxx / K;
        Syy = Syy / K;
        Sxy = Sxy / K;

        % Coherence (only frequencies in freqRange)
        C_spec(:,w) = abs(Sxy).^2 ./ (Sxx .* Syy);
    end

    % Select frequency range of interest
    C = C_spec(fInd, :);

    % Compute time vector centered in each window
    tvec = t_start + ((0:nWins-1)*stepSamples + winSamples/2)/fs;

    %% Plot
    if opts.Plot
        figure;
        imagesc(tvec, fvec, C);
        axis xy;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('Multitaper Coherence Spectrogram (MATLAB)');
        colorbar;
        caxis([0 1]);
    end
end

function out = crop_spikes(spk, t0, t1)
    if iscell(spk)
        out = cellfun(@(x) x(x>=t0 & x<=t1), spk, 'UniformOutput', false);
    else
        spk = spk(spk>=t0 & spk<=t1);
        out = spk;
    end
end
