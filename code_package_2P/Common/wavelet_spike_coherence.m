function [Cxy, t_wav, f, signifMask] = wavelet_spike_coherence( ...
    spikeTimes1, spikeTimes2, fs, timeRange, freqRange, varargin)
% WAVELET_SPIKE_COHERENCE_FINAL
% Computes wavelet-based time-frequency coherence between two spike trains:
%   - Proper cross-spectral averaging (fixes "all 1s" issue)
%   - Adaptive frequency resolution (CWT)
%   - Choice of wavelet (Morlet, Morse, Paul)
%   - Optional surrogate-based significance testing
%
% Usage:
%   [Cxy, t_wav, f, signifMask] = wavelet_spike_coherence_final(spk1, spk2, fs, timeRange, freqRange, ...
%        'Wavelet','amor','SmoothingTime',10,'SmoothingFreq',3,'NumSurrogates',50,'Alpha',0.05,'Plot',true)

    %% Parse inputs
    p = inputParser;
    addParameter(p,'Wavelet','amor');          % Morlet wavelet by default
    addParameter(p,'SmoothingTime',10);        % temporal smoothing
    addParameter(p,'SmoothingFreq',3);         % frequency smoothing
    addParameter(p,'Sigma',0.01);              % Gaussian smoothing of firing rate (s)
    addParameter(p,'NumSurrogates',0);         % no significance testing by default
    addParameter(p,'Alpha',0.05);              % significance level
    addParameter(p,'Plot',true);               % plot by default
    parse(p,varargin{:});
    opts = p.Results;

    %% Parameters
    t_start = timeRange(1);
    t_end   = timeRange(2);
    binWidth = 1/fs;

    %% Time vector
    t = t_start:binWidth:t_end;

    %% Bin spike trains into firing rate
    rate1 = histcounts(spikeTimes1(spikeTimes1 >= t_start & spikeTimes1 <= t_end), [t t_end+binWidth]) / binWidth;
    rate2 = histcounts(spikeTimes2(spikeTimes2 >= t_start & spikeTimes2 <= t_end), [t t_end+binWidth]) / binWidth;

    %% Smooth firing rates with Gaussian
    kernelSize = max(1, round(6 * opts.Sigma / binWidth));
    g = gausswin(kernelSize, 3)';   % use gausswin instead of fspecial
    g = g / sum(g);                 % normalize kernel
    rate1 = conv(rate1, g, 'same');
    rate2 = conv(rate2, g, 'same');

    %% Continuous Wavelet Transform
    [W1, f_all] = cwt(rate1, fs, opts.Wavelet);
    [W2, ~]     = cwt(rate2, fs, opts.Wavelet);

    %% Restrict to desired frequency band
    freq_idx = (f_all >= freqRange(1) & f_all <= freqRange(2));
    f = f_all(freq_idx);
    W1 = W1(freq_idx,:);
    W2 = W2(freq_idx,:);

    %% Cross- and auto-spectra
    Pxx = abs(W1).^2;
    Pyy = abs(W2).^2;
    Pxy = W1 .* conj(W2);

    %% === Proper Cross-Spectral Averaging ===
    kernel2D = fspecial('average', [ceil(opts.SmoothingFreq) ceil(opts.SmoothingTime)]);
    Pxx_s = conv2(Pxx, kernel2D, 'same');
    Pyy_s = conv2(Pyy, kernel2D, 'same');
    Pxy_real = conv2(real(Pxy), kernel2D, 'same');
    Pxy_imag = conv2(imag(Pxy), kernel2D, 'same');
    Pxy_s = Pxy_real + 1i*Pxy_imag;

    %% Magnitude-squared wavelet coherence
    epsVal = 1e-12;
    Cxy = abs(Pxy_s).^2 ./ ((Pxx_s .* Pyy_s) + epsVal);

    %% Initialize significance mask
    signifMask = ones(size(Cxy));

    %% === Optional Surrogate-Based Significance Testing ===
    if opts.NumSurrogates > 0
        surroVals = zeros([size(Cxy), opts.NumSurrogates]);
        for s = 1:opts.NumSurrogates
            shift = randi(length(rate2));
            rate2s = circshift(rate2, shift);
            [W2s, ~] = cwt(rate2s, fs, opts.Wavelet);
            W2s = W2s(freq_idx,:);
            Pyy_surr = abs(W2s).^2;
            Pxy_surr = W1 .* conj(W2s);

            % smooth surrogate
            Pyy_surr = conv2(Pyy_surr, kernel2D, 'same');
            Pxy_sr = conv2(real(Pxy_surr), kernel2D, 'same');
            Pxy_si = conv2(imag(Pxy_surr), kernel2D, 'same');
            Pxy_surr = Pxy_sr + 1i*Pxy_si;

            % compute coherence surrogate
            surroVals(:,:,s) = abs(Pxy_surr).^2 ./ ((Pxx_s .* Pyy_surr) + epsVal);
        end

        % threshold from surrogate distribution
        thresh = prctile(surroVals, 100*(1-opts.Alpha), 3);
        signifMask(Cxy < thresh) = NaN;
    end

    %% Time vector
    t_wav = t;

    %% === Plot ===
    if opts.Plot
        figure;
        imagesc(t_wav, f, Cxy .* signifMask);
        axis xy;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title(sprintf('Wavelet Coherence (%s)', opts.Wavelet));
        colorbar;
        caxis([0 1]);
    end
end
