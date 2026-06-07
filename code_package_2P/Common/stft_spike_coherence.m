function [Cxy, t_stft, f] = stft_spike_coherence(spikeTimes1, spikeTimes2, fs, timeRange, freqRange, doPlot)
% STFT_SPIKE_COHERENCE Computes time-frequency coherence between two neurons
% based on firing rate estimates and STFT with spectral smoothing.
%
% Usage:
%   [Cxy, t_stft, f] = stft_spike_coherence(spikeTimes1, spikeTimes2, fs, timeRange, freqRange, doPlot)
%
% Inputs:
%   spikeTimes1 - spike times of neuron 1 (s)
%   spikeTimes2 - spike times of neuron 2 (s)
%   fs          - sampling rate for firing rate signal (Hz)
%   timeRange   - [t_start t_end], analysis window (s)
%   freqRange   - [f_min f_max], frequency range of interest (Hz)
%   doPlot      - (optional) true to plot (default = true)
%
% Outputs:
%   Cxy         - coherence spectrogram (freq x time)
%   t_stft      - STFT time vector (s)
%   f           - frequency vector (Hz)

    if nargin < 6
        doPlot = true;
    end

    %% Parameters
    t_start = timeRange(1);
    t_end   = timeRange(2);
    binWidth = 1/fs;        % bin width from sampling rate
    sigma = 0.1;           % 10 ms Gaussian smoothing
    smoothWin = 5;          % smoothing across time windows for coherence

    %% Time vector
    t = t_start:binWidth:t_end;

    %% Filter spike times
    spk1 = spikeTimes1(spikeTimes1 >= t_start & spikeTimes1 <= t_end);
    spk2 = spikeTimes2(spikeTimes2 >= t_start & spikeTimes2 <= t_end);

    %% Bin spikes and estimate firing rates
    counts1 = histcounts(spk1, [t t_end+binWidth]);
    counts2 = histcounts(spk2, [t t_end+binWidth]);
    rate1 = counts1 / binWidth;
    rate2 = counts2 / binWidth;

    %% Smooth with Gaussian kernel
    g = fspecial('gaussian', [1 round(6*sigma/binWidth)], sigma/binWidth);
    rate1 = conv(rate1, g, 'same');
    rate2 = conv(rate2, g, 'same');

    %% STFT parameters
    win = hamming(fs*10);
    noverlap = fs*5;
    nfft = fs*30;

    %% Compute STFT
    [Sx, f_all, t_stft] = spectrogram(rate1, win, noverlap, nfft, fs);
    [Sy, ~, ~]          = spectrogram(rate2, win, noverlap, nfft, fs);

    %% Adjust time to analysis window
    t_stft = t_stft + t_start;

    %% Compute cross- and auto-spectra
    Pxx = abs(Sx).^2;
    Pyy = abs(Sy).^2;
    Pxy = Sx .* conj(Sy);

    %% Spectral smoothing across time windows
    Pxx_s = movmean(Pxx, smoothWin, 2);
    Pyy_s = movmean(Pyy, smoothWin, 2);
    Pxy_s = movmean(Pxy, smoothWin, 2);

    %% Coherence
    epsVal = 1e-12;
    Cxy_all = abs(Pxy_s).^2 ./ ((Pxx_s .* Pyy_s) + epsVal);

    %% Restrict to frequency range
    f_min = freqRange(1);
    f_max = freqRange(2);
    freq_idx = (f_all >= f_min & f_all <= f_max);
    f   = f_all(freq_idx);
    Cxy = Cxy_all(freq_idx, :);

    %% Optional plot
    if doPlot
        figure;
        imagesc(t_stft, f, Cxy);
        axis xy;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('STFT-Based Coherence (Firing Rate)');
        colorbar;
        caxis([0 1]);
    end
end
