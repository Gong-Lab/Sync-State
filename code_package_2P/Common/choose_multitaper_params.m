function [NW, K] = choose_multitaper_params(T, target_df)
% CHOOSE_MULTITAPER_PARAMS - Suggest NW and K for desired frequency resolution
%
% Inputs:
%   T         : Window length (s) used in multitaper analysis
%   target_df : Desired frequency resolution (Hz)
%
% Outputs:
%   NW : Suggested time-half bandwidth product
%   K  : Suggested number of DPSS tapers
%
% Formula: delta_f ≈ 2*NW / T
%   => NW ≈ (target_df * T) / 2
%   => K ≈ 2NW - 1 (Slepian's rule)

    % Compute NW from target delta_f
    NW = max(1.5, (target_df * T) / 2);  % keep NW >= 1.5
    NW = round(NW);                      % use integer
    
    % Compute K from NW
    K = max(3, round(2*NW - 1));         % at least 3 tapers

    fprintf('For window T = %.2f s and desired Δf = %.1f Hz:\n', T, target_df);
    fprintf('  Suggested NW = %d, K = %d\n', NW, K);
end
