function [Z, p, C, mu_null, sigma_null] = event_synchrony_index(A, B, tau, nperm)
% EVENT_SYNCHRONY_INDEX  Compute coincidence synchrony index (Z-score) between two event trains.
%
% Usage:
%   [Z, p, C, mu_null, sigma_null] = event_synchrony_index(A, B, tau, nperm)
%
% Inputs:
%   A      - vector of event times (sorted, same units as tau)
%   B      - vector of event times (sorted)
%   tau    - coincidence window (e.g. 0.005 for 5 ms)
%   nperm  - number of surrogate permutations (default 1000)
%
% Outputs:
%   Z          - Z-score of coincidence count relative to null distribution
%   p          - one-sided p-value (probability that null >= observed)
%   C          - observed coincidence count
%   mu_null    - mean of null coincidence counts
%   sigma_null - std of null coincidence counts

if nargin < 4 || isempty(nperm), nperm = 1000; end

A = A(:); B = B(:);

% duration for circular shift
Tmin = min([A;B]);
Tmax = max([A;B]);
Tdur = Tmax - Tmin;

% --- observed coincidences ---
C = count_pairs(A, B, tau);

% --- surrogate coincidences ---
counts_null = zeros(nperm,1);
for k = 1:nperm
    offset = rand * Tdur;  % random circular shift
    Bshift = mod(B - Tmin + offset, Tdur) + Tmin;
    counts_null(k) = count_pairs(A, sort(Bshift), tau);
end

mu_null    = mean(counts_null);
sigma_null = std(counts_null);

if sigma_null == 0
    Z = NaN;
else
    Z = (C - mu_null) / sigma_null;
end

p = (sum(counts_null >= C) + 1) / (nperm + 1);

end

% -----------------------------
function C = count_pairs(A,B,tau)
% count pairs within coincidence window
C = 0;
for i = 1:numel(A)
    % find B indices within window
    lb = find(B >= A(i)-tau, 1, 'first');
    if isempty(lb), continue; end
    j = lb;
    while j <= numel(B) && B(j) <= A(i)+tau
        C = C + 1;
        j = j + 1;
    end
end
end
