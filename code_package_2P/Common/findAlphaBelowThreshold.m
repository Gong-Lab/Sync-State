function alpha_star = findAlphaBelowThreshold(y_shfl, target_globalAlpha, alpha_min, alpha_max, tol, side, n_b)

if nargin < 3 || isempty(alpha_min), alpha_min = 1e-4; end
if nargin < 4 || isempty(alpha_max), alpha_max = 0.05; end
if nargin < 5 || isempty(tol), tol = 1e-6; end
if nargin < 6 || isempty(side), side = 'both'; end
if nargin < 7 || isempty(n_b), n_b = 3; end

% Sanity check
f_min=getShflFrBrk(y_shfl, alpha_min, side, n_b);
if  f_min>= target_globalAlpha
    error('Even smallest alpha gives getfr >= threshold');
end
f_max=getShflFrBrk(y_shfl, alpha_max, side, n_b);
if f_max < target_globalAlpha
    alpha_star = alpha_max;
    return;
end

% Bisection search
low = alpha_min;
high = alpha_max;

while (high - low) > tol
    mid = (low + high) / 2;
    f_mid = getShflFrBrk(y_shfl, mid, side, n_b);
    if f_mid < target_globalAlpha
        low = mid;   % we can go higher
    else
        high = mid;  % too high, go lower
    end
end

alpha_star = low;
end
