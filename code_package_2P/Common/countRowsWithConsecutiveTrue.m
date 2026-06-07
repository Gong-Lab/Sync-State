function [nRows, rowIdx] = countRowsWithConsecutiveTrue(A, N)
%COUNTRROWSWITHCONSECUTIVETRUE Count rows with at least N consecutive true values.
%
%   [nRows, rowIdx] = countRowsWithConsecutiveTrue(A, N)
%
%   Inputs:
%       A - logical or numeric 2D array
%       N - positive integer, number of consecutive true values to detect
%
%   Outputs:
%       nRows  - number of rows containing at least one sequence of N consecutive true values
%       rowIdx - logical row vector indicating those rows (true = has run)

    % Convert to logical if numeric
    if isnumeric(A)
        A = A ~= 0;
    elseif ~islogical(A)
        error('Input A must be logical or numeric.');
    end

    if ~isscalar(N) || N <= 0
        error('N must be a positive integer scalar.');
    end

    % Convolution kernel of length N
    kernel = ones(1, N);

    % For each row, check if any segment sums to N
    rowIdx = any(conv2(double(A), kernel, 'same') >= N, 2);

    % Count number of rows satisfying the condition
    nRows = sum(rowIdx);
end

