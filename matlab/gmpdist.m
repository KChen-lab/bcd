function [D] = gmpdist(X, S, th, r)
%GMPDIST(X, S, th, r) calculates generalized Mahalanobis pairwise distance
%   X: data matrix, each row is an instance and each column is a feature
%   S: Covariance matrix
%   th: threshold of eigenvalue. any eigenvalue smaller than th * max
%   eigenvalue will be discarded. Default value is 1e-6.
%   r: maximum number of eigenvalues. Only r largest eigenvalues, if still
%   available given the threshold, are kept.

    if ~exist('r', 'var')
        r = size(X, 2);
    end
    if ~exist('th', 'var')
        th = 1e-6;
    end
    
    largest = eigs(S, 1);
    
    [v, s] = eigs(S, r);
    s = diag(s);
    mask = s / largest > th;
    s = s(mask);
    v = v(:, mask);
    
    D = 0;
    temp = X * v;
    for ii = 1:length(s)
        D = D + (temp(:, ii) - temp(:, ii)') .^ 2 / s(ii);
    end
    D = real(sqrt(D - min(min(D))));
end

