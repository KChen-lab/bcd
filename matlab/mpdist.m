function [D] = mpdist(X, S)
%MPDIST(X, S) calculates Mahalanobis pairwise distance
%   X: data matrix, each row is an instance and each column is a feature
%   S: Covariance matrix
%   If you run into error because S is not strictly positive definite, try
%   gmpdist which generalizes Mahalanobis distance to positive SEMIdefinite
%   covariance matrices.
    X = X / chol(S);
    X2 = sum(X .* X, 2);
    XX = X * X';
    D = X2 + X2' - 2 * XX;
    D = real(sqrt(D - min(min(D))));
end

