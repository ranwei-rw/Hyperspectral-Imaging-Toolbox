% Usage:
%     X = mds(W, dims);
%     
% Description:
%   Computes the multidimensional scaling embedding of the adjacency matrix W.
%   The algorithm centers the matrix and performs an eigendecomposition
%   on the resulting matrix so as to recover the final embedding coordinates using
%   the eigenvectors weighted by their rank-ordered eigenvalues.
%     
% Inputs:
%   W: weights, a square symmetric matrix. This can correspond to a
%      binary or weighted graph 
%   dims: Number of dimensions returned after the embedding
% 
% Outputs:
%   X = Matrix of embedding coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = mds_(W, dims)

    [n, m] = size(W);
    
    if n ~= m
        disp('Input is not a square matrix. It will be downsize to a square one');
        n = min(n, m);
        W = W(1:n, 1:n);
    end
    
    if ~exist('dims', 'var')
        dims = n;
    else
        if dims > n
            dims = n;
        end
    end

    %   Compute grand and row wmeans
    column_wmean = mean(W, 2);
    wmean = mean(column_wmean, 1);
    %   Center the matrix
    WII = repmat(column_wmean, [1, n]);
    WII = WII .* WII;
    B = -(W .* W - WII - WII' + wmean^2)*0.5;

    [S, V] = eig(B);

    [V, I] = sort((diag(V)));
    X = zeros(dims, n);
    for i = 1:n
        for j = 1:dims
            X(j, i) = S(i, I(n-j+1))*sqrt(abs(V(I(n-j+1))));
        end
    end
end