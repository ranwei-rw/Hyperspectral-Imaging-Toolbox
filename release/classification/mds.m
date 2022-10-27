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

function X = mds(W, dims)
    
    if nargin == 2
        X = mds_(W, dims);
    else
        X = mds_(W);
    end
    
end