%   Usage:
%
%       data = gaussianclusters(clusters, n, var, spread)
%
%   Generates data for training and testing which comprises a number of
%   Gaussian clouds 
%
%   Inputs:
%       clusters: Number of Gaussian clusters
%       n:        Mean number of cluster elements about the cluster set
%       var:      Variance for the number of cluster elements across the cluster set
%       spread:   Mean spread (distance) between cluster centers
%
%   Outputs:
%       data [struct] Binary labeled training data points:
%           .X [2 x (Ncenter x Nperiphery x Nannulus)] Training vectors.
%           .Y [1 x (Ncenter x Nperiphery x Nannulus)] Labels (1 or 2).
%
%   Example:
%       data = gaussianclusters(4, 100, 25, 4);
%
% See also 
%  annulus, milkyway, gaussianclouds  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = gaussianclusters_(clusters, n, var, spread)


    %Compute the data points
    TOKENS = n + floor(randn(1, clusters)*var);
    CENTERS = randn(2, clusters)*spread;
    X = ones(TOKENS(1), 2)*diag(CENTERS(:, 1)) + randn(2, TOKENS(1))';
    
    for i = 2:clusters
        X = cat(1, X, ones(TOKENS(i), 2)*diag(CENTERS(:, i)) + randn(TOKENS(i), 2));
    end
    
    data.X = X;
    %Compute the labels
    Y = ones(TOKENS(1), 1);
    for i = 2:clusters
        Y = cat(1, Y, i*ones(TOKENS(i), 1));
    end
    data.Y = Y;

end
