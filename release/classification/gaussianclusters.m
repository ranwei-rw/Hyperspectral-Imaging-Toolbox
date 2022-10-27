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
%   See also 
%       annulus, milkyway, gaussianclouds  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = gaussianclusters(clusters, n, var, spread)

    data = gaussianclusters_(clusters, n, var, spread);

end
