%
%   Usage:
%       data = gaussianclouds(N1, N2, step);
%
%   Generates randomised data for svm training and testing which comprises
%   two intertwined semicircles in 2D space, similar to a two-arm milky way.
%
%   Inputs:
%       N1, N2: Number of data points for each of the Gaussian point clouds
%       step: separation between the clouds, it delivers data structure:
%           .X  Training vectors. Each row is a data unit.
%           .Y  Labels (1 for class 1 and 2 for class 2).
%
%   Example:
%       data = gaussianclouds(100, 100, 4);
%
%   See also 
%       annulus, milkway
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = gaussianclouds_(N1, N2, step)

    X = randn(N1, 2);
    data.X = cat(1, randn(N2, 2), X + step);
    data.Y = cat(1, ones(N1, 1), 2*ones(N2, 1));

end