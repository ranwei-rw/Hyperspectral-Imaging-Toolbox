%   Usage:
%       data  =  milkyway(point_num, spread);
%
%   Generates data for svm training and testing which comprises two
%   intertwined semicircles in 2D space, similar to a two-arm milky way. 
%
%	Inputs:
%   
%       point_num: Number of data points for each semicircles
%       spread: Spread of the semicircles
%
%	Outputs:
%       data [struct] Binary labeled training data points:
%           .X [2*point_num by 2] Training vectors.
%           .Y [2*point_num by 1] Labels (1 or 2).
%
%   Example:
%       data = milkyway(150, 0.5);
%
%   See also 
%       annulus, gaussianclouds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data  =  milkyway(point_num, spread)
    
    data  =  milkyway_(point_num, spread);
    
end