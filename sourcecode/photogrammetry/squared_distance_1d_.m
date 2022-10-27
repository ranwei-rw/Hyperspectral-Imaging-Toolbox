%%  Calculate squared_distance between two input positions of two elements in 1D vector.
%
%% Syntax
%   sd = squared_distance(p1, p2, height, width);
%
%% Description
%
%   Get squared distance between two 3D points in neuron network which is denoted by two elements in 1D vectors 
%    
%% Input:
%   p1 and p2: positions of two elements in 1D vectors. These elements are reshaped from a 3D neuron matrix.
%   height and width: 2D dimensions of of the 3D neuron matrix
%
%% Output:
%   sd: squared distance between two points in 3D neuron network. These two points are converted from given positions.
%   
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.0
% Last Update Date: 02 July 2014

function sd = squared_distance_1d_(p1, p2, height, width)

    if p1 == p2
        sd = 0;
        return;
    end
        
    z = ceil(p1/(height*width));
    r = mod(p1, height*width);
    x = ceil(r/height);
    if x == 0
        x = width;
    end
    y = mod(r, height);
    if y == 0
        y = height;
    end
    
    V1 = [y, x, z];
    z = ceil(p2/(height*width));
    r = mod(p2/height*width);
    x = ceil(r/height);
    if x == 0
        x = width;
    end
    y = mod(r, height);
    if y == 0
        y = height;
    end
    
    V2 = [y, x, z];

    sd = sum((V1 - V2).*(V1 - V2));
end
