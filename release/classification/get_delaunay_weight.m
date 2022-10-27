% Usage
%
%   W = get_delaunay_weight_(TRI, POINTS, h);
%
% Description:
%
%   Computes the weight matrix for the Delaunay triangulation with
%   simplexes TRI and node coordinates POINTS. It uses the constant h for
%   the graph weights computed using an exponential decay over the distance
%   for adjacent nodes.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = get_delaunay_weight(TRI, POINTS, h)

    W = get_delaunay_weight_(TRI, POINTS, h);

end