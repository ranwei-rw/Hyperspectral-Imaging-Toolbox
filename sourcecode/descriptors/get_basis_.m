%   Description
%                       N = get_basis_(span, point, deg, KNOTS)
%   Compute the non-zero basis function vectors at a parametric point point. 
% 
%   INPUT:
%
%       span:   knot span ( from function findspan() )
%       point:  parametric point
%       deg:    degree of the spline basis functions. Note: span should be
%               equal to or greater than deg. 
%       KNOTS:  knot sequence (row vector containing [1 x (m + 1)] knots where (m = n + deg + 1), n + 1 is the
%               number of control points).  
%
%   OUTPUT:
%
%       N:      Non-zero basis functions vector with degree deg, 
%               evaluated at point, stored in an array of size [1 x (deg+1)].
%
%   Explanation:
%
%       According to Property P2.1 (page 55 of "The NURBS book" 2nd edition), N_{span, deg}(point) = 0 if point
%       is outside the interval [u_{span}, u_{span+deg+1}) `.e. assuming that the knot span of point is 
%               s (u_{s} <= point < u_{s + 1}), 
%       then 
%               N_{span, deg}(point) = 0 if s < span or s >= span + deg + 1, 
%       So the algorithm below aims to compute N_{span, deg}(point) with span such that 
%               s - deg <= span <= s. 
%       The function is coded according to Algorithm A2.2 from 'The NURBS BOOK' (Page 70).
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2014 All Rights Reserved.
% Author: Ran Wei and Cong Phuoc Huynh
% Version 1.0.2
% Date: 16 Sep 2014
%
%   this function is modified from BasisFuncs.m

function N = get_basis_(span, point, deg, KNOTS)

    %   check whether KNOTS is a vector
    [h, w] = size(KNOTS);

    if h ~= 1 && w ~= 1
        error('Please use a row vector in KNOTS');
    end

    %   declare variables 
    LEFT  = zeros(1, deg+1);
    RIGHT = zeros(1, deg+1);
    N     = zeros(1, deg+1);
    N(1)  = 1.0; % In Matlab, vectors are indexed from 1.

    for j = 2:deg+1 % j is the iterating variable through the degree of the basis functions.
        LEFT(j)  = point - KNOTS(span+2-j);
        RIGHT(j) = KNOTS(span-1+j) - point;
        saved    = 0.0;

        for r = 1:j-1 % r runs through the index of the basis function. 
            temp  = N(r) / (RIGHT(r+1) + LEFT(j-r+1));
            N(r)  = saved + RIGHT(r+1) * temp;
            saved = LEFT(j-r+1) * temp;
        end

        N(j) = saved;
    end

    %   so far, haven't found the following case happens
%     %   check for nan
%     if sum(isnan(N))
%         N = zeros(size(N));
%         N(1) = 1;
%     end

%   end of function
end
