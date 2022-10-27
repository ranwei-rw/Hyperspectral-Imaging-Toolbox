function N = get_basis_funcs_(KNOT_SPAN, PARA_POINTS, DEG, KNOTS)

%   function N = get_basis_funcs_(KNOT_SPAN, PARA_POINTS, DEG, KNOTS)
%
%   This function is designed to compute all the non-zero basis functions of degrees from 0 to DEG,
%   evaluated at a parametric point PARA_POINTS. This algorithm is invoked by Algorithm A3.4 in 'The
%   NURBS Book' Page 99.
% 
%   INPUT:
%
%       KNOT_SPAN:   knot span  ( Use findspan.m to get this )
%       PARA_POINTS: the parametric point at which basis functions are evaluated.
%       DEG:         the maximum degree of the spline basis functions.
%   `   KNOTS:       knot sequence (vector containing [1 x (m + 1)] knots. where m = n + DEG + 1 and
%                    n + 1 is the number of control points).
%
%   OUTPUT:
%
%       N:           A matrix of size (DEG+1) x (DEG+1) storing the values of the B-Spline basis
%                    functions, where the jth row corresponds to all the non-zero basis functions of
%                    the jth-degree evaluated at the parameter value PARA_POINTS. N(j, r) is the
%                    value of the r^th basis function of degree j - 1. 
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Cong Phuoc Huynh
% Version 1.0.1
% Date: 2012.11.18

%
%   From Property P2.1 (in Page 55 of "The NURBS Book" 2nd edition), we have 
%           N_{r, DEG}(PARA_POINTS) = 0 
%   if PARA_POINTS is outside the interval [u_{r}, u_{r+DEG+1}). i.e. assuming that the knot span of
%   PARA_POINTS is KNOT_SPAN(u_{s} <= PARA_POINTS < u_{s + 1}), then N{r, DEG}(PARA_POINTS) = 0 if
%   KNOT_SPAN < r or KNOT_SPAN >= r + DEG + 1, So the algorithm below aims to compute N_{r,
%   DEG}(PARA_POINTS) with r such that KNOT_SPAN - DEG =< r <= KNOT_SPAN

%
%   This function is adopted from BasisFuncsAllDegrees.m

    %   Initialise the resulting array
    N = zeros(DEG+1, DEG+1); %   In Matlab, vectors are indexed from 1.

    %   the only non-vanishing B-spline basis function with degree 0 have a value of 1.
    N(1, 1) = 1.0; 

    %   Initialise the left and right arrays as instructed in the NURBS book (Page 69)
    left = zeros(DEG+1, 1);
    right = zeros(DEG+1, 1);

    %   left and right in the NURBS book are indexed from 0
    %   whereas the matlab implementations of these are index from 1.
    for j = 2:DEG+1 
        left(j)  = PARA_POINTS - KNOTS(1, KNOT_SPAN + 2-j);     %   +1 to compensate for the index implementation in matlab
        right(j) = KNOTS(1, KNOT_SPAN -1 + j) - PARA_POINTS;    %   -1 to compensate for the index implementation in matlab
    end

    %   j is the iterating variable through the degree of the basis functions.
    for j = 2:DEG+1 
        saved = 0.0;

        for r = 1:j-1 %   r runs through the index of the basis function. 
            temp = N(j-1, r) / (right(r+1) + left(j-r+1)); %   the index of left is incremented by 1 (in matlab)

            %   This is a basis function
            N(j, r) = saved + right(r+1) * temp;            

            saved = left(j-r+1) * temp;
        end

        N(j, j) = saved;        
    end

end

