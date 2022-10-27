%Function
%   POINTS = univar_bspline_(deg, CTR_PTS, KNOTS, PARAMS)
% 
%   This function is designed to- Evaluate a univariate B-Spline.
%    (Algorithm 3.1 in the NURBS book).
% 
% Description:
% 
%   Evaluate a univariate B-Spline. This function provides an interface to
%   a toolbox 'C' routine.
%
% Input:
% 
%       deg:     Degree of the B-Spline.
% 
%       CTR_PTS: Control Points, matrix of size n by dim. 
%
%       KNOTS:   Knot sequence, row vector of size nKnots.
% 
%       PARAMS:  Parametric evaluation points, row vector of size params_n.
% 
% Output:
%
%       POINTS:  Evaluated points, matrix of size (params_n, dim)
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2015 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.3
% Last Update Date: 10 Jun 2015

% Version: 1.0.2
% Last Update Date: 16 Sep 2014

function POINTS = univar_bspline_(deg, CTR_PTS, KNOTS, PARAMS)

    %   this function is modified from bspeval.m.

    n   = size(CTR_PTS, 1); %   number of control points (data points)
    dim = size(CTR_PTS, 2); %   dimension of each control (or data) point.

    if n == 1 && dim > 1
        %   CTR_PTS is used as a row vector, then transpose it
        CTR_PTS = CTR_PTS';
        n       = dim;
        dim     = 1;
    end

    [num, params_n, ~] = size(PARAMS);

    POINTS   = zeros(params_n, dim*num);

    %   for each parametric point PARAMS(index)
    if num > 1
        for i = 1:num
            for index = 1:params_n
                %   find the span of u(1, col)    
                s = findspan_(n-1, deg, PARAMS(i, index), KNOTS);
                B = get_basis_(s, PARAMS(i, index), deg, KNOTS); %  B(1, i) is the ith basis function of degree d.

                tmp1 = s - deg;
                for coord = 1:dim   %   coord is the coordinate index of a point.
                    POINTS(index, (i-1)*dim + coord) = B(1:1+deg)*CTR_PTS(tmp1:tmp1+deg, coord);
                end
            end
        end
    else
        for index = 1:params_n
            %   find the span of u(1, col)    
            s = findspan_(n-1, deg, PARAMS(index), KNOTS);
            B = get_basis_(s, PARAMS(index), deg, KNOTS); %  B(1, i) is the ith basis function of degree d.

            tmp1 = s - deg;
            for coord = 1:dim   %   coord is the coordinate index of a point.
                POINTS(index, coord) = B(1:1+deg)*CTR_PTS(tmp1:tmp1+deg, coord);
            end
        end
    end

end