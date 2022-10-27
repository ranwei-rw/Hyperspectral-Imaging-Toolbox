%   function POINTS = find_parapoints_(WAVES, error, degree, KNOTS, CP)
%
%   Find the independent parametric evaluation point at which the evaluation of a B-Spline curve
%   with the given degree (degree), knots (KNOTS) and control points (CP) does not differ from the
%   given vector of coordinates WAVES more than a threshold (error). The prior assumption is that x is a
%   monotonic function with respect to the parameter t.
%
%   Output:
%       POINTS:  the approximating parameter values that yield the given coordinates WAVES, stored
%                as a row vector. 
%
%   Input:
%   
%       WAVES:   the wavelength (a m x 1 column vector).
%       error:   error tolerance of the approximating values in X (abs(WAVES - X(i)) < error).
%       degree:  degree of the curve.
%       KNOTS:   the knot vector (in the form of a row (1 x m) vector).
%       CP:      control point coordinates (n_ x r) of the curve in the dimension of given WAVES.  
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 - 2015 All Rights Reserved.
% Author: Ran Wei and Cong Phuoc Huynh
% Version: 1.0.7
% Last Update Date: 10 June 2015

% Version: 1.0.6
% Last Update Date: 21 Aug 2014

%   Changes:
%   Version: 1.0.3
%   Add type check for wavelengths to make data compatible
%   Version: 1.0.2
%   Date: 6 Feb 2014

function POINTS = find_parapoints_(WAVES, error, degree, KNOTS, CP)

    bands       = length(WAVES);
    POINTS      = zeros(1, bands);
    LOWER_BOUND = zeros(1, bands);
    UPPER_BOUND = ones(1, bands);
    %WAVES = reshape(WAVES, bands, 1);
    
    %   make sure wavelength is double type
    if ~isfloat(WAVES)
        WAVES = double(WAVES);
    end
    
    [h, w] = size(WAVES);
    
    if h == 1 && w~= 1
        WAVES = WAVES';
    end;
    %   set up a counter
    iteration  = 1;
    %   Use binary search to search for the parameter values that map to WAVES.
    while true
        %   If the given WAVES does not lie within the interval [LOWER_BOUND, UPPER_BOUND], then
        %   increase the upper bound or decrease the lower bound 
        MID_BOUND = (LOWER_BOUND + UPPER_BOUND)/2;
        %   is 1% faster than the below option. 
        LOW_X     = univar_bspline_(degree, CP, KNOTS, LOWER_BOUND);
        HIGH_X    = univar_bspline_(degree, CP, KNOTS, UPPER_BOUND);
        MID_X     = univar_bspline_(degree, CP, KNOTS, MID_BOUND);
         
%         T = univar_bspline_(degree, CP, KNOTS, [LOWER_BOUND;UPPER_BOUND;MID_BOUND]);
%         LOW_X  = T(:, 1);
%         HIGH_X = T(:, 2);
%         MID_X  = T(:, 3);
         
        %   Clamp out-of-range values to [0, 1]
        LOGI_LOW_X  = logical(WAVES < LOW_X - error);
        LOGI_HIGH_X = logical(WAVES > HIGH_X + error);        
        
        POINTS(LOGI_LOW_X)  = 0;
        POINTS(LOGI_HIGH_X) = 1;                
        
        %   evaluate the x values at the current mid-point
        
        
        DIFF_WAVE = abs(MID_X - WAVES);
        LOGI_DIFF = ~LOGI_LOW_X & ~LOGI_HIGH_X;
        if iteration <= 50
            if max(DIFF_WAVE(LOGI_DIFF), [], 1) <= error
                POINTS(LOGI_DIFF) = MID_BOUND(LOGI_DIFF);
                break;
            else
                iteration = iteration + 1;
            end
        else
            %   when iteration reaches a certain amount and it hasn't converged, we increase the threshold by double.
            error = error * 2;
            iteration = 1;
        end

        %   If WAVES lies within the lower and upper bound values, then we compare its value with
        %   MID_X. A mask for which all MID_X with lower values than the required WAVES are
        %   assigned 1. 
        MASK = logical(MID_X < WAVES);
        
        %   Tighten the lower and upper bounds on POINTS
        LOWER_BOUND(MASK & ~LOGI_HIGH_X) = MID_BOUND(MASK & ~LOGI_HIGH_X);
        UPPER_BOUND(~MASK & ~LOGI_LOW_X) = MID_BOUND(~MASK & ~LOGI_LOW_X);        
    end
end
