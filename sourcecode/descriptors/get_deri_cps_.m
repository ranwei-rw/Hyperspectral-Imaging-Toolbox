function [REFLEC_CPS, WAVE_CPS] =...
         get_deri_cps_(BSP_DEG, KNOTS, CP_REFLEC, CP_WAVE, MAX_DEG, MIN_CP_INDEX, MAX_CP_INDEX)

%    
%   This function Compute the control points of all the derivatives of a hyperspectral image
%   reflectance and wavelength with respect to the independent parameter, with order up to and
%   including the MAX_DEGree^{th} derivative. 
%   
% Input:
%       BSP_DEG: B-spline degree.
%       
%       KNOTS: knot sequence, vector [1 x (m + 1)] (m = n + BSP_DEG + 1), 
%           where n is the number of control points - 1.
%           
%       CP_REFLEC: (height x width x nCtrls) the control point coordinates
%           in the reflectance dimension, where CP_REFLEC(h, w, :) are the 
%           control point coordinates across all the bands 
%           at pixel coordinate (h, w).
%   
%       CP_WAVE: (nCtrls x 1) - a row vector - the control point coordinates
%           in the wavelength dimension. These coordinates in this dimension
%           is the same for all the pixels in the given image.
%           nCtrls is the number of control points.    
%   
%       MAX_DEG: the maximum degree of the derivatives.
%       
%       MIN_CP_INDEX, MAX_CP_INDEX: the lower bound and upper bound of the
%           index of the control points in the original curves 
%           which are then used to generate control points of the
%           derivatives of higher orders according to the formula in
%           Equation 3.8 (1 <= MIN_CP_INDEX <= MAX_CP_INDEX <= n)
%           (since Matlab arrays are indexed from 1)
%              
% Output:
%
%       REFLEC_CPS: the reflectance (or radiance) 
%           coordinates of the control points of 
%           the derivative curves (which are also B-Splines), 
%           with dimension 
%           (height x width x (MAX_CP_INDEX - MIN_CP_INDEX + 1) x MAX_DEGree),       
%           where REFLEC_CPS(i, j, k, l) is reflectance (or radiance)
%           coordinate of the control points of the derivative curve of
%           order l(1 <= l <= BSP_DEG), 
%           for the pixel (i, j) and the k^th wavelength.
%       
%       WAVE_CPS: : the wavelength coordinates of the control points 
%           of the derivative curves (which are also B-Splines), 
%           with dimension ((MAX_CP_INDEX - MIN_CP_INDEX + 1) x MAX_DEGree), 
%           where WAVE_CPS(k, l) is wavelength
%           coordinate of the control points of the derivative curve 
%           of order l(1 <= l <= BSP_DEG) for the k^th sampled wavelength.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Version 1.0.1
% Author: Ran Wei
% Date: 2012.11.20

%   this function is modified from DerivCtrls_HyperImage.m.

    %   get the difference between index borders
    index_diff = MAX_CP_INDEX - MIN_CP_INDEX;    
    [height, width, nCtrls] = size(CP_REFLEC);
    
    %   Initialise the control points of the zeroth-degree curve.
    temp_REFLEC_CPS = zeros(height, width, index_diff+1, MAX_DEG + 1);    
    temp_REFLEC_CPS(:, :, :, 1) = CP_REFLEC(:, :, MIN_CP_INDEX:MAX_CP_INDEX);
    
    temp_WAVE_CPS = zeros(index_diff + 1, MAX_DEG + 1);    
    temp_WAVE_CPS(:, 1) = CP_WAVE(MIN_CP_INDEX:MAX_CP_INDEX, 1);
    

    %   Note: we add 1 to the actual degree k in each iteration to account for the fact that Matlab
    %   arrays are indexed from 1 (not 0 as in the original Algorithm A3.3)
    for k = 2:MAX_DEG+1
        temp = BSP_DEG - k + 2;
        for i = 1: index_diff - k + 2
            temp_REFLEC_CPS(:, :, i, k) = temp * ...
                (temp_REFLEC_CPS(:, :, i+1, k-1)  - temp_REFLEC_CPS(:, :, i, k-1)) ...
                /(KNOTS(MIN_CP_INDEX+i+BSP_DEG) - KNOTS(MIN_CP_INDEX+i+k-2));
            
            temp_WAVE_CPS(i, k) = temp * ...
                (temp_WAVE_CPS(i+1, k-1)  - temp_WAVE_CPS(i, k-1)) ...
                /(KNOTS(MIN_CP_INDEX+i+BSP_DEG) - KNOTS(MIN_CP_INDEX+i+k-2));
        end
    end    
    
    % Only retain the control points of the actual derivative curves (of degree >= 1)
    REFLEC_CPS = temp_REFLEC_CPS(:, :, :, 2:MAX_DEG+1);    
    WAVE_CPS = temp_WAVE_CPS(:, 2:MAX_DEG+1);
    
end