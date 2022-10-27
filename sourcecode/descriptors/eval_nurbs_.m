%   Function [R, WAVE_NEW] = eval_nurbs(KNOTS, WAVES, CP_REF, CP_WAVE, degree)
%
%   This function is designed to reconstruct the hyperspectral image using NURBS matrix 
%
%   Input:
%       KNOTS:   the final minimal knot vector - row column.
%       
%       CP_REF:  the minimal control point coordinates in the reflectance dimension - in the form of a 3D array of size
%                (height x width x ctrl_pts). 
% 
%       CP_WAVE: the minimal control point coordinates in the wavelength dimension - in the form of a row vector (1 x
%                ctrl_pts). 
%
%       degree:  the degree of the basis (polynominal) functions. by default degree = 3;
%
%       WAVES:   the wavelengths at the bands of the image.
% 
%   
%   Output:
%
%       R:         the input multispectral image, stored as a 3D array of size (height x width x band).
% 
%       WAVE_NEW:  the actual wavelength where the generated image is built on 
% 
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013-2014 All Rights Reserved.
%   Version 1.0.3
%   Date: 2 Sep 2014


%   Version 1.0.2
%   Date: 2013.1.30
%
%   This function is developed by Ran Wei, Cong Phuoc Huynh and Antonio Robles Kelly.


%   This algorithm was copied from Cong's TestCurveInterp3D.m. The idea is: firstly, we need to get estimations of
%   parameter t from given knots and wavelength, wavelength control points. After that, using t and other paramaters to
%   reconstruct the hyperspectral image
%

% Estimate the parametric point corresponding to the given bands.

function [R, WAVE_NEW] = eval_nurbs_(KNOTS, WAVES, CP_REF, CP_WAVE, degree)

    if ~exist('degree', 'var')
        degree = 3;
    end

    if size(WAVES, 2) > size(WAVES, 1)
        WAVES = reshape(WAVES, [length(WAVES) 1]);
    end

    POINTS        = find_parapoints_(WAVES, 0.00001, degree, KNOTS, CP_WAVE);
    [R, WAVE_NEW] = reconstruct_curve_(KNOTS, POINTS, CP_REF, CP_WAVE, degree);
end
