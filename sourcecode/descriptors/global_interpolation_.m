%   This is the function of Global interpolation through all the reflectance spectra of the given
%   hyperspectral image (denoted by R in the input parameter list). This function was adapted from
%   Algorithm A9.1 on page 369-370 of "The NURBS book". 
%   
%   Function [KNOTS, EVA_PT, CP_REFLEC, CP_WAVE] = global_interpolation_(R, WAVELENGTHS, deg, debug)
%
%   Input:
%       R:           the hyperspectral reflectance (or radiance) image (height x width x n). 
%       WAVELENGTHS: the wavelengths at which R was captured.
%       deg:         the degree of the B-spline basis function (by default: 3).
%
%   Output:
%       KNOTS:     knot vector (a column vector of size m x 1), same for all pixels.

%       EVA_PT:    parametric evaluation points corresponding to the bands (same for all pixels).
%   
%       CP_REFLEC: coordinates of control points in reflectance dimension, where CP_REFLEC(h, w, :)
%                  are the control point coordinates across all the bands at pixel coordinate (h,
%                  w). It is of the size: height x width x n 
%
%       CP_WAVE:   coordinates of control points in wavelength dimension (this is a 1 x n row
%                  vector). These coordinates in this dimension are the same for all the pixels in
%                  the given image.  
%   
%   Notes: 
%
%       All the reflectance spectrum (bands) correspond to the same set of bands indepedent
%       parameters (in the parameter domain). These parameters is computed by averaging the
%       independent parameter sets resulting from each individual spectrum. In this implementation,
%       the parameters are computed from the original distribution of the sampled
%       reflectance-wavelength pairs in each spectrum using the centripetal method of (Lee89)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2015 All Rights Reserved.
% Author: Ran Wei and Lin Gu
% Version: 1.0.4
% Last Update Date: 19 Lune 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version: 1.0.3
% Last Update Date: 22 Aug 2014

%   Version 1.0.2
%   Changes: added type check for wavelengths and R
%   Version 1.0.1
%
%   Date: 2012.11.10

%   this function is modified from Cong's CurveInterp3D.m

function [KNOTS, EVA_PT, CP_REFLEC, CP_WAVE] = global_interpolation_(R, WAVELENGTHS, deg, debug)

    if ~exist('debug', 'var') || numel(debug) ~= 1 || debug(1) < 0
        debug = 1;
    end

    if ~exist('deg', 'var') || numel(deg) ~= 1 || deg(1) < 0
        deg = 3;
    end
    
    if debug >= 2
        disp('Start to generate global interpolation');
    end
    
    %   n: the number of control points, also the number of wavelengths (band number)
    [height, width, n] = size(R);
    
    %   use deg to denote deg for convenience.     
    m = n + deg + 1;

    if debug >= 3
        s = sprintf('Input reflectance matrix size: height: %d, width: %d, bands: %d, Degree: %d, m = %d', height, width, n, deg, m);
        disp(s);
    end
    
    %   Compute the parameters passed to the B-spline functions. One parameter per control point.
    %   Use the centripetal method. by construction, d is float type.
    %   R is float naturally
    d = zeros(height, width);
    
    if isinteger(WAVELENGTHS)
        WAVELENGTHS = double(WAVELENGTHS);
    end
    
    for i = 1:n-1
        % centripetal method of initialising the parametric points.
        d = d + ((R(:, :, i+1) - R(:, :, i)).^2 + (WAVELENGTHS(i+1) - WAVELENGTHS(i)).^2).^(1/4);      
    end

    %   The independent parameters for each separate spectrum.
    uk = zeros(height, width, n);
    uk(:, :, 1) = 0;
    uk(:, :, n) = 1;
    for i = 2:n-1
        %   centipetal method of initialising the parametric points.
        uk(:, :, i) = uk(:, :, i-1) + ((R(:, :, i+1) - R(:, :, i)).* (R(:, :, i+1) - R(:, :, i)) + ...
                            (WAVELENGTHS(i+1) - WAVELENGTHS(i)).*(WAVELENGTHS(i+1) - WAVELENGTHS(i))).^(1/4)./d;        
    end
    
    %   Compute the shared parameters by averaging across the spatial domain
    EVA_PT = reshape(mean(mean(uk, 1), 2), [1 n]);
    clear uk;
    
    %   Compute the knots (shared for all the pixels)
    KNOTS = zeros(1, m);
    KNOTS(1:deg+1) = 0;
    KNOTS(m-deg:m) = 1; % m-deg == n + 1, should be m-deg-1? which is n?
    for i = 2:n-deg
        KNOTS(i+deg) = sum(EVA_PT(i:i+deg-1))/deg;
    end

    %   Calculate the matrix containing the evaluation of the B-spline functions at evaluation
    %   points EVA_PT(:). This is the coefficient matrix of the linear system with the control point
    %   coordinates as the variables.     
    A = zeros(n, n);
    for i = 1:n
        span = findspan_(n-1, deg, EVA_PT(i), KNOTS);
        %   This returns a vector with deg + 1 elements
        A(i, span-deg : span) = get_basis_(span, EVA_PT(i), deg, KNOTS);
    end

    %   Control-point coordinates in the wavelength dimension.
    %  A(isnan(A)) = 0;
    CP_WAVE = A\reshape(WAVELENGTHS,length(WAVELENGTHS),1);    
    
    %   Init the control point coordinates in the reflectance dimension.
    CP_REFLEC = zeros(height, width, n);
    
    %   Solve for the control point coordinates in the reflectance dimension. Perform the operation
    %   per blocks of block_height x block_width pixels. to prevent out of memory error.
    block_height = 200; 
    block_width  = 200; 
    if debug >= 2
        disp('Processing data...');
    end
    for i = 1:block_height:height
        for j = 1:block_width:width
            
            block_rowmax    = min(i+block_height-1, height);
            block_colmax    = min(j+block_width-1, width);                        
            BLOCK_R         = R(i:block_rowmax, j:block_colmax, :);
            [cur_block_height, cur_block_width, ~] = size(BLOCK_R); 
            blockpixels   = cur_block_height * cur_block_width;

            BLOCKR_2D = reshape(BLOCK_R, [blockpixels, n])';
            sol = A\BLOCKR_2D;
            
            for b = 1:n
                CP_REFLEC(i:block_rowmax, j:block_colmax, b) ...
                    = reshape(sol(b, :), [cur_block_height, cur_block_width]);
            end
            
        end
    end
    
    if debug >= 2
        disp('finish generating global interpolation');
    end
end