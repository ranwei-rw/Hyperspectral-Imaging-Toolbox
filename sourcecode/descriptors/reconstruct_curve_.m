
% Description: 
%       [R, WAVELENGTH] = reconstruct_curve_(KNOTS, T, REF_CP, WAVES_CP, deg)
% 
%   This function is designed to recontruct the original hyperspectral image or the original set of
%   spectra from a common knot vector (knots) for all the spectra, a common set of independent
%   parameter values (T) corresponding to the sampled reflectance-wavelength pairs in all the
%   spectra, and the given control point coordinates in the reflectance dimension(REF_CP) and the
%   wavelength dimension (WAVES_CP). (Based on Algorithm 3.1 in the NURBS book).
% 
% Input:
%
%       KNOTS      : Knot sequence, row vector of size nKnots.
%
%       T          : Independent parametric evaluation points, stored in a row vector of size band, where 
%                    band is the number of bands to be reconstructed.
% 
%       REF_CP     : control point coordinates in the reflectance dimension, a 3D array of size (height x
%                    width x cp_num), where cp_num is the number of control points in the refletance
%                    (radiance) spectrum at each pixel. 
% 
%       WAVES_CP   : control point coordinates in the wavelength dimension, a column vector of size (band x 1).
%
%       deg        : Degree of the B-Spline.
%   
% Output:
%
%       R          : the reconstructed image, of size (height x width x band)
%
%       WAVELENGTH : the wavelengths corresponding to the bands in the image.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013-2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.7
% Last Update Date: 17 Sep 2015

% Version: 1.0.6
% Last Update Date: 4 June 2014




%   this function was modified from cong's ReconstructNURBS_HyperImg.m 
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Auther: Ran Wei and Cong Phuoc Huynh
% Version 1.0.1
% Date: 2013.1.14
%
function [R, WAVELENGTH] = reconstruct_curve_(KNOTS, T, REF_CP, WAVES_CP, deg)

    %   Note that the bands correspond to the parametric evaluation points.
    bands = length(T); % bands: number of bands
    
    [height, width, tn] = size(REF_CP);
    
    if tn == 1 && width ~= 1
        %   reconstruct REF_CP
        REF_CP = reshape(REF_CP, [height, 1, width]);
        width = 1;
    end
    [cp_num, cp_width] = size(WAVES_CP); % number of control points.
    if cp_num == 1 && cp_width > 1
        % we have a row vector here.
        WAVES_CP = WAVES_CP';
        cp_num = cp_width;
    end
    
    %   if there is something wrong and cp_num is still 1 or less
    %   then we'll have a problem to use findspan_ cause in there will be 
    %   something more to deal with
    
    % Calculate the matrix containing the evaluation 
    % of the B-spline functions at the arguments T(:)
    % This is the coefficient matrix of the linear system with the control
    % point coordinates as the variables.    
    A = zeros(bands, cp_num);
    for i = 1:bands
        span = findspan_(cp_num-1, deg, T(i), KNOTS);
        % This returns a (deg + 1)- element vector
        A(i, span-deg : span) = get_basis_(span, T(i), deg, KNOTS);
    end

    % Evaluate the wavelengths.
    WAVELENGTH = A * WAVES_CP;
    
    % Evaluate the spectral reflectance (or radiance)
    % Perform the operation per blocks of patch_height x patch_width pixels.
    % to prevent out of memory error.
    patch_height = 200; patch_width = 200; 
    R = zeros(height, width, bands);
    
    for i = 1:patch_height:height
        for j = 1:patch_width:width            
            block_row_max  = min(i+patch_height-1, height);
            block_col_max  = min(j+patch_width-1, width);
                        
            BLOCK_REF_CP    = REF_CP(i:block_row_max, j:block_col_max, :);
            [current_block_height, current_block_width, ~] = size(BLOCK_REF_CP); 
            block_size     = current_block_height * current_block_width;
            BLOCK_REF_CP2D  = zeros(cp_num, block_size); 
            for b = 1:cp_num
                BLOCK_REF_CP2D(b, :) = reshape(BLOCK_REF_CP(:, :, b), [block_size 1]);
            end

            % recover the reflectance (or radiance)
            R2D = A * BLOCK_REF_CP2D;            
            for b = 1:bands
                R(i:block_row_max, j:block_col_max, b) ...
                    = reshape(R2D(b, :), [current_block_height, current_block_width]);
            end
        end
    end    
end