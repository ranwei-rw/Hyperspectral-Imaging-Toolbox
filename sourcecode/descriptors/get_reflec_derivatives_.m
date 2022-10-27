function DERI = get_reflec_derivatives_(R, WAVELENGTHS, DEG, MAX_DEG)

%   Function: DERI = get_reflec_derivatives_(R, WAVELENGTHS, DEG, MAX_DEG)
%   This function is designed to compute the derivatives of the reflectance or radiance (given as
%   the input variable R) with respect to the wavelength, using the B-Spline interpolation to the
%   spectrum at every pixel. The derivatives are computed at every sampled wavelength of the image
%
%   Parameters
% 
%   Input:
%
%       R:           hyperspectral reflectance (or radiance) image (height x width x bands).
%
%       WAVELENGTHS: the wavelengths at which the hyperspectral image R was captured.
%
%       DEG:         the degree of the B-spline basis functions used for interpolation.
%
%       MAX_DEG:     the maximum degree of computed derivatives.
%
%
%   Output:
%
%       DERI: a 4D vector whose first three dimensions are the same as those of the input
%             hyperspectral image R, and the fourth dimension corresponds to the order of the
%             derivative. i.e. DERI(:, :, :, k) is the k^th derivative of the image at every given
%             wavelength. 
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
%
%   Version 1.0.1
%
%   Date: 2012.11.29

%   this function is modified based on DerivsWrtWavelength_HyperImage.m

    [height, width, bands] = size(R);
        
    %   The maximal degree of the resulting derivatives.
    MAX_DEG = min(MAX_DEG, DEG);
    
    DERI = zeros(height, width, bands, MAX_DEG);
    
    %   Take the derivative with respect to the independent parameter.
    [first_derives, wave_derivs] = get_para_derivatives_(R, WAVELENGTHS, DEG, 1);

    % Compute the first derivative of the reflectance (radiance) with respect to wavelengths.
    for b = 1:bands
        DERI(:, :, b, 1) = first_derives(:, :, b, 1) ./ wave_derivs(b, 1);
    end
    
    % For higher-order derivatives, we use the chain rule 
    % involving the derivative of the immediately lower order.    
    for k = 2:MAX_DEG
        % Take the derivative of the previous dR/\dlambda
        % with respect to the independent parameter
        [previous_paras, wave_derivs] = ...
            get_para_derivatives_(DERI(:, :, :, k-1), WAVELENGTHS, DEG-1, 1);
        
        % Now use the chain rule
        for b = 1:bands
            DERI(:, :, b, k) = previous_paras(:, :, b, 1) ./ wave_derivs(b, 1);
        end
    end
    
end