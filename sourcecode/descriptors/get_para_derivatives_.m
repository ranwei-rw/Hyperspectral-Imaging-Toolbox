function [R_DERIVS, WAVE_DERIVS] = get_para_derivatives_(R, WAVELENGTHS, DEG, MAX_DEG)
%   function [R_DERIVS, WAVE_DERIVS] = get_para_derivatives_(R, waves, DEG, MAX_DEG)
%
%   Compute the derivatives of the reflectance or radiance and wavelengths with respect to the
%   independent parameters, using an B-Spline interpolation to the spectrum at every pixel.
% 
%   The derivatives are computed at the independent parameter values corresponding to every wavelength. 
% 
%   Input:
%
%       R:           hyperspectral reflectance (or radiance) image (height x width x bands).
%       WAVELENGTHS: the wavelengths at which the hyperspectral image R was captured.
%       DEG:         degree of the B-spline basis functions used for interpolation.
%       MAX_DEG:     maximum degree allowed for the derivatives computed.
%
%   Output:
%
%       R_DERIVS:    derivatives of the reflectance or radiance with respect to the independent
%                    parameter, Its sizes are: height x width x bands x degree where degree =
%                    min(MAX_DEG, DEG). 
%   
%       WAVE_DERIVS: the derivative of the wavelength with respect to the independent parameter,
%                    where degree = min(MAX_DEG, DEG). These derivatives are the same for all the
%                    pixels in the given image. Size: bands x degree
%   
%   Note: 
%       - the knot vectors (in the parameter domain) is common for all pixels.
%       - the control point coordinates in the wavelength dimension are the same for all pixels.
%       - all the image reflectance spectra correspond to the same set of bands indepedent
%         parameters (in the parameter domain). These parameters are computed by averaging the
%         independent parameter sets resulting from each individual spectrum. In this
%         implementation, the parameters are computed from the original distribution of the sampled 
%         reflectance-wavelength pairs in each spectrum using the centripetal method of Lee (1989).

%   this function is modified based on DerivsWrtParameter_HyperImage.m

% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
%
%   Version 1.0.1
%
%   Date: 2012.11.18

    %   get the dimensions of input data
    [height, width, bands] = size(R);

    %   The maximal degree of the resulting derivatives.
    MAX_DEG = min(MAX_DEG, DEG);
    
    %   Initialise the derivative arrays.
    R_DERIVS    = zeros(height, width, bands, MAX_DEG);
    WAVE_DERIVS = zeros(bands, MAX_DEG);

    %   params are the independent parameter values corresponding to wavelengths.
    [knots, tk, ctrl_ref, ctrl_wave] = global_interpolation_(R, WAVELENGTHS, DEG);
    
    %   number of control points per spectrum.
    cp_num = size(ctrl_ref, 3); 
  
    
    %   Now to compute the Independent parameters corresponding to wavelengths. 
    
    %   tolerance of error in wavelength while searching for the parametric point corresponding to a
    %   wavelength.  
    wave_tolerance = 1e-4; 
    
    params = find_parapoints_(WAVELENGTHS, wave_tolerance, DEG, knots, ctrl_wave);
    
    for i = 1:bands
        %   The independent parameter corresponding to the current
        %   wavelength.
        param = params(i);
        
        %   Find the span of the current parameter
        span = findspan_(cp_num-1, DEG, param, knots);
        
        %   Compute the basis functions up to degree DEG
        N = get_basis_funcs_(span, param, DEG, knots);
        
        %   Compute the control points of all the derivative curves up to and including the
        %   d_degree^{th} derivative. 
        [R_deriv_ctrls, wave_deriv_ctrls] = ...
            get_deri_cps_(DEG, knots, ctrl_ref, ctrl_wave, MAX_DEG, span - DEG, span);

        for k = 1:MAX_DEG
            temp = N(DEG - k + 1, 1:DEG-k+1)';
            
            R_temp = R_deriv_ctrls(:, :, 1:DEG-k+1, k);
            wave_temp = wave_deriv_ctrls(1:DEG-k+1, k);
            
            %   Use the Equation 3.8 in the NURBS book to compute the
            %   derivatives
            for j = 1:DEG-k+1
                R_DERIVS(:, :, i, k) = R_DERIVS(:, :, i, k) + temp(j, 1) .* R_temp(:, :, j, k);
            end
            
            WAVE_DERIVS(i, k) =  sum(temp .* wave_temp, 1);
        end
    end
end
