%   Function R = eval_gaussian_mixture_(MIXCOEFF, MEAN, STD, WAVELENGTH)
%   
%   This function is used to reconstruct the hyperspectral spectrum from Gaussian descriptors. These descriptors include
%   the mixing coefficients, MEAN and standard deviations of the Gaussian components fitted to each spectrum.
%
%   Input:
%
%       MIXCOEFF, MEAN, STD: Mixture coefficients, means and standard deviation (of size height x width x M)
%
%       WAVELENGTH: the wavelength vector at which the reconstructed reflectance is evaluated
%
%   Output:
%
%       R: the reconstructed image, of size (height x width x band), where band is the number of wavelengths in WAVELENGTH.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei, Cong Phuoc Huynh and Antonio Robles Kelly.
% Version: 1.0.6
% Last Update Date: 16 June 2014

%   This function was adopted from Cong's ReconstructGauss_HyperImg.m 

function R = eval_gaussian_mixture_(MIXCOEFF, MEAN, STD, WAVELENGTH)

    band = length(WAVELENGTH);
    [height, width, ~] = size(MIXCOEFF);
    n = ndims(MIXCOEFF);
    
    if n > 3
        error('Error in the size of MIXCOEFF. Could not handle matrix whose dimensions are more than 3D');
    elseif n > 2
        R = zeros(height, width, band);
    else
        R = zeros(height, band);
    end
    
    %Avoid zero mixing coefficient vetors and negative coefficients
    MIXCOEFF(MIXCOEFF < 0) = 0;
    if sum(MIXCOEFF) == 0
        MIXCOEFF = ones(size(MIXCOEFF));
    end
    %Avoid std's out of range
    STD(STD <= 0) = 0.01;
    for i = 1:band
        %   lh denotes the likelihood of band i being generated from each component. lh is computed for all the pixels in
        %   the current band and for all components.
        lh = normpdf(WAVELENGTH(i), MEAN, STD);
        if n > 2
            R(:, :, i) = sum(MIXCOEFF .* lh, 3); % sum the joint pdf across components
        else
            R(:, i) = sum(MIXCOEFF .* lh, 2); % sum the joint pdf across components
        end
    end

%   function ends
end
