function [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(I, WAVELENGTH, M, DEBUG)

%   Function [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(I, WAVELENGTH, M, DEBUG)
%
%   This function generates Gaussian features for a multispectral reflectance image (as the normalisation result of the
%   given radiance image by the given light).
%
%   The method used in this function is the EM algorithm for fitting. In each iteration it re-estimate the components of
%   the spectra at all pixels together. The features are the mixture coefficients, the mean and standard deviation of
%   the Gaussian components in the wavelength domain.
%
%   Input:
%
%       I:          the 3D data matrix. (height x width x band)
%       WAVELENGTH: the wavelengths at which the image was captured.
%       M:          The number of Gaussian components.
%       DEBUG:      the level of debugging information to be displayed. Default to 1. Max is 3
%
%   Output:
%
%       MIXCOEFF:  mixture coefficients of the Gaussian components (height x width x M).
%       MEAN:      the means of the Gaussian components (height x width x M).
%       STD:       the standard deviations of the Gaussian components (height x width x M).
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Version 1.0.1
% Author: Ran Wei, Cong Phuoc Huynh and Antonio Robles Kelly.
% Date: 2013.1.31

%   this function is modified from Cong's GaussFeatureGenEM_HyperImg.m

if ~exist('DEBUG', 'var')
    DEBUG = 1;
end

[height, width, band] = size(I);
%I=I/max(max(max(I)));
%   normalise the reflectance spectrum at each pixel so that its integration over the wavelength range of the image is 1.
sumR    = sum(I, 3); % sum of reflectance across all the bands.
n_input = I;% ./ sumR(:, :, ones(band, 1));


%   Initialise the Gaussian components. 
%   It is assumed that the means of the first and last Gaussian components are the first and last band. Variable
%   gaussian_cover_range determines the number of bands covered by each Gaussian component (on average) since M-1
%   components are already enough to cover the whole spectral range. 
MIXCOEFF = zeros(height, width, M);
MEAN     = zeros(height, width, M);
STD      = zeros(height, width, M);
gaussian_cover_range = floor(band / M);

for i=1:M
    MIXCOEFF(:, :, i) = 1/M;                                                                   %   average weights
    MEAN(:, :, i)     = mean(WAVELENGTH((i-1)*gaussian_cover_range+1:i*gaussian_cover_range)); %   mean wavelength
    STD(:, :, i)      = std(WAVELENGTH((i-1)*gaussian_cover_range+1:i*gaussian_cover_range));  %   standard deviation.
end

%   Compute the Gaussian feature by fitting a mixture of Gaussian components (in the probabilistic sense) to the
%   reflectance spectrum at each pixel. 

%   p is a 4D matrix whose element p(h, w, b, j) is the expectation that the spectral band i at pixel location (h, w) is
%   associated with component j. 
p = zeros(height, width, band, M);

%   The average change in the mixing coefficients, MEAN and STD across pixels and Gaussian components (in the previous
%   iteration) 
previous_avg_mix_diff  = Inf;
previous_avg_mean_diff = Inf;
previous_avg_std_diff  = Inf;

iter = 0;
max_iteration = 20; % maximum allowable number of iterations.

while (true)
    iter = iter + 1;
    if (iter >= max_iteration)
        break;
        if DEBUG >= 2
            fprintf('Maximum iteration number reached./n');
        end
    end
    
    last_mixcoef = MIXCOEFF;
    last_mean    = MEAN;
    last_std     = STD;
    
    if DEBUG >= 2
        fprintf('Iteration %d\n', iter);
    end
    %   E-Step
    for i = 1:band 
        %   possibility is likelihood of band i being generated from each component. It is computed for all the pixels.
        possibility = normpdf(WAVELENGTH(i), MEAN, STD);
        
        %   joint distribution of P(WAVELENGTH{i}, component(j)) weighted by % the reflectance values.
        
        %   weighted_joint_pdf is a temporary variable storing the joint pdf of % the reflectance value at a pixel and a
        %   band and a particular % Gaussian component.
        weighted_joint_pdf = MIXCOEFF .* possibility;
        sum_joint_pdf      = sum(weighted_joint_pdf, 3);
        p(:, :, i, :)      = weighted_joint_pdf./sum_joint_pdf(:, :, ones(M, 1));
    end
    
    %   M-step
    for j = 1:M
        MIXCOEFF(:, :, j) = sum(n_input .* p(:, :, :, j), 3);
        
        %   compute the variance then take the square root. 
        %   Implementation in this way to save memory, but it takes longer time.
        band_vari_sum     = zeros(height, width);
        band_wrefwave_sum = zeros(height, width);
        band_wref_sum     = zeros(height, width);
        for b = 1:band
            temp = n_input(:, :, b).* p(:, :, b, j) .* WAVELENGTH(b);
            band_wrefwave_sum = band_wrefwave_sum + temp;
            
            temp = n_input(:, :, b) .* p(:, :, b, j) .* (WAVELENGTH(b) - last_mean(:, :, j)).^2;
            band_vari_sum = band_vari_sum + temp;
            
            temp = n_input(:, :, b) .* p(:, :, b, j);
            band_wref_sum = band_wref_sum + temp;
        end
        band_wref_sum(find(band_wref_sum==0))=1; 
        MEAN(:, :, j) = band_wrefwave_sum ./ band_wref_sum;
        STD(:, :, j) = sqrt(band_vari_sum ./ band_wref_sum);
    end
    
    %   Check if the components do not change significantly.
    %   absolute difference between the mixture coefficients in this iteration and previous iteration. 
    mix_change  = abs(last_mixcoef - MIXCOEFF); 
    %   absolute difference between the MEAN in this iteration and previous iteration.
    mean_change = abs(last_mean - MEAN); 
    %   absolute difference between the standard deviations in this iteration and previous iteration.
    std_change  = abs(last_std - STD); 
    avg_mix_change  = mean(mix_change(:));
    avg_mean_change = mean(mean_change(:));
    avg_std_change  = mean(std_change(:));
    
    if DEBUG == 3
        fprintf('Mean change in mixing coeffs %g, in MEAN %g, in STD %g\n', avg_mix_change, avg_mean_change, avg_std_change);
    end
    
    if (isnan(avg_mix_change)  || ...
        isnan(avg_mean_change) || ...
        isnan(avg_std_change)  || ...
        previous_avg_mix_diff  < avg_mix_change  || ...
        previous_avg_mean_diff < avg_mean_change || ...
        previous_avg_std_diff  < avg_std_change)
        % diverging, change to the parameters estimated in the previous iteration.
        MIXCOEFF = last_mixcoef;
        MEAN     = last_mean;
        STD      = last_std;
        break;
    end
    
    previous_avg_mix_diff  = avg_mix_change;
    previous_avg_mean_diff = avg_mean_change;
    previous_avg_std_diff  = avg_std_change;
    
end
clear p;
%Get the scale to reconstruct the spectra accurately
R = eval_gaussian_mixture_(MIXCOEFF, MEAN, STD, reshape(WAVELENGTH,band,1));
C=sum(I./R,3)/band;
MIXCOEFF = MIXCOEFF.*C(:,:,ones(M,1));

%   function ends
end