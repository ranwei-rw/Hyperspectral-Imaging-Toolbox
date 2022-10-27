% White patch method for estimating the illuminant colour from an image 
% by taking the maximal image brightness in each channel 
% as the illuminant spectral band value.
%
% Usage:
%   L = light_white_patch(I, drate)
%
% Input: 
%   I:     a spectral (or RGB) image cube with size height x width x band
%   drate: the downsampling rate.
% 
% Output: 
%   L:     the illumination power spectrum with size band x 1
% 
%% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly.
% Version: 1.0.1
% Last Update Date: 15 Aug 2014 

function L = light_white_patch_(I, drate)

    [~, ~, band] = size(I);
    
    L = zeros(band, 1);

    % Estimate L by taking the brightest value in each band
    for n = 1:band
        temp    = I(:, :, n);
        downI   = temp(1:drate:end, 1:drate:end);
        L(n, 1) = max(downI(:));
    end
    
    clear temp;
end