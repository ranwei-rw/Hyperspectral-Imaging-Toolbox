%% Gray Edge method for estimating the illuminant colour from an image 
% by taking the Minkowski mean of the image gradient magnitudes 
% in each channel, where an order of the mean is taken as input.
%
%% Usage:
%   L = light_gray_edge(I, drate, order)
% 
%% Input: 
%   I:      a spectral (or RGB) image with size height x width x bands
%   drate:  the downsampling rate.
%   order:  the order of the Minkowski mean. Default order == 2;
% 
%% Output: 
%   L:      the illumination power spectrum with size bands x 1
% 
%% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly.
% Version: 1.0.2
% Last Update Date: 5 Nov 2014
% Added default order value

function L = light_gray_edge_(I, drate, order)

    if ~exist('order', 'var')
        order = 2;
    end
    
    [~, ~, band] = size(I);    
    L = zeros(band, 1);

    % Estimate L by the shade of grey method
    for n = 1:band
        temp     = I(:, :, n);        
        downI    = temp(1:drate:end, 1:drate:end);
        [Ix, Iy] = gradient(downI);
        L(n, 1)  = (mean(abs(Ix(:)).^order)^(1/order) + mean(abs(Iy(:)).^order)^(1/order))/2;        
    end
    
    clear downI;    
    clear temp;
    
end