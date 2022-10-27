% Shade of Gray method for estimating the illuminant colour from an image 
% by taking the L^p mean of the image brightness in each channel, where p 
% is an order taken as input.
%
% Usage:
%   L = light_shade_gray(I, drate, order)
% 
% Input: 
%   I:     a spectral (or RGB) image with size height x width x band
%   drate: the downsampling rate.
%   order: the order p of the L^p mean. By default, order = 2;
%
% Output: 
%   L: the illumination power spectrum with size band-by-1
% 
%% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly.
% Version: 1.0.2
% Last Update Date: 5 Nov 2014

%   added default value for order

% Version: 1.0.1
% Last Update Date: 15 Aug 2014

function L = light_shade_gray_(I, drate, order)

    if ~exist('order', 'var')
        order = 2;
    end
    
    [~, ~, band] = size(I);
    
    L = zeros(band, 1);

    % Estimate L by the shade of grey method
    for n = 1:band
        temp    = I(:, :, n);        
        downI   = temp(1:drate:end, 1:drate:end);
        L(n, 1) = mean(downI(:).^order).^(1/order);
    end
    
    clear temp;
    
end