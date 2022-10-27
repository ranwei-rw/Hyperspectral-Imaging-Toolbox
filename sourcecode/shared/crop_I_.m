%% Resize a hyperspectral image cube
%
%% Syntax:
%     Q = crop_I(I, rect);
% 
%% Description:
%     Crops an image cube comprised of n wavelength indexed bands.
% 
%% Input:
%     I: Image cube
%     rect: rect is a four-element position vector[xmin ymin width height] 
%           that specifies the size and position of the crop rectangle. 
% 
%% Output:
%     Q: Cropped image cube.
%
%% Example
%   
%     Crop an image cube based on a rect of width 150 and height of 200 pixels.
%
%   I = FLAread('.\shared\samples\face.fla');
%
%   Q = crop_I(I.I, [80 70 150 200]);
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.2
% Last Update Date: 12 Aug 2014

function Q = crop_I_(I, rect)

    if ndims(I)==3
        [~, ~, bands] = size(I);    
        for i = 1:bands
             Q(:, :, i) = double(I(round(rect(2)):round(rect(2)+rect(4)-1), ...
                                   round(rect(1)):round(rect(1)+rect(3)-1), i));
        end
    else
        Q(:, :) = double(I(round(rect(2)):round(rect(2)+rect(4)-1), ...
                           round(rect(1)):round(rect(1)+rect(3)-1)));
    end

end
