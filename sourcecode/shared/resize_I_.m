%% Resize a hyperspectral image cube
%
%% Syntax:
%     Q = resize_I(I, scale);
%     Q = resize_I(I, rows, cols);
% 
%% Description:
%     Resizes an image cube comprised of n wavelength indexed bands to be of
%     dimensions rows x cols x n.
% 
%% Input:
%     I: Image cube
%     scale: Scaling factor for the image cube. This should be a value
%           greater than 0, where values between 0 and 1 imply the output image 
%           cube is smaller in size and values greater than one yield a
%           larger image cube.
%     rows, cols: New image cube dimensions
% 
%% Output:
%     Q: Resized image cube.
%
%% Example
%   
%     Resize an image cube so as to be 256 by 320 pixels.
%
%   I = FLAread('.\shared\samples\face.fla');
%   Q = resize_I(I.I, 256, 320);
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.1
% Last Update Date: 15 Jan 2014

function Q = resize_I_(I, rows, cols)

[~,~,bands]=size(I);
for i=1:bands
    if exist('rows','var')
        if ~exist('cols','var')
            Q(:,:,i)=imresize(I(:,:,i),rows,'nearest');
%            Q(:,:,i)=imresize(I(:,:,i),rows);
        else
            Q(:,:,i)=imresize(I(:,:,i),[rows cols],'nearest');
%            Q(:,:,i)=imresize(I(:,:,i),[rows cols]);
        end
    else
        error('Please provide a scaling value or a valid image size');
    end
end
