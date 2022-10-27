%% Convert a trichromatic image into an LFA data structure
%
%% Syntax:
%     I = rgb2fla(IMAGE);
% 
%% Description
%   Converts an RGB image into an FLA (flat) image structure.
%
%% Input:
%     IMAGE: RGB of gray-scale image contained in a cols x rows x 3 array.
%
%% Output:
%     I: Structure containing the flat image. This is designed to be 
%        consistent with that delivered by FLAread. 
%             
%% Example:
%   IMAGE = imread('test.jpg'); 
%   I     = rgb2fla(IMAGE);
%   HSZ   = Scyllarus(I);
%
%% See also:
%   Scyllarus, FLAwrite, FLAread
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Antonio Robles-Kelly. 
% Version: 1.0.6
% Last Update Date: 26 July 2014

function I = rgb2fla_(Image)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get the image structure and check the input 
    % corresponds to color or RGB data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n m t] = size(Image);

    if t == 1
            fprintf('Converting a gray-scale image. This will yield a single-band FLA data structure\n');
            I.HDR.wavelength = 550;
            I.HDR.description = 'Gray-scale image converted into an FLA - Scyllarus Matlab Toolbox';
            I.I = double(Image);
    elseif t == 3
            I.HDR.wavelength = reshape([450,560,640], 3, 1);
            I.HDR.description = 'RGB image converted into an FLA - Scyllarus Matlab Toolbox';
            I.I = double(cat(3,Image(:, :, 3),Image(:, :, 2),Image(:, :, 1)));
    else
        error('The image is neither gray-scale or RGB. Now exiting...\n');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Finish off with the header paramerers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I.HDR.samples = m;
    I.HDR.lines = n;
    I.HDR.bands = t;
    I.HDR.interleave = 'BSQ';
    I.HDR.header_offset = 0;
    I.HDR.byteorder = 0;
    I.HDR.reflectance_scaler = 0;

end