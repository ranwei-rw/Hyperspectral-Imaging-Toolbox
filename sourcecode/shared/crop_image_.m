%% Crop an HSZ or hyperspectral image
%
%% Syntax:
%     HSZ = crop_image(HS, rect);
%     I = crop_image(Im, rect);
%
%% Description:
%     Crops an HSZ or hyperspectral image.
% 
%% Input:
%     I: Image data structure
%     HS: Scyllarus hyperspectral data structure
%     rect: rect is a four-element position vector[xmin ymin width height] 
%           that specifies the size and position of the crop rectangle. 
% 
%% Output:
%     I: Cropped image data structure.
%     HSZ: Cropped Scyllarus data structure.
%
%% Example
%   
%     Crop an image based on a rect of width 150 and height of 200 pixels.
%
%   I = FLAread('.\shared\samples\face.fla');
%   Q = crop_image(I.I, [80 70 150 200]);
%
%% See also
%
%     crop_I_, resize_image, crop_I_
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.1
% Last Update Date: 15 Jan 2014

function Q = crop_image_(I, rect)

% Setup the variables to start with
Q.HDR = I.HDR;

%crop an image
if isfield(I, 'I')
    Q.I = crop_I_(I.I, rect);
    [Q.HDR.lines, Q.HDR.samples, ~] = size(Q.I);
else
    %Do the HSZ
    %
    %Start with the reflectance
    if isfield(I, 'S')
        Q.S = I.S;
        Q.S.Factor = crop_I_(I.S.Factor, rect);
        if I.HDR.IndexedS==0
            Q.S.Elements = crop_I_(I.S.Elements, rect);
            [Q.HDR.lines, Q.HDR.samples, ~] = size(Q.S.Elements);
        else
            Q.S.ElementAbundanceIndexes = crop_I_(I.S.ElementAbundanceIndexes, rect);
            Q.S.ElementAbundances = crop_I_(I.S.ElementAbundances, rect); 
            [Q.HDR.lines, Q.HDR.samples, ~] = size(Q.S.ElementAbundances);
        end
        if I.HDR.EndmemberIndexedS==1
            Q.S.EndmemberAbundanceIndexes = crop_I_(I.S.EndmemberAbundanceIndexes, rect);
            Q.S.EndmemberAbundances = crop_I_(I.S.EndmemberAbundances, rect);    
        end
    end
    %Do the illuminant
    if isfield(I, 'L')
        Q.L = I.L;
        if I.HDR.IndexedL==0
            Q.L.Elements = crop_I_(I.L.Elements, rect);
        else
            Q.L.ElementAbundanceIndexes = crop_I_(I.L.ElementAbundanceIndexes, rect);
            Q.L.ElementAbundances = crop_I_(I.L.ElementAbundances, rect);    
        end
        if I.HDR.EndmemberIndexedL==1
            Q.L.EndmemberAbundanceIndexes = crop_I_(I.L.EndmemberAbundanceIndexes, rect);
            Q.L.EndmemberAbundances = crop_I_(I.L.EndmemberAbundances, rect);    
        end
    end
    %Do the specularities
    if isfield(I, 'K')
        Q.K = I.K;
        Q.K.Factor = crop_I_(I.K.Factor, rect);
        if I.HDR.IndexedK==0
            Q.K.Elements = crop_I_(I.K.Elements, rect);
        else
            Q.K.ElementAbundanceIndexes = crop_I_(I.K.ElementAbundanceIndexes, rect);
            Q.K.ElementAbundances = crop_I_(I.K.ElementAbundances, rect);    
        end
    end
end
