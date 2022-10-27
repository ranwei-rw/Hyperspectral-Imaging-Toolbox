%% Resize an HSZ or hyperspectral image
%
%% Syntax:
%     HSZ = resize_image(HS, rows, cols);
%     HSZ = resize_image(HS, scale);
%     I   = resize_image(Im, rows, cols);
%     I   = resize_image(Im, scale);
%
%% Description:
%     Resizes an HSZ or hyperspectral image.
% 
%% Input:
%     I:            Image data structure
%     HS:           Scyllarus hyperspectral data structure
%     rows, cols:   New image cube dimensions
%     scale:        Scale up to which the image is to be resized.
% 
%% Output:
%     I: Resized image data structure.
%     HSZ: Resized Scyllarus data structure.
%
%% See also
%
%     crop_I, crop_image, resize_I_
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.1
% Last Update Date: 15 Jan 2014

function Q = resize_image_(I, rows, cols)

% Setup the variables to start with
Q.HDR = I.HDR;

if exist('cols', 'var')
    param = [rows cols];
else
    param = rows;
end

%Resize an image
if isfield(I, 'I')
    Q.I = resize_I_(I.I, param);
    [Q.HDR.lines, Q.HDR.samples, ~] = size(Q.I);
else
    %Do the HSZ
    %
    %Start with the reflectance
    if isfield(I, 'S')
        Q.S = I.S;
        Q.S.Factor = resize_I_(I.S.Factor, param);
        if I.HDR.IndexedS == 0
            Q.S.Elements = resize_I_(I.S.Elements, param);
            [Q.HDR.lines, Q.HDR.samples, ~] = size(Q.S.Elements);
        else
            Q.S.ElementAbundanceIndexes = floor(resize_I_(I.S.ElementAbundanceIndexes, param));
            Q.S.ElementAbundances = resize_I_(I.S.ElementAbundances, param); 
            [Q.HDR.lines, Q.HDR.samples, ~] = size(Q.S.ElementAbundances);
        end
        if I.HDR.EndmemberIndexedS == 1
            Q.S.EndmemberAbundanceIndexes = floor(resize_I_(I.S.EndmemberAbundanceIndexes, param));
            Q.S.EndmemberAbundances = resize_I_(I.S.EndmemberAbundances, param);    
        end
    end
    %Do the illuminant
    if isfield(I, 'L')
        Q.L = I.L;
        if I.HDR.IndexedL == 0
            Q.L.Elements = resize_I_(I.L.Elements, param);
        else
            Q.L.ElementAbundanceIndexes = floor(resize_I_(I.L.ElementAbundanceIndexes, param));
            Q.L.ElementAbundances = resize_I_(I.L.ElementAbundances, param);    
        end
        if I.HDR.EndmemberIndexedL == 1
            Q.L.EndmemberAbundanceIndexes = floor(resize_I_(I.L.EndmemberAbundanceIndexes, param));
            Q.L.EndmemberAbundances = resize_I_(I.L.EndmemberAbundances, param);    
        end
    end
    %Do the specularities
    if isfield(I, 'K')
        Q.K = I.K;
        Q.K.Factor = resize_I_(I.K.Factor, param);
        if I.HDR.IndexedK == 0
            Q.K.Elements = resize_I_(I.K.Elements, param);
        else
            Q.K.ElementAbundanceIndexes = floor(resize_I_(I.K.ElementAbundanceIndexes, param));
            Q.K.ElementAbundances = resize_I_(I.K.ElementAbundances, param);    
        end
    end
end
