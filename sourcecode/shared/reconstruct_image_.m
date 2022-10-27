%% Reconstruct a FLA data structure from a HSZ structure
%
%% Syntax:
%
%     I = reconstruct_image(HSZ);
% 
%% Description:
%     Recover a hyperspectral image structure (the image cube and the header) 
%     from an HSZ structure which can be imported from HSZ files by calling 
%     function HSZread 
% 
%% Input:
%     HSZ:  HSZ structure
%     
%% Output:
%     I:    Structure containing the image cube (I.I) and the corresponding header (I.HDR).
% 
%% See also
%
%     Scyllarus
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.6
% Last Update Date: 24 July 2014


function I = reconstruct_image_(HSZ)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Evaluate the HSZ structure if necesary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HSZ = eval_HSZ_(HSZ, HSZ.HDR.wavelength);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reconstruct the reflectance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    R = reconstruct_reflectance_(HSZ);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reconstruct the illuminant
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L = reconstruct_illuminant_(HSZ);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reconstruct the specularities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K = reconstruct_specularity_(HSZ);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reconstruct the image and create the structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I.I = L.*(R+K);

    I.HDR = HSZ.HDR;
end