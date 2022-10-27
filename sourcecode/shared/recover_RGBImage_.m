% Syntax
%   IMAGE = recover_RGBImage_(DATA);
%   IMAGE = recover_RGBImage_(DATA, CMF);
%   IMAGE = recover_RGBImage_(DATA, CMF, enable_endmember);
%
% Inputs:
%
%     Data:             HSZ data struct or FLA image to be processed.
%     CMF:              Colour matching function (CMF) to be used. This is a 4xbands array whose first column contains
%                       the wavelength vector. The remaining columns account for the BGR components of the matching
%                       function. If no CMF is provided, the matching function for the Nikon D70 camera is used. 
%     'use_endmembers': Option which determines whether the endmembers embedded in the file should be used for the color
%                       rendering rather than the materials. By default, the materials are used.
%     'tonemap':        Option which determines whether the image should be tonemapped for display purposes. The default
%                       output is without tonemapping
%
%
% Outputs:
%     Irgb: Pseudo-colour image.
%
% Example:
%
%   I = FLAread('.\shared\samples\face.fla', 'fullheader');
%   HSZ = Scyllarus(I);
%   LibCMF = SLZread('StilesBurch1955CMFs2deg.slz');
%   CMF = LibCMF.Endmember(1:3,:);
%   Irgb = recover_RGBImage(HSZ, [LibCMF.HDR.wavelength; CMF(3,:); CMF(2,:); CMF(1,:)]);
%     
% See also:
%
%   FLAread, HSZread, Scyllarus
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013-2015 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.9
% Last Update Date: 9 June 2015
%   increased speed.

% Version: 1.0.8
% Last Update Date: 25 July 2014

function IMAGE = recover_RGBImage_(DATA, CMF, enable_endmember)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup the variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('enable_endmember', 'var') || (enable_endmember ~= 1 && enable_endmember ~= 0)
        enable_endmember = 0;
    end
    
    if ~exist('CMF', 'var') 
        s = sprintf('Colour matching function not provided. Using spectral response for the Nikon D70 instead.');
        disp(s);
        load('NikonD70.mat');
        CMF = NikonD70';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Load the mat files containig the responses for the Nikon D70
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CMF = CMF';
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Resample the CMF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [KNOTS, CP_REF, CP_WAVE] = get_nurbs_(reshape(CMF(:, 2:4)', [3 1 length(CMF(:, 1))]),...
                                         reshape(CMF(:, 1), [length(CMF(:, 1)) 1]));
    
    [WL, IND] = wavelength_subset_(DATA.HDR.wavelength, CMF(:, 1));
    IND = IND';
    
    if isempty(WL)
        error('There are no common bands between the image and the colour matching function... exiting');
    end
    
    Q = eval_nurbs_(KNOTS, WL, CP_REF, CP_WAVE, 2);
    C = cat(2, WL', reshape(Q, [3 length(WL)])');

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Verify the data structures and remove the
    %illuminant from the DATA data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(DATA, 'S')
        DATA = eval_HSZ_(DATA, DATA.HDR.wavelength);
        %DATA.L.Elements = ones(size(DATA.L.Elements));
        if enable_endmember == 0
            if DATA.HDR.EndmemberIndexedL == 1
                DATA.HDR.EndmemberIndexedL = 0;
            end
            if DATA.HDR.EndmemberIndexedS == 1
                DATA.HDR.EndmemberIndexedS = 0;
            end
        end
        I = reconstruct_image_(DATA);
        I.I = I.I(:, :, IND);
    elseif isfield(DATA, 'I')
        I.I = DATA.I(:, :, IND);
    else
        error('This is not a valid image or DATA structure ... exiting\height');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup the rest of the variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [height, width, bands] = size(I.I);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Normalise factors for the CMF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BGR = 1./(C(:, 2:4)'*ones(length(WL), 1));
    BGR = [C(:,4)*BGR(3)/BGR(2), C(:,3), C(:,2)*BGR(1)/BGR(2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recover the pseudocolor image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IMAGE = reshape(reshape(I.I, [height*width, bands])*BGR, [height, width, 3]);
%     IMAGE = zeros(height, width, 3);
%     for i = 1:height
%         for j = 1:width
%             IMAGE(i, j, 1) = dot(reshape(I.I(i, j, :), 1, bands), C(:,4))*BGR(3)/BGR(2);
%             IMAGE(i, j, 2) = dot(reshape(I.I(i, j, :), 1, bands), C(:,3));
%             IMAGE(i, j, 3) = dot(reshape(I.I(i, j, :), 1, bands), C(:,2))*BGR(1)/BGR(2);
%         end
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Normalise the image accordingly using a sin
    %function to avoid highlight blowout
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    maxt = reshape(max(max(IMAGE, [], 1), [], 2), 1, 3);
    for i = 1:3
        IMAGE(:, :, i) = sin(IMAGE(:, :, i) /maxt(i)*pi/2);
    end

end