% Encoding routine for the HSZ data structure
%
% Syntax:
%     HSZ = encode_HSZ_(HSPipeline)
%     HSZ = encode_HSZ_(HSPipeline, options)
%
% Description:
%     Encodes an HSZ structure so as to convert it from RAW to a Gaussian mixture model 
%     (GMM) or a spline (NURBS).
%
% Input:
%     HSPipeline: Input RAW encoded HSZ structure
%     options: Structure with the following fields
%         IlluminantEncoding: Determines the manner in which the spectra on 
%             HSZ.L.Elements and HSZ.L.Endmembers 
%             is encoded. The default is 'NURBS'
%         MaterialEncoding: Determines the manner in which the spectra on 
%             HSZ.S.Elements and HSZ.S.Endmembers 
%             is encoded. The default is 'NURBS'
%         SpecularityEncoding: Encoding scheme used for the highlights at 
%             output. The default is 'NURBS'
%         numGMMsL: Number of mixtures used for the GMM encoding of the
%              illuminant
%         numGMMsS: Number of mixtures used for the GMM encoding of the
%              reflectance
%         numGMMsK: Number of mixtures used for the GMM encoding of the
%              specularity
%
% Output:
%   HSZ: An HSZ structure encoded according to the options above.
% 
% See also
%       eval_HSZ, get_nurbs_, get_gmm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2014, All Rights Reserved.
% Author: Antonio Robles-Kelly. 
% Version: 1.0.6
% Last Update Date: 28 Jan 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HSZ = encode_HSZ_(HSPipeline, options)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check the input variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bands = length(HSPipeline.HDR.wavelength);

    if ~exist('options', 'var') || ~isfield(options, 'MaterialEncoding') ||...
            (~strcmpi(options.MaterialEncoding, 'RAW') && ~strcmpi(options.MaterialEncoding, 'GMM'))
        options.MaterialEncoding = 'NURBS';    %Encode the materials using raw measurements
    end
    if ~isfield(options, 'IlluminantEncoding') || (~strcmpi(options.IlluminantEncoding, 'RAW') && ...
            ~strcmpi(options.IlluminantEncoding, 'GMM'))
        options.IlluminantEncoding = 'NURBS';    %Encode the illumiants using raw measurements
    end
    if ~isfield(options, 'SpecularityEncoding') || (~strcmpi(options.SpecularityEncoding, 'RAW') && ...
            ~strcmpi(options.SpecularityEncoding, 'GMM'))
        options.SpecularityEncoding = 'NURBS';    %Encode the specularities using raw measurements
    end
    if strcmpi(options.MaterialEncoding, 'GMM') && ( ~isfield(options, 'numGMMsS') || ...
            options.numGMMsS<1 || options.numGMMsS>bands )
         options.numGMMsS = min(bands, 4);
    end
    if strcmpi(options.SpecularityEncoding, 'GMM') && ( ~isfield(options, 'numGMMsK') ||...
            options.numGMMsK<1 || options.numGMMsK>bands )
         options.numGMMsK = min(bands, 4);
    end
    if strcmpi(options.IlluminantEncoding, 'GMM') || ( ~isfield(options, 'numGMMsL') || ...
            options.numGMMsL<1 || options.numGMMsL>bands )
         options.numGMMsL = min(bands, 4);
    end
    if ~isfield(options, 'DEBUG') || numel(options.DEBUG)~=1 || ...
            options.DEBUG > 5 || options.DEBUG < 0
        options.DEBUG = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the encoding for non-RAW formats
    %Note that, for the NURBS, the last row of the Elements/Endmember matrix
    %correponds to the control points on the wavelength space.
    %For the GMMs, the ordering is weights, means and standard deviations.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Encode the variables using NURBS
    %Do the following if the variables haven't been encoded

    %Do the indexed cases
    if strcmpi(options.MaterialEncoding, 'NURBS') 
        if HSPipeline.HDR.IndexedS == 1
            [mats, ~] = size(HSPipeline.S.Elements);
            [HSPipeline.S.ElementKnots, HSPipeline.S.Elements, HSPipeline.S.ElementCP] = get_nurbs_(reshape(HSPipeline.S.Elements, [mats 1 bands]), HSPipeline.HDR.wavelength, 2);
            HSPipeline.S.Elements  = reshape(HSPipeline.S.Elements,mats,length(HSPipeline.S.ElementCP));
            HSPipeline.S.ElementCP = HSPipeline.S.ElementCP';
            HSPipeline.HDR.degreeNURBSS = 2;
        else
            [HSPipeline.S.ElementKnots, CP_REF, CP_WAVE] = get_nurbs_(HSPipeline.S.Elements, HSPipeline.HDR.wavelength, 2);
            HSPipeline.S.Elements  = CP_REF;
            HSPipeline.S.ElementCP = CP_WAVE';
        end
        if HSPipeline.HDR.EndmemberIndexedS == 1
            [mats, ~] = size(HSPipeline.S.Endmembers);
            [HSPipeline.S.EndmemberKnots, HSPipeline.S.Endmembers, HSPipeline.S.EndmemberCP] = get_nurbs_(reshape(HSPipeline.S.Endmembers, ...
                [mats 1 length(HSPipeline.HDREndmembersS.wavelength)]), HSPipeline.HDREndmembersS.wavelength, 2);
            HSPipeline.S.Endmembers  = reshape(HSPipeline.S.Endmembers,mats,length(HSPipeline.S.EndmemberCP));
            HSPipeline.S.EndmemberCP = HSPipeline.S.EndmemberCP';
            HSPipeline.HDREndmembersS.degreeNURBS = 2;
        end
    end
    
    %Encode the Specularity
    %Do the following if the variables haven't been encoded
    %Do the indexed cases
    if strcmpi(options.SpecularityEncoding, 'NURBS')  
        if HSPipeline.HDR.IndexedK == 1
            [mats, ~] = size(HSPipeline.K.Elements);
            [HSPipeline.K.ElementKnots, CP_REF, CP_WAVE] = get_nurbs_(reshape(HSPipeline.K.Elements, [mats 1 bands]), HSPipeline.HDR.wavelength,2);
            [~, ~, cps] = size(CP_REF);
             HSPipeline.K.Elements  = reshape(CP_REF, [mats cps]);
             HSPipeline.K.ElementCP = CP_WAVE';
             HSPipeline.HDR.degreeNURBSK = 2;
        else
            [HSPipeline.K.ElementKnots, CP_REF, CP_WAVE] = get_nurbs_(HSPipeline.K.Elements, HSPipeline.HDR.wavelength,2);
            HSPipeline.K.Elements  = CP_REF;
            HSPipeline.K.ElementCP = CP_WAVE';
        end
    end
    %Encode the Illuminant
    %Do the following if the variables haven't been encoded
    if strcmpi(options.IlluminantEncoding, 'NURBS')
        if HSPipeline.HDR.IndexedL == 1
            [mats, ~] = size(HSPipeline.L.Elements);
            [HSPipeline.L.ElementKnots, CP_REF, CP_WAVE] = get_nurbs_(reshape(HSPipeline.L.Elements, [mats 1 bands]), HSPipeline.HDR.wavelength, 2);
            [~, ~, cps] = size(CP_REF);
            HSPipeline.L.Elements  = reshape(CP_REF, [mats cps]);
            HSPipeline.L.ElementCP = CP_WAVE';
            HSPipeline.HDR.degreeNURBSL = 2;
        else
            [HSPipeline.L.ElementKnots, CP_REF, CP_WAVE] = get_nurbs_(HSPipeline.L.Elements, HSPipeline.HDR.wavelength, 2);
            HSPipeline.L.Elements  = CP_REF;
            HSPipeline.L.ElementCP = CP_WAVE';
        end
        if HSPipeline.HDR.EndmemberIndexedL == 1
            [mats, ~] = size(HSPipeline.L.Endmembers);
            [HSPipeline.L.EndmemberKnots, CP_REF, CP_WAVE] = get_nurbs_(reshape(HSPipeline.L.Endmembers, ...
                [mats 1 length(HSPipeline.HDREndmembersL.wavelength)]), HSPipeline.HDREndmembersL.wavelength, 2);
            [~, ~, cps] = size(CP_REF);
            HSPipeline.L.Endmembers  = reshape(CP_REF, [mats cps]);
            HSPipeline.L.EndmemberCP = CP_WAVE';
        end
        
    end
        
    
    %Encode the variables using GMMs    
    %Do the following if the variables haven't been encoded

    %Do the indexed cases
    if strcmpi(options.MaterialEncoding, 'GMM') 
        HSPipeline.HDR.numGMMsS = options.numGMMsS;
        if HSPipeline.HDR.IndexedS == 1
            [mats, ~] = size(HSPipeline.S.Elements);
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(reshape(HSPipeline.S.Elements, [mats 1 bands]), HSPipeline.HDR.wavelength, HSPipeline.HDR.numGMMsS);
            HSPipeline.S.Elements    = reshape(MIXCOEFF, [mats HSPipeline.HDR.numGMMsS]);
            HSPipeline.S.ElementMean = reshape(MEAN, [mats HSPipeline.HDR.numGMMsS]);
            HSPipeline.S.ElementStd  = reshape(STD, [mats HSPipeline.HDR.numGMMsS]);
        else
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(HSPipeline.S.Elements, HSPipeline.HDR.wavelength, ...
                HSPipeline.HDR.numGMMsS);
            HSPipeline.S.Elements    = MIXCOEFF;
            HSPipeline.S.ElementMean = MEAN;
            HSPipeline.S.ElementStd  = STD;
        end
        if HSPipeline.HDR.EndmemberIndexedS == 1
            [mats, ~] = size(HSPipeline.S.Endmembers);
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(reshape(HSPipeline.S.Endmembers, ...
                [mats 1 length(HSPipeline.HDREndmembersS.wavelength)]), HSPipeline.HDREndmembersS.wavelength, HSPipeline.HDR.numGMMsS);
            HSPipeline.S.Endmembers    = reshape(MIXCOEFF, [mats HSPipeline.HDR.numGMMsS]);
            HSPipeline.S.EndmemberMean = reshape(MEAN, [mats HSPipeline.HDR.numGMMsS]);
            HSPipeline.S.EndmemberStd  = reshape(STD, [mats HSPipeline.HDR.numGMMsS]);
        end        
    end
    
    %Encode the Specularity
    %Do the following if the variables haven't been encoded
    %Do the indexed cases
    if strcmpi(options.SpecularityEncoding, 'GMM')
        if HSPipeline.HDR.IndexedK == 1
            HSPipeline.HDR.numGMMsK = options.numGMMsK;
            [mats, ~] = size(HSPipeline.K.Elements);
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(reshape(HSPipeline.K.Elements, [mats 1 bands]), HSPipeline.HDR.wavelength, HSPipeline.HDR.numGMMsK);
            HSPipeline.K.Elements    = reshape(MIXCOEFF, [mats HSPipeline.HDR.numGMMsK]);
            HSPipeline.K.ElementMean = reshape(MEAN, [mats HSPipeline.HDR.numGMMsK]);
            HSPipeline.K.ElementStd  = reshape(STD, [mats HSPipeline.HDR.numGMMsK]);
        else
            HSPipeline.HDR.numGMMsK = options.numGMMsK;
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(HSPipeline.K.Elements, HSPipeline.HDR.wavelength, HSPipeline.HDR.numGMMsK);
            HSPipeline.K.Elements    = MIXCOEFF;
            HSPipeline.K.ElementMean = MEAN;
            HSPipeline.K.ElementStd  = STD;
        end
    end
    %Encode the Illuminant
    %Do the following if the variables haven't been encoded
    %Do the indexed cases
    if strcmpi(options.IlluminantEncoding, 'GMM')
        if HSPipeline.HDR.IndexedL == 1
            HSPipeline.HDR.numGMMsL = options.numGMMsL;
            [mats, ~] = size(HSPipeline.L.Elements);
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(reshape(HSPipeline.L.Elements, [mats 1 bands]), HSPipeline.HDR.wavelength, HSPipeline.HDR.numGMMsL);
            HSPipeline.L.Elements    = reshape(MIXCOEFF, [mats HSPipeline.HDR.numGMMsL]);
            HSPipeline.L.ElementMean = reshape(MEAN, [mats HSPipeline.HDR.numGMMsL]);
            HSPipeline.L.ElementStd  = reshape(STD, [mats HSPipeline.HDR.numGMMsL]);
        else
            HSPipeline.HDR.numGMMsL = options.numGMMsL;
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(HSPipeline.L.Elements, HSPipeline.HDR.wavelength, ...
                HSPipeline.HDR.numGMMsL);
            HSPipeline.L.Elements    = MIXCOEFF;
            HSPipeline.L.ElementMean = MEAN;
            HSPipeline.L.ElementStd  = STD;
        end
        if HSPipeline.HDR.EndmemberIndexedL == 1
            HSPipeline.HDR.numGMMsL = options.numGMMsL;
            [mats, ~] = size(HSPipeline.L.Endmembers);
            [MIXCOEFF, MEAN, STD] = get_gaussian_mixture_(reshape(HSPipeline.L.Endmembers, ...
                [mats 1 length(HSPipeline.HDREndmembersL.wavelength)]), HSPipeline.HDREndmembersL.wavelength, ...
                HSPipeline.HDR.numGMMsL);
            HSPipeline.L.Endmembers    = reshape(MIXCOEFF, [mats HSPipeline.HDR.numGMMsL]);
            HSPipeline.L.EndmemberMean = reshape(MEAN, [mats HSPipeline.HDR.numGMMsL]);
            HSPipeline.L.EndmemberStd  = reshape(STD, [mats HSPipeline.HDR.numGMMsL]);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Add variables to the header
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    HSPipeline.HDR.EncodingS = options.MaterialEncoding;
    HSPipeline.HDR.EncodingK = options.SpecularityEncoding;
    HSPipeline.HDR.EncodingL = options.IlluminantEncoding;
    
    HSZ = HSPipeline;
    
end