% Syntax:
%   HSZ = eval_HSZ_(HSPipeline);
%   HSZ = eval_HSZ_(HSPipeline, WAVELENGTH);
%
% Description:
%   Evaluates an NICTA pipeline structure so as to return another HSZ struct encoded in RAW measurements over the
%   wavelengths on WAVELENGTH. Note that this routine has no effect if the HSZ structure has been RAW encoded.
%
% Input:
%   WAVELENGTH: Vector of target wavelengths. These have to be within the scope of
%               HSPipeline.HDR.wavelength. If this parameter is not provided, 
%               HSPipeline.HDR.wavelength is used instead. 
%   HSPipeline: Data structure to be evaluated
%
% Output:
%   HSZ: Data structure encoded in RAW format over the wavelength vector.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2014, All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.6
% Last Update Date: 17 July 2014

function HSZ = eval_HSZ_(HSPipeline, WAVELENGTH)

    %Check the input parameters
    sub_wavelength = 0;
    if ~exist('WAVELENGTH', 'var')
        if isfield(HSPipeline.HDR, 'wavelength')
            WAVELENGTH = HSPipeline.HDR.wavelength;
        else
            error('Wavelength information is missing, please check');
        end
        sub_wavelength = 1;
    end
     
    %Make a copy of HSPipeline
    HSZ = HSPipeline;
    
    if sub_wavelength == 0
        %Assure the evaluation is given by a row vector
        WAVELENGTH = WAVELENGTH((WAVELENGTH<=max(HSPipeline.HDR.wavelength) & WAVELENGTH>=min(HSPipeline.HDR.wavelength)));
        if isempty(WAVELENGTH)
            WAVELENGTH = min(HSPipeline.HDR.wavelength) : (max(HSPipeline.HDR.wavelength) - min(HSPipeline.HDR.wavelength))/33 : ...
                         max(HSPipeline.HDR.wavelength);
        end
    end
    
    bands = length(WAVELENGTH);
    wavelength = reshape(WAVELENGTH, [bands 1]);
    HSZ.HDR.wavelength = reshape(WAVELENGTH, [1, bands]);
    HSZ.HDR.bands = bands;
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start the evaluation. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the NURBS first
    %Start with the reflectance
    if strcmpi(HSPipeline.HDR.EncodingS, 'NURBS')
        if HSPipeline.HDR.IndexedS == 1
            [mats, cols] = size(HSPipeline.S.Elements);
            Q = eval_nurbs(HSPipeline.S.ElementKnots, wavelength, reshape(HSPipeline.S.Elements,[mats 1 cols]), HSPipeline.S.ElementCP', HSPipeline.HDR.degreeNURBSS); 
            HSZ.S.Elements = reshape(Q, [mats bands]);
            %Update the header
            HSZ.HDR.EncodingS = 'RAW';
            HSZ.S = rmfield(HSZ.S, 'ElementCP');
            HSZ.S = rmfield(HSZ.S, 'ElementKnots');
        else
            HSZ.S.Elements = eval_nurbs(HSPipeline.S.ElementKnots, wavelength, HSPipeline.S.Elements, HSPipeline.S.ElementCP', HSPipeline.HDR.degreeNURBSS);
            %Update the header
            HSZ.HDR.EncodingS = 'RAW';
            HSZ.S = rmfield(HSZ.S, 'ElementCP');
            HSZ.S = rmfield(HSZ.S, 'ElementKnots');
        end
        if HSPipeline.HDR.EndmemberIndexedS == 1
            [mats, cols] = size(HSPipeline.S.Endmembers);
            Q = eval_nurbs(HSPipeline.S.ElementKnots, wavelength, reshape(HSPipeline.S.Endmembers,[mats 1 cols]), HSPipeline.S.ElementCP', HSPipeline.HDR.degreeNURBSS); 
            HSZ.S.Endmembers = reshape(Q, [mats bands]);
            HSZ.S = rmfield(HSZ.S, 'EndmemberCP');
            HSZ.S = rmfield(HSZ.S, 'EndmemberKnots');
        end
        HSZ.HDR = rmfield(HSZ.HDR, 'degreeNURBSS');
    end

     %Encode the Specularity
    if strcmpi(HSPipeline.HDR.EncodingK, 'NURBS')
        if HSPipeline.HDR.IndexedK == 1
            [mats, cols] = size(HSPipeline.K.Elements);
            Q = eval_nurbs(HSPipeline.K.ElementKnots, wavelength, reshape(HSPipeline.K.Elements, [mats 1 cols]), HSPipeline.K.ElementCP', HSPipeline.HDR.degreeNURBSK);   
            HSZ.K.Elements = reshape(Q, [mats bands]);
            %Update the header
            HSZ.HDR.EncodingK = 'RAW';
            HSZ.K = rmfield(HSZ.K, 'ElementCP');
            HSZ.K = rmfield(HSZ.K, 'ElementKnots');
        else            
            HSZ.K.Elements = eval_nurbs(HSPipeline.K.ElementKnots, wavelength, HSPipeline.K.Elements, HSPipeline.K.ElementCP', HSPipeline.HDR.degreeNURBSK); 
            %Update the header
            HSZ.HDR.EncodingK = 'RAW';
            HSZ.K = rmfield(HSZ.K, 'ElementCP');
            HSZ.K = rmfield(HSZ.K, 'ElementKnots');
        end
        HSZ.HDR = rmfield(HSZ.HDR, 'degreeNURBSK');
    end
    
     %Encode the Illuminant
    if strcmpi(HSPipeline.HDR.EncodingL, 'NURBS')
        if HSPipeline.HDR.IndexedL == 1
            [mats, cols] = size(HSPipeline.L.Elements);
            %clear HSZ.L.Elements;
            Q = eval_nurbs(HSPipeline.L.ElementKnots, wavelength, reshape(HSPipeline.L.Elements, [mats 1 cols]), HSPipeline.L.ElementCP', HSPipeline.HDR.degreeNURBSL);   
            HSZ.L.Elements = reshape(Q, [mats bands]);
            %Update the header
            HSZ.HDR.EncodingL = 'RAW';
            HSZ.L = rmfield(HSZ.L, 'ElementCP');
            HSZ.L = rmfield(HSZ.L, 'ElementKnots');
        else
            HSZ.L.Elements = eval_nurbs(HSPipeline.L.ElementKnots, wavelength, HSPipeline.L.Elements, HSPipeline.L.ElementCP', HSPipeline.HDR.degreeNURBSL); 
            %Update the header
            HSZ.HDR.EncodingL = 'RAW';
            HSZ.L = rmfield(HSZ.L, 'ElementCP');
            HSZ.L = rmfield(HSZ.L, 'ElementKnots');
        end
        if HSPipeline.HDR.EndmemberIndexedL == 1
            HSZ.HDR.EncodingL = 'RAW';
            [mats, cols] = size(HSPipeline.L.Endmembers);
            clear HSZ.L.Endmembers;
            Q = eval_nurbs(HSPipeline.L.EndmemberKnots, wavelength, reshape(HSPipeline.L.Endmembers, [mats 1 cols]), HSPipeline.L.EndmemberCP', HSPipeline.HDR.degreeNURBSL);  
            HSZ.L.Endmembers = reshape(Q, [mats bands]);
            HSZ.L = rmfield(HSZ.L, 'EndmemberCP');
            HSZ.L = rmfield(HSZ.L, 'EndmemberKnots');
        end
        HSZ.HDR = rmfield(HSZ.HDR, 'degreeNURBSL');
     end
     
    %Encode the variables using GMMs
    %Encode the reflectance
     
    if strcmpi(HSPipeline.HDR.EncodingS, 'GMM')
        if HSPipeline.HDR.IndexedS == 1
            [mats, ~]=size(HSPipeline.S.Elements);
            %clear HSZ.S.Elements;
            Q = eval_gaussian_mixture_(HSPipeline.S.Elements, HSPipeline.S.ElementMean, HSPipeline.S.ElementStd, wavelength);
            HSZ.S.Elements = reshape(Q, [mats bands]);
            HSZ.S = rmfield(HSZ.S, 'ElementStd');
            HSZ.S = rmfield(HSZ.S, 'ElementMean');
        else
            HSZ.S.Elements = eval_gaussian_mixture_(HSPipeline.S.Elements, HSPipeline.S.ElementMean, HSPipeline.S.ElementStd, wavelength);
            HSZ.S = rmfield(HSZ.S, 'ElementStd');
            HSZ.S = rmfield(HSZ.S, 'ElementMean');
        end
        
        if HSPipeline.HDR.EndmemberIndexedS == 1
            [mats, ~]=size(HSPipeline.S.Endmembers);
            Q = eval_gaussian_mixture_(HSPipeline.S.Endmembers, HSPipeline.S.EndmemberMean, HSPipeline.S.EndmemberStd, wavelength);
            HSZ.S.Endmembers = reshape(Q, [mats bands]);
            HSZ.S = rmfield(HSZ.S, 'EndmemberStd');
            HSZ.S = rmfield(HSZ.S, 'EndmemberMean');
        end
        
        HSZ.HDR.EncodingS = 'RAW';
        HSZ.HDR = rmfield(HSZ.HDR, 'numGMMsS');
    end

    %Encode the Specularity
    if strcmpi(HSPipeline.HDR.EncodingK, 'GMM')
        HSZ.HDR.EncodingK = 'RAW';
        if HSPipeline.HDR.IndexedK == 1
            [mats, ~]=size(HSPipeline.K.Elements);
            %clear HSZ.K.Elements;
            Q = eval_gaussian_mixture_(HSPipeline.K.Elements, HSPipeline.K.ElementMean, HSPipeline.K.ElementStd, wavelength);      
            HSZ.K.Elements = reshape(Q, [mats bands])
        else
            HSZ.K.Elements = eval_gaussian_mixture_(HSPipeline.K.Elements, HSPipeline.K.ElementMean, HSPipeline.K.ElementStd, wavelength);                                              
        end
        
        HSZ.K = rmfield(HSZ.K, 'ElementStd');
        HSZ.K = rmfield(HSZ.K, 'ElementMean');
        HSZ.HDR = rmfield(HSZ.HDR, 'numGMMsK');
    end
    
    %Encode the Illuminant
    if strcmpi(HSPipeline.HDR.EncodingL, 'GMM')
        HSZ.HDR.EncodingL = 'RAW';
        if HSPipeline.HDR.IndexedL == 1
            [mats, ~]=size(HSPipeline.L.Elements);
            Q = eval_gaussian_mixture_(HSPipeline.L.Elements, HSPipeline.L.ElementMean, HSPipeline.L.ElementStd, wavelength);
            HSZ.L.Elements = reshape(Q, [mats bands]);
            HSZ.L = rmfield(HSZ.L, 'ElementStd');
            HSZ.L = rmfield(HSZ.L, 'ElementMean');
        else
            HSZ.L.Elements = eval_gaussian_mixture_(HSPipeline.L.Elements, HSPipeline.L.ElementMean, HSPipeline.L.ElementStd, wavelength);
            HSZ.L = rmfield(HSZ.L, 'ElementStd');
            HSZ.L = rmfield(HSZ.L, 'ElementMean');
        end
        if HSPipeline.HDR.EndmemberIndexedL == 1
            
            [mats, ~]=size(HSPipeline.L.Endmembers);
            clear HSZ.L.Endmembers;
            Q = eval_gaussian_mixture_(HSPipeline.L.Endmembers, HSPipeline.L.EndmemberMean, HSPipeline.L.EndmemberStd, wavelength);
            HSZ.L.Endmembers = reshape(Q, [mats bands]);
            HSZ.L = rmfield(HSZ.L, 'EndmemberStd');
            HSZ.L = rmfield(HSZ.L, 'EndmemberMean');
        end
        HSZ.HDR = rmfield(HSZ.HDR, 'numGMMsL');
    end
    %   end of function eval_HSZ_.m
 end