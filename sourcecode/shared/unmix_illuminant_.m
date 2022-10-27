%% Illuminant unmixing routine
%
%% Syntax:
%   HSZ = unmix_illuminant_(HS, Endmembers);
%   HSZ = unmix_illuminant_(HS, Endmembers, options);
% 
%% Description:
%     Unmixes the illuminants in an HSZ file so as to refer to the library in 
%     Endmembers. If the HSZ is non-indexed, it converts it into an indexed file.
% 
%% Input:
%     HS: Input Scyllarus data structure
%     options: Struct containing the following fields
%           numCanonicalIlluminants.: Number of endmembers used for unmixing each material
%           numIlluminants: Number of illuminants used for the indexing of the
%                   light in the image.
%           PSFFactor: Factor used for the point spread function. This is
%                   applied to enforce smoothness on the coefficients recovered by the
%                   L-2 unmixing method.
%
%% Output:
%     HSZ: Scyllarus data structure indexed and unmixed to the endmember 
%         library. This is RAW encoded. For encodings other than RAW, use the
%         encode_HSZ routine.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.6
% Last Update Date: 28 July 2014


function HSZ = unmix_illuminant_(HSPipeline, Endmembers, options)

    bands       = HSPipeline.HDR.bands;
    [rows cols] = size(HSPipeline.S.Factor);
    old_hdr     = HSPipeline.HDR;
    [wl, ~]     = wavelength_subset_(HSPipeline.HDR.wavelength, Endmembers.HDR.wavelength);
    HSPipeline  = translate_HSZ_(HSPipeline, wl);
    options.EndmemberIndexed = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Verify input parameters. Note that the variables options.NumIlluminants
    %   and options LightIndexed have not been implemented yet!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('options', 'var') || ~isfield(options,'numMaterials') || options.numMaterials > bands || ...
            options.numIlluminants < 1 || numel(options.numIlluminants) ~= 1
        options.numIlluminants = min(bands,5);    %Find five materials per pixel by default
    end
    
    if ~isfield(options,'numCanonicalIlluminants.')   
        options.numCanonicalIlluminants = min(bands,4);
    elseif  options.numCanonicalIlluminants > bands || options.numCanonicalIlluminants < 1
            options.numCanonicalIlluminants = min(bands, 4);    %Find four endmembers per material
    end
    
    if ~isfield(options,'IndexingMethod') || ~strcmp(options.IndexingMethod, 'KM')
         options.IndexingMethod = 'DA';
    end
    
    if ~isfield(options,'PSFFactor')|| options.PSFFactor<=0 || numel(options.PSFFactor) ~= 1
         options.PSFFactor = min([10, rows, cols]);
    end
    
    if ~isfield(options,'DEBUG') || numel(options.DEBUG)~=1 || options.DEBUG > 5 || options.DEBUG < 0
        options.DEBUG = 1;
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the material indexation if necesary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if HSPipeline.HDR.IndexedL == 0
        %Setup the parameters
        HSPipeline.L.Reflectance       = HSPipeline.L.Elements;
        [rows cols bands]              = size(HSPipeline.L.Reflectance);
        options.numMaterials           = min(bands, options.numIlluminants);
        [~, HSPipeline.L.Elements, EA] = recover_materials(HSPipeline.L.Reflectance, options.IndexingMethod, options.DEBUG); 
        [materials,~]                  = size(HSPipeline.L.Elements);
        
        %Check that there are enough highlights in the scene. Update options.numMaterials otherwise
       
        fprintf('Performing the indexation to %d illuminants.\n', materials);
        %Get the abundance and indexes matrices
        [~, HSPipeline.L.ElementAbundanceIndexes] = sort(EA, 3, 'descend');
        HSPipeline.L.ElementAbundanceIndexes      = HSPipeline.L.ElementAbundanceIndexes(:, :, 1:materials);
        %Get the least squares solution for the most abundant materials
        for i = 1:rows
            for j = 1:cols
                HSPipeline.L.ElementAbundances(i, j, :)=lsqnonneg(...
                    HSPipeline.L.Elements(reshape(HSPipeline.L.ElementAbundanceIndexes(i, j, :), 1, materials), :)',...
                    reshape(HSPipeline.L.Reflectance(i, j, :),bands, 1));
            end 
        end
        
        %Compensate for zero-values
        PSF = fspecial('gaussian', ceil(min(cols ,rows)/options.PSFFactor), ceil(min(cols, rows)/options.PSFFactor*3));
        for k = 1:materials
            Blurred = imfilter(reshape(HSPipeline.L.ElementAbundances(:,:,k), [rows cols]),PSF,'circular','conv');
            HSPipeline.L.ElementAbundances(:,:,k) = Blurred;
        end
        %Finish with the process
        [~, Indx] = sort(HSPipeline.L.ElementAbundances,3,'descend');
        HSPipeline.S.ElementAbundanceIndexes = HSPipeline.L.ElementAbundanceIndexes(:, :, Indx(1:options.numIlluminants));
        HSPipeline.S.ElementAbundances       = HSPipeline.L.ElementAbundances(:, :, Indx(1:options.numIlluminants));
        HSPipeline.HDR.numElementsL          = options.numIlluminants;
        HSPipeline.L                         = rmfield(HSPipeline.L, 'Reflectance');
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the final umixing    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HSPipeline.HDR.numEndmembersL = min([options.numCanonicalIlluminants bands Endmembers.HDR.numEndmembers]);
    [HSPipeline.L.EndmemberAbundances, HSPipeline.L.EndmemberAbundanceIndexes, HSPipeline.L.Endmembers, ~] = ...
        L2_unmixing(HSPipeline.L.Elements, HSPipeline.HDR.wavelength, Endmembers.Endmember, ...
        Endmembers.HDR.wavelength, HSPipeline.HDR.numEndmembersL);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Add the variables to the header 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HSPipeline.HDR.IndexedL              = 1;
    HSPipeline.HDR.EndmemberIndexedL     = 1;
    HSPipeline.HDR.numEndmembersL        = options.numCanonicalIlluminants;
    HSPipeline.HDREndmembersL            = Endmembers.HDR;
    HSPipeline.HDREndmembersL.wavelength = HSPipeline.HDR.wavelength;
    HSPipeline.HDR.bands                 = length(HSPipeline.HDR.wavelength);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the encoding if necesary  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HSPipeline = encode_HSZ(HSPipeline, struct('IlluminantEncoding',old_hdr.EncodingL, ...
            'MaterialEncoding', old_hdr.EncodingS, 'SpecularityEncoding', old_hdr.EncodingK));
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Transfer the struct to the final variable and exit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HSZ = HSPipeline;
    
end