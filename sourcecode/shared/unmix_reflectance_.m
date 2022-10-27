% Syntax:
%   HSZ = unmix_reflectance_(HSPipeline, Endmembers);
%   HSZ = unmix_reflectance_(HSPipeline, Endmembers, options);
% 
% Description:
%   This function is used to unmixe the materials in an HSZ file so as to refer to the library in Endmembers.
%   If the HSZ is non-indexed, it converts it into an indexed file.
% 
% Input:
%     HSPipeline: Input NICTA pipeline file
%     options: Struct containing the following fields
%       Endmembers: Input struct containing the library
%       numEndmembers: Number of endmembers used for unmixing each material
%       PSFFactor: Factor used for the point spread function. This is
%               applied to enforce smoothness on the coefficients recovered by the
%               L-2 unmixing method.
% 
% Output:
%     HSZ: NICTA pipeline file which has been indexed and unmixed to the endmember 
%        library. This is RAW encoded. For encodings other than RAW, use the
%          encode_HSZ routine.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.6
% Last Update Date: 28 July 2014

function HSZ = unmix_reflectance_(HSPipeline, Endmembers, options)
    
    bands        = HSPipeline.HDR.bands;
    [rows, cols] = size(HSPipeline.S.Factor);
    old_hdr      = HSPipeline.HDR;
    Endmembers.HDR.wavelength = reshape(Endmembers.HDR.wavelength, [1, length(Endmembers.HDR.wavelength)]);
    if ~isequal(HSPipeline.HDR.wavelength, Endmembers.HDR.wavelength);
        wl = HSPipeline.HDR.wavelength;
    else
        [wl, ~] = wavelength_subset_(HSPipeline.HDR.wavelength, Endmembers.HDR.wavelength);
    end
    
    HSPipeline  = translate_HSZ_(HSPipeline, wl);
    options.EndmemberIndexed = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Verify input parameters. Note that the variables options.NumIlluminants
    %   and options LightIndexed have not been implemented yet!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('options', 'var') || ...
       ~isfield(options, 'numMaterials') || ...
       options.numMaterials > bands || ...
       options.numMaterials < 1 || ...
       numel(options.numMaterials) ~= 1
        options.numMaterials = min(bands, 5);    %Find five materials per pixel by default
    end
    
    if ~isfield(options, 'numEndmembers')
        options.numEndmembers = min(bands, 4);
    elseif  options.numEndmembers > bands || options.numEndmembers < 1
            options.numEndmembers = min(bands, 4);    %Find four endmembers per material
    end
    
    if ~isfield(options, 'IndexingMethod') || ~strcmpi(options.IndexingMethod, 'KM')
         options.IndexingMethod = 'DA';
    end
    
    if ~isfield(options, 'PSFFactor')|| options.PSFFactor <= 0 || numel(options.PSFFactor) ~= 1
         options.PSFFactor = min([10, rows, cols]);
    end
    
    if ~isfield(options, 'DEBUG') || numel(options.DEBUG) ~= 1 || options.DEBUG > 5 || options.DEBUG < 0
        options.DEBUG = 1;
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the material indexation if necesary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if HSPipeline.HDR.IndexedS == 0
        %Setup the parameters
        HSPipeline.S.Reflectance       = HSPipeline.S.Elements;
        [rows, cols, bands]              = size(HSPipeline.S.Reflectance);
        options.numMaterials           = min(bands, options.numMaterials);        
        %Recover the Material basis
        [~, HSPipeline.S.Elements, EA] = recover_materials(HSPipeline.S.Reflectance, options.IndexingMethod, options.DEBUG); 
        [materials, ~]                 = size(HSPipeline.S.Elements);
        %Check that there are enough highlights in the scene. Update options.numMaterials otherwise       
        fprintf('Performing the indexation to %d materials.\n', materials);
        %Get the abundance and indexes matrices
        [~, HSPipeline.S.ElementAbundanceIndexes] = sort(EA, 3, 'descend');
        HSPipeline.S.ElementAbundanceIndexes      = HSPipeline.S.ElementAbundanceIndexes(:, :, 1:materials);
        
        %Get the least squares solution for the most abundant materials
        for i = 1:rows
            for j = 1:cols
                HSPipeline.S.ElementAbundances(i, j, :)=lsqnonneg(...
                    HSPipeline.S.Elements(reshape(HSPipeline.S.ElementAbundanceIndexes(i, j, :), 1, materials), :)', ...
                    reshape(HSPipeline.S.Reflectance(i, j, :), bands, 1));
            end 
        end
        
        %Compensate for zero-values
        PSF = fspecial('gaussian', ceil(min(cols, rows)/options.PSFFactor), ceil(min(cols, rows)/options.PSFFactor*3));
        
        for k = 1:materials
            Blurred = imfilter(reshape(HSPipeline.S.ElementAbundances(:, :, k), [rows cols]), PSF, 'circular', 'conv');
            HSPipeline.S.ElementAbundances(:, :, k) = Blurred;
        end
        
        %Finish with the process
        [~, Indx]                            = sort(HSPipeline.S.ElementAbundances, 3, 'descend');
        HSPipeline.S.ElementAbundanceIndexes = HSPipeline.S.ElementAbundanceIndexes(:, :, Indx(1:options.numMaterials));
        HSPipeline.S.ElementAbundances       = HSPipeline.S.ElementAbundances(:, :, Indx(1:options.numMaterials));
        HSPipeline.HDR.numElementsS          = options.numMaterials;
        HSPipeline.S                         = rmfield(HSPipeline.S, 'Reflectance');
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the final umixing    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HSPipeline.HDR.numEndmembersS = min([options.numEndmembers bands Endmembers.HDR.numEndmembers]);
    [HSPipeline.S.EndmemberAbundances, HSPipeline.S.EndmemberAbundanceIndexes, HSPipeline.S.Endmembers, ~] = ...
        L2_unmixing_(HSPipeline.S.Elements, HSPipeline.HDR.wavelength, Endmembers.Endmembers, ...
        Endmembers.HDR.wavelength, HSPipeline.HDR.numEndmembersS);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Add the variables to the header 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HSPipeline.HDR.IndexedS              = 1;
    HSPipeline.HDR.EndmemberIndexedS     = 1;
    HSPipeline.HDR.numEndmembersS        = options.numEndmembers;
    HSPipeline.HDREndmembersS            = Endmembers.HDR;
    HSPipeline.HDREndmembersS.wavelength = HSPipeline.HDR.wavelength;
    HSPipeline.HDR.bands                 = length(HSPipeline.HDR.wavelength);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the encoding if necesary  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HSZ = encode_HSZ_(HSPipeline, struct('IlluminantEncoding', old_hdr.EncodingL, ...
            'MaterialEncoding', old_hdr.EncodingS, 'SpecularityEncoding', old_hdr.EncodingK));

end