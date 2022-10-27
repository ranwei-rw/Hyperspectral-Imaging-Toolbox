% Scyllarus main routine script.
%
% Syntax:
%   HSZ = Scyllarus(I);
%   HSZ = Scyllarus(I, options);
%   HSZ = Scyllarus(filename);
%   HSZ = Scyllarus(filename, options);
%   HSZ = Scyllarus(I, options, Endmembers);
%   HSZ = Scyllarus(I, options, Endmembers, CanonicalIlluminants);
%   HSZ = Scyllarus(HS);
%   HSZ = Scyllarus(HS, options);
%   HSZ = Scyllarus(HS, options, Endmembers);
%   HSZ = Scyllarus(HS, options, Endmembers, CanonicalIlluminants);
%
% Description:
% This is the implementation of the Scyllarus pipeline for the dichromatic model, 
% which comprises the recovery of the specularities, shading factor, reflectance
% and illuminant power spectrum.
%
% Input:
%
%   I: hyperspectral image stored as a 3D array. 
%   HS: An HSZ data structure (see below)
%   options: Structure with the following fields
%         L: A structure defining the options for the illuminant recovery step
%               bitdepth: Is the data type for the spectral cube, i.e. number of bits per
%                   spectral measurement. By fault this is 16.
%               method: String specifying the method to be used. It can be any of the following
%                       'HRK': Employs the method of Huynh and Robles-Kelly (A Solution of the 
%                           Dichromatic Model for Multispectral Photometric Invariance, 
%                           International Journal of Computer Vision 2010).
%                       'FS': Uses the method of Finlayson and Schaefer (Convex and Non-convex 
%                           Illuminant Constraints for Dichromatic Colour Constancy, CVPR 2001).
%                       'GW': Uses the Grey World method.
%                       'SG': Uses the Shade of Grey method.
%                       'WP': Uses the White Patch method.
%                       '1stOGE': Uses the 1st order Grey Edge method.
%                       '2ndOGE': Uses the 2nd order Grey Edge method.
%               drate: Image downsampling rate for the Grey World, Shade of Grey, 
%                       White Patch and Grey Edge methods. The default is 1, 
%                       i.e. no downsampling.
%               shadeOfGreyOrder: The order of the L^p mean used for the Shade of Grey method.
%                       The default is 1.
%               alpha:   The value for the regularisation term used for the HRK
%                   (Huynh and Robles-Kelly) method. The default for this is 50. 
%               patches: Pre-selected patches. This could be a set of geometry data of patches
%                   with a format of (Top_left_y, top_left_x, height, 
%                   width). This can be left empty.
%               PSFFactor: The factor used for the PSF employed to 'soften' the distribution of abundances for
%                   the illuminants. The default is 10.
%               DEBUG: Defines the level of debugging information shown at execusion time (DEBUG<6).
%                   the default is 3.
%         R: A structure defining the options for the photometric parameter recovery step
%               method: String specifying the method to be used. This can be 'LS'
%                   for the linear least squares solver, 'KL' for the
%                   energy.
%               neigbourhoodsize (optional, used by the LS option): a threshold used to ignore the processing of regions with a small
%                   number of pixels. Those with a number of pixels smaller than or equal to this are considered
%                   diffuse. (default value is 5) ignoreThresh = 5;
%               grayThresh (optional, used by the LS option): a threshold used to determine whether a material is a shade of gray. If
%                   the reflectance spectra within a cluster does not deviate from a uniform spectrum (a flat
%                   line) beyond this threshold, then we will assume that it is a shade of gray and purely
%                   diffuse. (default value is 2) grayThresh = 2;
%               numclusters (option, used by the KL option): the number of clusters used for the K-means, the
%                   default is 20.
%               DEBUG: Defines the level of displaying debugging information. Default is 1, the least
%                   information will be given
%         IndexingMethod: Method used to index to materials. This can be either
%             'DA' (deterministic annealing) or 'KM' (k-means clustering).
%             The default is 'KM'
%         IlluminantEncoding: Determines the manner in which the spectra on HSZ.L.Elements and HSZ.L.Endmembers 
%             is encoded. The default is 'NURBS', but 'GMM' and 'RAW' can
%             also be used.
%         IlluminantIndexed: Determines whether the pixel wise illuminant is indexed to a set of 
%             illuminants "representative" of the lights in the scene. The default is 0 (no indexation)
%         CanonicalIlluminantIndexed: Turns on and off the indexing of the pixel-wise illuminant
%             power spectrum to canonical illuminants. The default is 0 (off).
%         numIlluminants: Number of lights used to represent each pixel. The default value 
%             is min(maxnumIlluminants, bands).
%         maxnumIlluminants: Number of target representative lights in the scene.
%             Note that this must be greater or equal to numIlluminants.  
%             The default value is 1, i.e. a single illuminant in the scene.
%         numCanonicalIlluminants: Number of target cannonical illuminants in which decompose the
%             each representative lights.
%         MaterialEncoding: Determines the manner in which the spectra on HSZ.S.Elements and HSZ.S.Endmembers 
%             is encoded. The default is 'NURBS', but 'GMM' and 'RAW' can
%             also be used.
%         MaterialIndexed: Determines whether the reflectance is to be decomposed into
%             the materials in the scene. The default is 1 (material indexed).
%         EndmemberIndexed: Turns on and off the indexation of the scene materials to the endmembers on
%             the Endmembers matrix.
%         maxnumMaterials: Maximum number of materials used for the indexation 
%             of the image as a whole. Note that this must be greater or equal to 
%             numMaterials. The default value is min(bands, 25).
%         numMaterials: Number of materials used to express the reflectance for each pixel.
%         numEndmembers: Number of endmembers that should be used to express each material.
%         SpecularityEncoding: Encoding scheme used for the highlights at output. 
%             The default is 'NURBS', but 'GMM' and 'RAW' can also be used.
%         SpecularityIndexed: Turn on and off the indexation of the highlights to a statistically
%             consistent set of spectral vectors. The default is 1 (on)
%         maxnumHighlights: Maximum number of spectral vectors used for the
%             indexation of the image as a whole. Note that this must be
%             greater or equal to numHighlights. The default is min(bands, 25).
%         numHighlights: Number of spectral vectors used for the indexation. The default is 1, which, 
%             by definition, sets the output HSZ.K.Elements to a uniform (all ones) spectral
%             vector. This also implies that the highlights are not wavelength dependant.
%         PSFFactor: The factor used for the PSF employed to 'soften' the distribution of abundances for
%             the materials, highlights and illuminants. The default is 10.
%   Endmembers: A structure with the following fields.
%         Endmembers: An array containing the endmember spectra
%         HDR: A structure containing the following fields
%               wavelength: vector of wavelengths.
%               Matxx: Material names
%   CanonicalIlluminants: A structure with the following fields.
%         Endmembers: An array containing the canonical illuminant spectra
%         HDR: A structure containing the following fields
%               wavelength: vector of wavelengths.
%               Matxx: Material names
%
% Output:
%
%   HSZ: A structure containing the following fields
%       HDR: A structure containing the fields corresponding to I.HDR and, 
%             in addition, the following fields
%             EncodingL: Encoding scheme for the data on HSZ.L.Elements and HSZ.L.Endmembers. This is 
%                 equivalent to options.IlluminantEncoding
%             IndexedL: Field equivalent to options.IlluminantIndexed.
%             EndmemberIndexedL: Field equivalent to options.CanonicalIlluminantIndexed.
%             numElementsL: Number of lights comprising the illuminant mixture for each pixel. This is, 
%                 ideally, the same as options.numIlluminants.
%             numEndmembersL: Number of cannonical illuminants comprising each of the representative
%                 lights in the scene.
%             degreeNURBSL: Degree of the polynomial used for the spline
%                 function used to represent the illuminant spectra. This has been set to 2.
%             numGMMsL: Number of Gaussians comprising the mixture used to represent 
%                 the illuminant spectra.
%             EncodingS: Encoding scheme for  HSZ.S.Elements and HSZ.S.Endmembers. This is 
%                 equivalent to options.MaterialEncoding.
%             IndexedS: Field equivalent to options.MaterialIndexed.
%             EndmemberIndexedS: Field equivalent to options.EndmemberIndexed.
%             numElementsS: Number of materials comprising each pixel reflectance.
%             numEndmembersS: Number of endmembers comprising each of the materials in the scene.
%             degreeNURBSS: Degree of the polynomial used for the spline
%                 function used to represent the reflectance spectra. This has been set to 2.
%             numGMMsS: Number of Gaussians comprising the mixture used to represent 
%                   the reflectance spectra.
%             EncodingK: Encoding scheme used for the specularities in the image.
%             IndexedK: Turn on and off the indexation of the highlights to a statistically
%                 consistent set of spectral vectors.
%             numElementsK: Number of spectral vectors used for the indexation. 
%             numGMMsK: Number of Gaussians comprising the mixture used to 
%                 represent the highlight spectra. 
%             degreeNURBSK: Degree of the polynomial used for the spline
%                 function used to represent the highlights. This has been set to 2.
%       L: A structure containing one or more of the following fields corresponding to the
%          illuminant.
%             Elements: Spectra for the representative illuminants in the
%             scene. 
%             ElementAbundances: Abundance matrix for the spectra on Elemements.
%             ElementAbundanceIndexes: Index matrix for the abundances of the element spectra. 
%             Endmembers: Spectra for the canonical illuminant library.
%             EndmemberAbundances: Canonical illuminant abundance matrix for the representative
%                 illuminants. This is a step further in the indexation of the spectra.
%             EndmemberAbundanceIndexes: Index matrix for the abundances of the canonical
%                 illuminant spectra.
%       R: A structure containing one or more of the following fields corresponding to the
%          reflectance.
%             Elements: Spectra for the materials in the scene.
%             ElementAbundances: Abundance matrix for the spectra on Elemements.
%             ElementAbundanceIndexes: Index matrix for the abundances of the element spectra. 
%             Endmembers: Spectra for the endmember library.
%             EndmemberAbundances: Endmember abundance matrix for the representative
%                 illuminants. This is a step further in the indexation of the spectra.
%             EndmemberAbundanceIndexes: Index matrix for the abundances of the endmember spectra.
%             Factor: Matrix of pixelwise scaling factors for the reflectance
%       K: A structure containing one or more of the following fields corresponding to the
%          specular highlights
%             Elements: Spectra for the set of basis spectral vectors spanning the
%                  specular highlights. This is only present if IndexedK ==
%                  1.
%             ElementAbundances: Abundance matrix for the spectra on Elemements.
%             ElementAbundanceIndexes: Index matrix for the abundances of the element spectra. 
%             Factor: Matrix of pixelwise specular coefficients.
% 
%  Examples
% 
% Read an FLA file and process it using Scyllarus with material indexation and 
% 5 materials per pixel
% 
%   I = FLAread('.\shared\samples\apples_small.fla');
%   options = struct('numMaterials', 5, 'DEBUG', 3);
%   HSZ = Scyllarus(I, options);
%   or
%   options = struct('numMaterials', 5, 'DEBUG', 3);
%   HSZ = Scyllarus('.\shared\samples\apples_small.fla', options);
% 
% Read an FLA file rescaling it to half is size and an SLZ file (end-member library file) 
% and process it using Scyllarus with both, material indexation and unmixing.
%
%   I = FLAread('.\shared\samples\apples_small.fla', 0.5);
%   Endmembers = SLZread('ISSA72.slz');
%   options = struct('MaterialIndexed', 1, 'numIlluminants', 3);
%   HSZ = Scyllarus(I, options, Endmembers);
%
% Read an FLA file and an SLZ file (end-member library file) and process it using
% Scyllarus with both, material indexation and unmixing.
%
%   I = FLAread('.\shared\samples\apples_small.fla');
%   Endmembers = SLZread('ISSA72.slz');
%   CanonicalIlluminants = SLZread('SFUIlluminants.slz');
%   options = struct('MaterialIndexed', 1, 'numIlluminants', 3);
%   HSZ = Scyllarus(I, options, Endmembers, CanonicalIlluminants);
%   HSZwrite('.\shared\samples\face.hsz', HSZ);
%
% Read an FLA file and an SLZ file (end-member library file), process it using
% Scyllarus with both, material indexation and unmixing and encode it using 
% Gaussian mixtures.
%
%   I = FLAread('.\shared\samples\apples_small.fla', 500, 620);
%   Endmembers = SLZread('XRiteCC154.slz');
%   CanonicalIlluminants = SLZread('SFUIlluminants.slz');
%   options = struct('MaterialIndexed', 1, 'numIlluminants', 3, 'MaterialEncoding', ...
%                    'GMM', 'SpecularityEncoding', 'GMM', 'IlluminantEncoding', 'GMM');
%   HSZ = Scyllarus(I, options, Endmembers, CanonicalIlluminants);
%   Irgb = recover_RGBImage(HSZ);
%
% Read an FLA file and process it using Scyllarus without material indexation and 
% RAW output by recovering a single illuminant with the generalised gray
% world method ('GWW'). Once the HSZ is computed, the radiance with the
% specularities, shading and illuminant removed is visualised as a trichromatic
% image and the global illuminant in the scene is ploted.
%
%   I = FLAread('.\shared\samples\apples_small.fla');
%   options = struct('MaterialIndexed', 0, 'SpecularityIndexed', 0, 'MaterialEncoding', ...
%                    'RAW', 'SpecularityEncoding', 'RAW', 'IlluminantEncoding', 'RAW', 'numIlluminants', 1);
%   options.L.method = 'GWW';
%   HSZ = Scyllarus(I, options);
%   I.I = HSZ.S.Elements;
%   Irgb = recover_RGBImage(I);
%   imtool(Irgb); 
%   plot(HSZ.HDR.wavelength, HSZ.L.Elements);
%
%
% See also
%       recover_RGBImage, HSZwrite, HSZread, SLZread, FLAread 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013-2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly. 
% Version: 0.7.0
% Date: 21 Nov 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HSPipeline = Scyllarus_(I, options, Endmembers, CanonicalIlluminants)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Set the method to that corresponding to the dichromatic model. This is 
    %The only method supported so far
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(I, 'HDR')
        %   when HDR is not a filed of S, we then need to test this is a filename containing a fla file
        if ischar(I)
            %   I is given as a string. Try to see whether it's a file location
            filename = I;
            if exist(filename, 'file') == 2
               %    when it's a file name, use a different 
               [~, ~, ext] = fileparts(filename);
               disp('A source file name is given. Importing data');
               switch lower(ext)
                   case '.fla'
                       I = FLAread(filename);
                       disp('Data is imported.');
                   case '.hsz'
                       I = HSZread(filename);
                       disp('Data is imported.');
                   case '.slz'
                       I = SLZread(filename);
                       disp('Data is imported.');
                   otherwise
                       error('Unknown parameter I is given');
               end
            else
                error('Unknown parameter I is given');
            end
        end
    end
    
    if isfield(I, 'S')               %Recover a raw image in case the input is an HSZ data struct
        disp('Reconstructing image.');
        I = reconstruct_image_(I);
    end
    
    HSPipeline.HDR        = I.HDR;
    options.method        = 'Dichromatic';
    [cols, rows, bands]   = size(I.I);
    HSPipeline.HDR.method = options.method;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Verify input parameters. Note that the variables options.NumIlluminants
    %   and options LightIndexed have not been implemented yet!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    if ~exist('options', 'var') || ~isfield(options, 'MaterialIndexed') || numel(options.MaterialIndexed) ~= 1
        options.MaterialIndexed = 1;    %This is the default, i.e. material indexation
    end
    
    if ~exist('Endmembers', 'var')
        Endmembers = [];
    end
    
    if ~exist('CanonicalIlluminants', 'var')
        CanonicalIlluminants = [];
    end
    
    if ~isfield(options, 'IlluminantIndexed') || options.IlluminantIndexed ~= 0 || numel(options.IlluminantIndexed) ~= 1
        options.IlluminantIndexed = 1;    %This is the default, i.e. the specularity is not wavelength dependant
    end
    
    if ~isfield(options, 'numIlluminants') || options.numIlluminants < 1 || options.numIlluminants > bands ...
                || numel(options.numIlluminants) ~= 1
        options.numIlluminants = 1;    %This is the default, i.e. a single illuminant for the whole image
    end
    
    if ~isfield(options, 'maxnumIlluminants') || options.maxnumIlluminants < 1 || ...
                numel(options.maxnumIlluminants) ~= 1 ||...
                options.numIlluminants > options.maxnumIlluminants 
            %options.maxnumIlluminants > bands || 
        options.maxnumIlluminants = min(options.numIlluminants, bands); 
    end
    
    if ~isfield(options, 'SpecularityIndexed') || options.SpecularityIndexed ~= 0 ...
                || numel(options.SpecularityIndexed) ~= 1
        options.SpecularityIndexed = 1;    %This is the default, i.e. the specularity is 
                                           %not wavelength dependant
    end
    
    if ~isfield(options, 'MaterialEncoding') || (~strcmpi(options.MaterialEncoding, 'RAW') && ...
                ~strcmpi(options.MaterialEncoding, 'GMM'))
        options.MaterialEncoding = 'NURBS';    %Encode the materials using raw measurments
    end
    
    if ~isfield(options, 'IlluminantEncoding') || (~strcmpi(options.IlluminantEncoding, 'RAW') && ...
                ~strcmpi(options.IlluminantEncoding, 'GMM'))
        options.IlluminantEncoding = 'NURBS';    %Encode the illumiants using raw measurments
    end
    
    if ~isfield(options, 'SpecularityEncoding') || (~strcmpi(options.SpecularityEncoding, 'RAW') && ...
                ~strcmpi(options.SpecularityEncoding, 'GMM'))
        options.SpecularityEncoding = 'NURBS';    %Encode the specularities using raw measurments
    end
    
    if ~isfield(options, 'numMaterials') || options.numMaterials>bands || ...
                options.numMaterials<1 || numel(options.numMaterials) ~= 1
        options.numMaterials = min(bands, 5);    %Find five materials per pixel by default
    end
    
    if ~isfield(options, 'maxnumMaterials') || options.maxnumMaterials<options.numMaterials || ...
                options.maxnumMaterials<1 || numel(options.maxnumMaterials) ~= 1
        options.maxnumMaterials = min(25, bands);    %Find five materials per pixel by default
    end
    
    if ~isfield(options, 'numHighlights') || options.numHighlights>bands || ...
                options.numHighlights<1 || numel(options.numHighlights) ~= 1
        options.numHighlights = min(bands, 1);    %Find a single speculairty per pixel by default
    end
    
    if ~isfield(options, 'maxnumHighlights') || options.maxnumHighlights<options.numHighlights || ...
                options.maxnumHighlights<1 || numel(options.maxnumHighlights) ~= 1
        options.maxnumHighlights = min(bands, options.numHighlights);    %Find five materials per pixel by default
    end
    
    if ~isfield(options, 'CanonicalIlluminantIndexed') || numel(options.CanonicalIlluminantIndexed) ~= 1 || ...
            options.CanonicalIlluminantIndexed ~= 1
        if isempty(CanonicalIlluminants)
            options.CanonicalIlluminantIndexed = 0;    %No endmember indexation by default
        else
            options.CanonicalIlluminantIndexed = 1;
        end
    end
    
    if ~isfield(options, 'EndmemberIndexed') ||  numel(options.EndmemberIndexed) ~= 1 || options.EndmemberIndexed ~= 1 
        if isempty(Endmembers)
            options.EndmemberIndexed = 0;    %No endmember indexation by default
        else
            options.EndmemberIndexed = 1;
        end
    end
    
    if ~isfield(options, 'numCanonicalIlluminants')
        if options.CanonicalIlluminantIndexed == 0 ||  numel(options.CanonicalIlluminantIndexed) ~= 1
            options.numCanonicalIlluminants = 0;    %No endmember indexation by default
        elseif options.CanonicalIlluminantIndexed == 1
            options.numCanonicalIlluminants = min(bands, 4);
        end
    elseif options.numCanonicalIlluminants>bands || options.numCanonicalIlluminants<1
        options.numCanonicalIlluminants = min(bands, 4);    %Find four cannonical illuminants per light by default
    end
    
    if ~isfield(options, 'numEndmembers')   
        if numel(options.EndmemberIndexed)~=1 || options.EndmemberIndexed == 0
            options.numEndmembers = 0;    %No endmember indexation by default
        elseif options.EndmemberIndexed == 1 
            options.numEndmembers = min(bands, 4);
        end
    elseif  options.numEndmembers>bands || options.numEndmembers<1
        options.numEndmembers = min(bands, 4);    %Find four endmembers per material
    end
    
    %Setup the number of Gaussians used for the GMM
    if strcmpi(options.MaterialEncoding, 'GMM') && (~isfield(options, 'numGMMsS') || ...
              options.numGMMsS<1 || options.numGMMsS>bands)
    	options.numGMMsS = min(bands, 4);
    end
    
    if strcmpi(options.SpecularityEncoding, 'GMM') && (~isfield(options, 'numGMMsK') || ...
              options.numGMMsK<1 || options.numGMMsK>bands)
        options.numGMMsK = min(bands, 4);
    end
    
    if strcmpi(options.IlluminantEncoding, 'GMM') && (~isfield(options, 'numGMMsL') || ...
              options.numGMMsL<1 || options.numGMMsL>bands) 
        options.numGMMsL = min(bands, 4);
    end
    
    if ~isfield(options, 'IndexingMethod') || ~strcmpi(options.IndexingMethod, 'KM')
         options.IndexingMethod = 'DA';
    end
    
    if ~isfield(options, 'PSFFactor')|| options.PSFFactor<=0 ...
         || numel(options.PSFFactor) ~= 1
         options.PSFFactor = min([10, rows, cols]);
    end
    
    if ~isfield(options, 'DEBUG') || numel(options.DEBUG)~=1 || ...
            options.DEBUG > 5 || options.DEBUG < 0
        options.DEBUG = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recover a common wavelength basis if the HSZ file is to be library
    %indexed. Reduce the image wavelength support if necesary and, if no common
    %bands are found, exit.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Find the common band vector
    wl = I.HDR.wavelength;
    if options.EndmemberIndexed == 1 && ~isempty(Endmembers)
        [wlem, ~] = wavelength_subset(I.HDR.wavelength, Endmembers.HDR.wavelength);
    else
        wlem = wl;
    end
    
    if options.CanonicalIlluminantIndexed == 1 && ~isempty(CanonicalIlluminants)
        [wlci, ~] = wavelength_subset_(I.HDR.wavelength, CanonicalIlluminants.HDR.wavelength);
    else 
        wlci = wl;
    end
    
    if length(wlem)>length(wlci)
        [wl, ~] = wavelength_subset_(wlem, wlci);
    else
        if length(wlem)>length(wlci)
            [wl, ~] = wavelength_subset_(wlci, wlem);
        end
    end
    
    if isempty(wl)
        error('No common bands between the libraries and the image can be found... now exiting\n');
    end
    
    %Adjust the image band span accordingly
    bands = length(wl);
    Bindx = zeros(1, length(wl));
    t=1;
    for i=1:bands
        for j=1:length(wl)
            if I.HDR.wavelength(i)==wl(j) && Bindx(i)==0; %This condintion avoids extra bands in repeated wavelength vectors
                Indx(t)  = i;
                Bindx(i) = 1;
                t = t+1;
            end
        end
    end
    I.I = I.I(:, :, Indx);
    I.HDR.wavelength = wl;
    HSPipeline.HDR.wavelength = wl;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recover the illuminant. For a single illuminant, 
    %the class lookup matrix has a single entry, 
    %with the light abundances being all ones (a single illuminant accounts for
    %the illumination at each pixel). This has been done to make the approach
    %consistent with a spatially varying, canonical illuminant indexed setting.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

    %fprintf('Performing the illuminant power spectrum recovery ...\n');
    if strcmpi(options.method, 'Dichromatic') 
        if options.maxnumIlluminants == 1
            if ~isfield(options, 'L') || ~isstruct(options.L)
                HSPipeline.L.Elements = recover_global_illuminant_(I.I, struct('DEBUG', options.DEBUG));
            end
            if isfield(options, 'L') && isstruct(options.L)
                HSPipeline.L.Elements = recover_global_illuminant_(I.I, options.L);
            end
            %Avoid cases where the illuminant is undefined
            if isnan(HSPipeline.L.Elements)
                HSPipeline.L.Elements = ones(1, bands);
            end
            HSPipeline.L.ElementAbundances = ones(cols, rows, 1);
            HSPipeline.L.ElementAbundanceIndexes = ones(cols, rows, 1);
        else
            if options.maxnumIlluminants > 1
                options.L.ngmm = options.maxnumIlluminants;
                [HSPipeline.L.Elements, ModelL] = recover_multi_illuminant_LG_(I.I, options.L);
                [HSPipeline.L.ElementAbundances, HSPipeline.L.ElementAbundanceIndexes] = sort(ModelL.pLx_prob, 3, 'descend');
                options.maxnumIlluminants = ModelL.q;
                if options.numIlluminants > options.maxnumIlluminants;
                    options.numIlluminants = options.maxnumIlluminants;
                end
                HSPipeline.L.ElementAbundances = HSPipeline.L.ElementAbundances(:, :, 1:options.numIlluminants);
                HSPipeline.L.ElementAbundanceIndexes = HSPipeline.L.ElementAbundanceIndexes(:, :, 1:options.numIlluminants);
                T = sum(HSPipeline.L.ElementAbundances, 3);
                HSPipeline.L.ElementAbundances =  HSPipeline.L.ElementAbundances./...
                     repmat(T, [1 1 options.numIlluminants]);
                clear T;
                clear ModelL;
            else
                 error('The number of illuminants should be greater or equal to one... now exiting\n');
            end
        end
        %Recover the pixel-wise illuminant
        HSPipeline.Light = mix_spectra_(reshape(HSPipeline.L.ElementAbundances, [rows*cols options.numIlluminants]), ...
                                        reshape(HSPipeline.L.ElementAbundanceIndexes, [rows*cols options.numIlluminants]), ...
                                        HSPipeline.L.Elements);
        HSPipeline.Light = reshape(HSPipeline.Light, [cols rows bands]);     
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Set the illuminant variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %If the illuminant is not indexed, transfer the pixelwise illuminant to 
    %HSPipeline.L.Elements    
    if options.IlluminantIndexed == 0
        HSPipeline = rmfield(HSPipeline, 'L');
        HSPipeline.L.Elements = HSPipeline.Light;
    else
        fprintf('Finished with the pixelwise illuminant indexing\n');
    end
    %Add the variables to the header
    HSPipeline.HDR.EncodingL         = options.IlluminantEncoding;
    HSPipeline.HDR.IndexedL          = options.IlluminantIndexed;
    HSPipeline.HDR.EndmemberIndexedL = options.CanonicalIlluminantIndexed;
    HSPipeline.HDR.numElementsL      = options.numIlluminants;
    HSPipeline.HDR.numEndmembersL    = options.numCanonicalIlluminants;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the final, normalised reflectance and scale the shading factor
    %accordingly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('Performing the photometric parameter recovery ...\n');    
    if strcmpi(options.method, 'Dichromatic')
        if ~isfield(options, 'R') || ~isstruct(options.R)
            [HSPipeline.K.Factor, HSPipeline.S.Factor, HSPipeline.K.Elements, HSPipeline.S.Elements] = ...
                recover_dichromatic_parameters_(I.I, HSPipeline.Light, struct('DEBUG', options.DEBUG));
        end
        if isfield(options, 'R') && isstruct(options.R)
            [HSPipeline.K.Factor, HSPipeline.S.Factor, HSPipeline.K.Elements, HSPipeline.S.Elements] = ...
                recover_dichromatic_parameters_(I.I, HSPipeline.Light, options.R);
        end
        %Avoid nans and divisions by zero
        HSPipeline.S.Elements(isnan(HSPipeline.S.Elements)) = 0; 
        HSPipeline.S.Elements(isinf(HSPipeline.S.Elements)) = 0;
        %Setup the other variables to be consistent with a non wavelengthdependant specularity
        options.numHighlights = 1;
    end
    fprintf('done\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Remove unused variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Clear the image to save memory
    clear I;
    %Do likewise with the Light field
    HSPipeline = rmfield(HSPipeline, 'Light'); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Decompose into materials and recover the abundances and class lookup
    %information accordingly. It also clears the raw reflectances at output.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if options.MaterialIndexed==1 || options.EndmemberIndexed==1
        %Recover the Material basis
        Reflectance = HSPipeline.S.Elements;
        [HSPipeline.S.Elements, ~,  EA] = recover_materials(HSPipeline.S.Elements, ...
                                                            options.IndexingMethod, ...
                                                            options.DEBUG, [], ...
                                                            options.maxnumMaterials); 
        [materials, ~] = size(HSPipeline.S.Elements);
        %Check that there are enough highlights in the scene. Update 
        %options.numMaterials otherwise
        options.numMaterials = min(materials, options.numMaterials);
        materials = min(materials, bands);
        fprintf('Performing the indexation to %d materials.\n', options.numMaterials);
        %Get the abundance and indexes matrices
        [~, HSPipeline.S.ElementAbundanceIndexes] = sort(EA, 3, 'descend');
        HSPipeline.S.ElementAbundanceIndexes = HSPipeline.S.ElementAbundanceIndexes(:, :, 1:materials);
        %Get the least squares solution for the most abundant materials
        for i=1:cols
            for j=1:rows
                HSPipeline.S.ElementAbundances(i, j, :)=lsqnonneg(...
                    HSPipeline.S.Elements(reshape(HSPipeline.S.ElementAbundanceIndexes(i, j, :), 1, materials), :)', ...
                    reshape(Reflectance(i, j, :), bands, 1));
            end 
        end
        
        clear Reflectance;
        %Compensate for zero-values
        PSF = fspecial('gaussian', ceil(min(cols, rows)/options.PSFFactor), ...
                       ceil(min(cols, rows)/options.PSFFactor*3));
        for k=1:materials
            Blurred = imfilter(reshape(HSPipeline.S.ElementAbundances(:, :, k), ...
                               [cols rows]), PSF, 'circular', 'conv');
            HSPipeline.S.ElementAbundances(:, :, k) = Blurred;
        end
        
        %Finish with the process
        [HSPipeline.S.ElementAbundances, Indx] = sort(HSPipeline.S.ElementAbundances, 3, 'descend');
        
        eas = reshape(HSPipeline.S.ElementAbundanceIndexes, cols*rows, materials);
        ix2 = reshape(Indx, cols*rows, materials);
        
        ix1 = repmat([1:cols*rows]', [1 materials]);
        ix = sub2ind(size(eas), ix1, ix2);
        HSPipeline.S.ElementAbundanceIndexes = reshape(eas(ix), cols, rows, materials);
        HSPipeline.S.ElementAbundanceIndexes = HSPipeline.S.ElementAbundanceIndexes(:, :, 1:options.numMaterials);
        HSPipeline.S.ElementAbundances       = HSPipeline.S.ElementAbundances(:, :, 1:options.numMaterials);
     end
    %Add the variables to the header
    HSPipeline.HDR.EncodingS         = options.MaterialEncoding;
    HSPipeline.HDR.IndexedS          = options.MaterialIndexed;
    HSPipeline.HDR.EndmemberIndexedS = options.EndmemberIndexed;
    HSPipeline.HDR.numElementsS      = options.numMaterials;
    HSPipeline.HDR.numEndmembersS    = options.numEndmembers;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup the variables for the specularities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if options.SpecularityIndexed == 1 
        fprintf('Performing the specularity indexation ...\n');
        if options.numHighlights == 1
            HSPipeline.K.Elements = ones(1, bands);
            HSPipeline.K.ElementAbundanceIndexes = ones(cols, rows, 1);
            HSPipeline.K.ElementAbundances = ones(cols, rows, 1);
        else
        %Recover the highlight basis
        [HSPipeline.K.Elements, ~, EA] = recover_materials_(HSPipeline.K.Elements, ...
            options.IndexingMethod, options.DEBUG, [], options.maxnumHighlights); 
        [materials, ~] = size(HSPipeline.K.Elements);
        %Check that there are enough highlights in the scene. Update 
        %options.numHighlights otherwise
        options.numHighlights = min(materials, options.numHighlights);
        materials = min(materials, bands);
        fprintf('Performing the indexation to %d highlights.\n', options.numHighlights);
        %Get the abundance and indexes matrices
        [~, HSPipeline.K.ElementAbundanceIndexes] = sort(EA, 3, 'descend');
        HSPipeline.K.ElementAbundanceIndexes = HSPipeline.K.ElementAbundanceIndexes(:, :, 1:materials);
        %Get the least squares solution for the most abundant highlights
        for i=1:cols
            for j=1:rows
                HSPipeline.K.ElementAbundances(i, j, :)=lsqnonneg(...
                    HSPipeline.K.Elements(reshape(HSPipeline.K.ElementAbundanceIndexes(i, j, :), 1, materials), :)', ...
                    reshape(HSPipeline.K.Elements(i, j, :), bands, 1));
            end 
        end
        %Compensate for zero-values
        PSF = fspecial('gaussian', ceil(min(cols, rows)/options.PSFFactor), ...
                       ceil(min(cols, rows)/options.PSFFactor*3));
        for k = 1:materials
            Blurred = imfilter(reshape(HSPipeline.K.ElementAbundances(:, :, k), ...
                               [rows cols]), PSF, 'circular', 'conv');
            HSPipeline.K.ElementAbundances(:, :, k) = Blurred;
        end
        %Finish with the process
        [~, Indx] = sort(HSPipeline.K.ElementAbundances, 3, 'descend');
        HSPipeline.K.ElementAbundanceIndexes = HSPipeline.K.ElementAbundanceIndexes(:, :, Indx(1:options.numHighlights));
        HSPipeline.K.ElementAbundances = HSPipeline.K.ElementAbundances(:, :, Indx(1:options.numHighlights));
        fprintf('done\n');
        end
    end
    %Add the variables to the header
    HSPipeline.HDR.numElementsK = options.numHighlights;
    HSPipeline.HDR.EncodingK    = options.SpecularityEncoding;
    HSPipeline.HDR.IndexedK     = options.SpecularityIndexed;
    HSPipeline.HDR.numElementsK = options.numHighlights;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the library indexing if required
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the reflectance
    if (options.EndmemberIndexed==1)
        HSPipeline.HDR.numEndmembersS = min([options.numEndmembers bands Endmembers.HDR.numEndmembers]);
        [HSPipeline.S.EndmemberAbundances, HSPipeline.S.EndmemberAbundanceIndexes, HSPipeline.S.Endmembers, wavelength] = ...
            L2_unmixing_(HSPipeline.S.Elements, HSPipeline.HDR.wavelength, Endmembers.Endmembers, ...
            Endmembers.HDR.wavelength, HSPipeline.HDR.numEndmembersS);
        %Setup the header values
        HSPipeline.HDREndmembersS = Endmembers.HDR;
        HSPipeline.HDREndmembersS.wavelength=wavelength;
    end
    
    %Do the illuminants
    if (options.CanonicalIlluminantIndexed==1)
        HSPipeline.HDR.numEndmembersL = min([options.numCanonicalIlluminants bands CanonicalIlluminants.HDR.numEndmembers]);
        [HSPipeline.L.EndmemberAbundances, HSPipeline.L.EndmemberAbundanceIndexes, HSPipeline.L.Endmembers, wavelength] = ...
            L2_unmixing_(HSPipeline.L.Elements, HSPipeline.HDR.wavelength, CanonicalIlluminants.Endmembers, ...
            CanonicalIlluminants.HDR.wavelength, HSPipeline.HDR.numEndmembersL);
        %Setup the header values
        HSPipeline.HDREndmembersL = CanonicalIlluminants.HDR;
        HSPipeline.HDREndmembersL.wavelength=wavelength;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the encoding for non-RAW formats
    %Note that, for the NURBS, the last row of the Elements/Endmember matrix
    %correponds to the control points on the wavelength space.
    %For the GMMs, the ordering is weights, means and standard deviations.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HSPipeline = encode_HSZ_(HSPipeline, options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Remove the raw reflectance and the specularities before exiting.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    HSPipeline.HDR.bands = length(HSPipeline.HDR.wavelength);
end
