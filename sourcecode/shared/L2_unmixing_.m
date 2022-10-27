% L2 unmixing routine for spectral data
%
% Syntax:
%       [UNMIX_WEIGHT_MAP, UNMIX_INDEX_MAP, ENDMEMBER, WAVELENGTH] = L2_unmixing(S, S_WL, BASIS, BASIS_WL, num_endmember)
%
% Description:
%       This function is designed to unmix spectral images into combinations of reference spectra
%       or end-members. This routine uses NURBS to perform unmixing of spectra with respect to 
%       end-members with different spectral resolutions. Note that the value of input S and BASIS will be
%       squared-normalised to 1 per pixel.
%
% Input:
%
%     S:        Hyperspectral image data cube. This can be either a 3-D array 
%               whose third dimension is the wavelength domain or a 2D array 
%               whose second dimension is indexed to wavelength
%     S_WL:     Vector of wavelengths used for the evaluation of the spectra in S.
%     BASIS:    The input library of endmember spectra. This is a 2D array whose
%               second dimension is the wavelength domain.
%     BASIS_WL: Vector of wavelengths used for the evaluation of the end-member spectra.
%     num_endmember:   Number of end-members used for the unmixing operation.
%
% Output:
%
%     UNMIX_WEIGHT_MAP: The proportion of the end-member corresponding to every pixel in the input reflectance image,
%                       the index of BASIS is given in UNMIX_INDEX_MAP. It will be nomralised to 1.
%     UNMIX_INDEX_MAP:  The index of the end-member for UNMIX_WEIGHT_MAP
%     ENDMEMBER:        Matrix of endmembers evaluated over the common wavelength
%                       range for both, the image data and the library.
%     WAVELENGTH:       Vector of common wavelength values for the image data and
%                       the end-member library.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.11
% Last Update Date: 28 Aug 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [UNMIX_WEIGHT_MAP, UNMIX_INDEX_MAP, ENDMEMBER, WAVELENGTH] = L2_unmixing_(S, S_WL, BASIS, BASIS_WL, num_endmember)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Setup the variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BASIS_nd = ndims(BASIS);

    if BASIS_nd == 3
        [BASIS_rowsold, BASIS_colsold, BASIS_bands] = size(BASIS);
        BASIS = reshape(BASIS, [BASIS_rowsold*BASIS_colsold 1 BASIS_bands]);
    else
        [BASIS_rowsold, BASIS_bands] = size(BASIS);
        BASIS = reshape(BASIS, [BASIS_rowsold 1 BASIS_bands]);
    end
    
    I_nd = ndims(S);
    if I_nd == 2
        [I_rowsold, I_colsold] = size(S);
        S = reshape(S, [I_rowsold 1 I_colsold]);
    end
    [rows, cols, ~] = size(S);
    [n, ~, BASIS_bands] = size(BASIS);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Assure the wavelengths are compatible and given by row vectors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % make sure all wavelengths are in vector format
    S_WL = reshape(S_WL, length(S_WL), 1);
    BASIS_WL = reshape(BASIS_WL, length(BASIS_WL), 1);
    use_orig_wavelength = 0;
    common_index = [];
    if ~exist('BASIS_WL', 'var')
        BASIS_WL = S_WL;
        wl = S_WL;
        use_orig_wavelength = 1;
        if isempty(wl)
            error('Spectra wavelength is empty');
        end
    else
        [wl, common_index] = wavelength_subset_(S_WL, BASIS_WL);
        if isempty(wl)
            error('There are no common wavelengths between the library and the spectra');
        end
    end    

    if ~exist('num_endmember', 'var')
        num_endmember = length(wl);
    end

    bands = length(wl);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Start the evaluation. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isequal(BASIS_WL, S_WL)
        [BASIS_KNOTS, BASIS_CP_REF, BASIS_CP_WAVE] = get_nurbs_(BASIS, reshape(BASIS_WL, [BASIS_bands 1]), 2);
        Q = eval_nurbs_(BASIS_KNOTS, reshape(wl, [bands 1]), BASIS_CP_REF, BASIS_CP_WAVE, 2);   
    else
        Q = BASIS;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Recover the spectra in an equal wavelength basis to that in wl
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if use_orig_wavelength == 0
        S     = S(:, :, common_index);
        bands = length(common_index);
    end

    %   to make values in the same scale, we here normalise S and Q. 
    S           = S ./ repmat(sqrt(sum(S.^2, 3)), [1 1 bands]);
    S(isnan(S)) = 0;
    Q           = Q ./ repmat(sqrt(sum(Q.^2, 3)), [1 1 bands]);
    Q(isnan(Q)) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Do the actual unmixing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   Only use those materials that are more "similar" to each spectrum
    D         = dist_(reshape(Q, [n bands]), reshape(S, [rows*cols bands])');
    xm        = min(min(D(D>0)));  %   Avoid divisions by zero
    D(D==0)   = xm;
    Scale     = sum(D, 2);
    EA        = D./Scale(:, ones(rows*cols, 1));
    materials = min(num_endmember, bands);
    
    s = sprintf('Performing the unmixing to %d end members.', materials);
    disp(s);
    %   Get the abundance and indexes matrices
    EA           = reshape(EA', [rows cols n]);
    [~, EA_Indx] = sort(EA, 3, 'ascend');
    EA_Indx      = EA_Indx(:, :, 1:materials);
    
    UNMIX_WEIGHT_MAP = zeros(rows, cols, materials);
    %   Get the least squares solution for the most abundant highlights
    %   set options. 
    optionsoptim = optimset('MaxIter', 100);
    for i = 1:rows
        for j = 1:cols
            %   get normalised initial values for x
            id = reshape(EA_Indx(i, j, :), [1, materials]);
            CC = reshape(Q(id, 1, :), [materials bands])';
            dd = reshape(S(i, j, :), [bands, 1]);
            UNMIX_WEIGHT_MAP(i, j, :) = lsqnonneg(CC, dd, optionsoptim);
        end 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Reshape and exit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if I_nd == 2
        UNMIX_WEIGHT_MAP = reshape(UNMIX_WEIGHT_MAP, [I_rowsold materials]);
        temp = sum(UNMIX_WEIGHT_MAP, 2);
        temp(temp == 0) = 1;
        for m = 1:materials
            UNMIX_WEIGHT_MAP(:, m) = UNMIX_WEIGHT_MAP(:, m)./temp;
        end
        UNMIX_INDEX_MAP  = reshape(EA_Indx, [I_rowsold materials]);
    else
        UNMIX_INDEX_MAP  = EA_Indx;
        temp = sum(UNMIX_WEIGHT_MAP, 3);
        temp(temp == 0) = 1;
        for m = 1:materials
            UNMIX_WEIGHT_MAP(:, :, m) = UNMIX_WEIGHT_MAP(:, :, m) ./ temp;
        end
    end

    ENDMEMBER = reshape(Q, [n bands]);
    WAVELENGTH = wl;

%   end of function L2_unmixing_
end