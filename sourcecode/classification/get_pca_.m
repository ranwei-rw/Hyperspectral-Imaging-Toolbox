% Use principal Components Analysis (PCA) to a hyperspectral image file
%
% Syntax:
%      D = get_pca(I, n)
%
% Description:
%   This function is designed to conduct PCA on given a hyperspectral file.
%   It will return a data structure (the same as input I) which contains
%   the first n principal components generated from I. 
%
%   Inputs:
%
%       I: Structure containing the image cube and the header. Or, a
%          hyperspectral image cube of size height-by-width-by-band.
%       n: Number of principal components to be kept of PCA
%
%   Output:
%
%       D: The data structure containing n principal components. It shares
%          the same structure as I. D will has the same value range as I
%          does.
%
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.0
% Last Update Date: 15 Sep 2014

function D = get_pca_(I, n)

    %   check the data structure;
    has_struct = 0;
    has_header = 0;
    if ~exist('n', 'var')
        n = 3;
    end
    
    if isfield(I, 'I')
        has_struct = 1;
    end
    
    if isfield(I, 'HDR')
        has_header = 1;
    end
    %   copy origin data - OD (Original Data)
    if has_struct
        OD = I.I;
    else
        OD = I;
    end
    
    [height, width, bands] = size(OD);
    if bands <= 1
        clear OD;
        error('this function is designed to process hyperspectral image cube');
    end
    %   reshape to 2D format - entry number by bands
    OD     = reshape(OD, [height*width, bands]);
    OD_M   = mean(OD);      %   get its mean
    od_max = max(OD(:));       %   get the max value of OD
    
    %   substract mean from data
    for i = 1:bands
        OD(:, i) = OD(:, i) - OD_M(i);
    end
     
    OD        = OD/od_max; %   get unitilised
	C         = cov(OD);   %   get covariance
    [V, EIGD] = eig(C);    %   get eigenvectors and eigenvalues
    MAXD      = max(EIGD);    %   sort eigenvalues descendingly
    [~, IND]  = sort(MAXD, 'descend');

    
    
    if has_header
        D.HDR       = I.HDR;
        D.HDR.bands = n;
        if isfield(I.HDR, 'wavelength')
            D.HDR.wavelength = 0;
        end
    end
    
    VV = zeros(bands, bands);%   sort eigenvectors accordingly
    if n > bands
        warning('Scyllarus:shared:get_pca', 'Requred number of Principal Components is more than band number');
    end
    for i = 1:n
        VV(:, i) = V(:, IND(i));
    end


    VVS = VV(:, 1:n);   %   get the first n PCs

    COVA = (VVS'*OD')';

    PCA = reshape(COVA, [height width n])*od_max;
    
    if has_struct
        D.I = PCA;
    else
        D = PCA;
    end
    
    clear COVA;
    clear VVS;
    clear PCA;
    clear OD;
    
end