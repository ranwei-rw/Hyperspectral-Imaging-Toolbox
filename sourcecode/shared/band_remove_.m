%   Description:
%       this is a function used removed unwanted bands from given
%       hyperspectral image data cube. This function will estimate the
%       illuminant power spectrum using given data and remove all bands
%       that have a lower or equal illuminant value than the given threshold. 
%
%   [ND, Ind] = band_remove_(D, threshold, bottom)
%  
%   Input: 
%   
%       D:         Hyperspectral image data structure or cube. If a HDR
%                  structure is contained in D, this HDR structure will
%                  also be modified accordingly.
%       threshold: the value of illuminant below which bands will be
%                  removed. 0 is for a band being removed and 1 is for a
%                  band kept.
%       Bottom:    indicates if this function is going to remove dark bands
%                  (<= threshold) or saturated bands (> threshold). By
%                  default it's to remove dark bands.
%
%   Output: 
%   
%       NI:        image data or structure after bands removed
%       Ind:       the indexes of bands removed
%
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.1
% Last Update Date: 13 Oct 2014
%


function [NI, Ind] = band_remove_(I, threshold, bottom)

    if ~exist('bottom', 'var')
       bottom = 1;
    end
    
    % check the structure of input data D
    has_struct = 0;
    has_header = 0;
    
    if isfield(I, 'I')
        has_struct = 1;
    end
    
    if isfield(I, 'HDR')
        has_header = 1;
    end
    %   copy origin data
    if has_struct
        OD = I.I;
    else
        OD = I;
    end
    
    [~, ~, b] = size(OD);
    if b <= 1
        clear OD;
        error('this function is designed to process hyperspectral image cube');
    end
    
    %   Now to estimate the illuminants using toolbox
    L = recover_global_illuminant(OD);
    
    Ind = zeros(size(L));
    
    if bottom
        Ind = L > threshold;
    else
        Ind = L <= threshold;
    end
    
    OD = OD(:, :, Ind > 0);
    
    if has_header
        if isfield(I.HDR, 'wavelength')
            I.HDR.wavelength  = I.HDR.wavelength(Ind > 0);
        end
        if isfield(I.HDR, 'bands')
            I.HDR.bands = sum(Ind);
        end
    end
    
    if has_struct
        NI.I = OD;
    else
        NI = OD;
    end
    
    if has_header
        NI.HDR = I.HDR;
    end
    
%   end of function    
end