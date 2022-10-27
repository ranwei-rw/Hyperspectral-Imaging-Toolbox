% Mix end-members so as to obtain mixed spectra.
%
% Syntax
%     Q = mix_spectra(ABUNDANCES, INDEXES, ENDMEMBERS)
%     Q = mix_spectra(ABUNDANCES, INDEXES, ENDMEMBERS, material_num)
%     Q = mix_spectra(ABUNDANCES, INDEXES, ENDMEMBERS, material_num, normalise)
% 
% Description
%
%   Mix the spectra given a set of abundances, indexes and endmembers. 
%   This is, effectively, the inverse operation to unmixing.
%
% Inputs:
%
%     ABUNDANCES:   Matrix containign the abundance coefficients per endmember. This is a 
%                   cols x rows x endmembers matrix, where the output matrix Q is cols x rows x bands.
%     INDEXES:      Matrix containing the indexes for the endmembers.
%     ENDMEMBERS:   Matrix whose rows are indexed to the endmembers and columns to the 
%                   wavelength, i.e. Endmember(i, :) contains the ith endmember spectrum
%     material_num: Number of materials used for the recovery of the spectra. If 
%                   this is not provided, the size of 3rd dimension of the abundance matrix is used.
%     normalise:    If normalise is unity, the output spectra is normalised accordingly, 
%                    i.e. norm(Q(i, i), 3) is equivalent to 1. The default is 1.
% 
% Outputs:
%    Q: Matrix of mixed spectra
% 
% See also
%   L2_unmixing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2015 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly.
% Version: 1.1.0
% Last Update Date: 23 July 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Q = mix_spectra_(ABUNDANCES, INDEXES, ENDMEMBERS, material_num, normalise)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Check inputs and set the default to a band-normalised spectra matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 3
        error('Not enough input arguments.');
    end
    
    if ~exist('normalise', 'var') || normalise ~=  0
        normalise = 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup the parameters before starting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~, bands] = size(ENDMEMBERS);
    
    isnot_3d = ndims(ABUNDANCES)<3;
    if isnot_3d
        %   not a 3D matrix
        [rows, weights] = size(ABUNDANCES);
        Q = zeros(rows, bands);
    else
        [rows, cols, weights] = size(ABUNDANCES);
        Q = zeros(rows, cols, bands);
    end

    if ~exist('material_num', 'var')
        material_num = weights;
    else
        material_num = min(weights, material_num);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the mixing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isnot_3d
        for r = 1:rows
            if material_num > 1
                for w = 1:material_num
                    Q(r, :) = Q(r, :)+reshape(ABUNDANCES(r, w)*ENDMEMBERS(INDEXES(r, w), :), 1, bands);
                end
            else
                Q(r, :) = reshape(ABUNDANCES(r)*ENDMEMBERS(INDEXES(r), :), 1, bands);
            end
        end
    else
        for r = 1:rows
            for c = 1:cols
                if material_num > 1
                    for w = 1:material_num
                        Q(r, c, :) = Q(r, c, :)+reshape(ABUNDANCES(r, c, w)*ENDMEMBERS(INDEXES(r, c, w), :), 1, 1, bands);
                    end
                else
                    Q(r, c, :) = reshape(ABUNDANCES(r, c)*ENDMEMBERS(INDEXES(r, c), :), 1, 1, bands);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Normalise the matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if normalise ==  1
        if isnot_3d
            Norm = sqrt(sum(Q .* Q, 2));
            Norm(Norm == 0) = 1;
            Q = Q ./ Norm(:, ones(bands, 1));
        else
            Norm = sqrt(sum(Q .* Q, 3));
            Norm(Norm == 0) = 1;
            Q = Q ./ Norm(:, :, ones(bands, 1));
        end
    end

%   end of function mix_spectra_
end
