%% Recover the common set of wavelength values from two wavelength vectors
%
%% Syntax:
%     [MATCHED, IDX] = wavelength_subset(A, B)
%
%% Description:
%     Return the indices and values for the subset of wavelengths in
%     A that are within the range of wavelengths in B, i.e. bound by B.
%
%% Input:
%     A, B: Vectors of wavelength values
%
%% Output:
%     IDX: Array of wavelength indeces for the values in A which are within the range in B. If
%        there are no values of A in the range of B, IDX is an empty array
%     MATCHED: Subset of wavelength values in A which are within the range in B. If
%        there are no values of A in the range of B, MATCHED is an empty array
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.6
% Last Update Date: 25 July 2014

function [MATCHED, IDX] = wavelength_subset_(A, B)

    if ndims(A) ~= 2 || ndims(B) ~= 2
        error('Wavelength vectors must be two-dimensional arrays');
    end
    
    IDX = find(A >= min(B) & A <= max(B));
    
    if ~isempty(IDX)
        MATCHED = A(IDX);
        MATCHED = reshape(MATCHED, [1 length(MATCHED)]);
    else
        MATCHED = [];
    end
    
end