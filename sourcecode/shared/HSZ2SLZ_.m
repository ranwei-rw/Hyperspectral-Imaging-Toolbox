%% Create a library from an indexed HSZ data structure
%
%% Syntax
%
%   SLZ = HSZ2SLZ(HSZ, MaterialList, field);
%   SLZ = HSZ2SLZ(HSZ, [], field);
%
%% Input
%
%   HSZ:          Scyllarus data structure.
%   MaterialList: Cell array containing the names of the materials or lights to be
%                 stored on the library (SLZ)
%   field:        Determines whether the materials or lights are to be stored. If
%                 they are to be saved, field = 'S', otherwise field = 'L'.
%
%% Output
%
%   SLZ: Scyllarus library data structure
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.1.0
% Last Update Date: 23 July 2014


function SLZ = HSZ2SLZ_(HSZ, MaterialList, field)

    if nargin < 3
        error('Not enough input arguments');
    end
    
    switch upper(field)
        case 'L'
            if HSZ.HDR.IndexedL ~= 1
                error('The HSZ data structure L is not indexed');
            end
            switch upper(HSZ.HDR.EncodingL)
                case 'NURBS'
                    [mats, cols] = size(HSZ.L.Elements);
                    Q = eval_nurbs(HSZ.L.ElementKnots, HSZ.HDR.wavelength, ...
                        reshape(HSZ.L.Elements, [mats 1 cols]),  HSZ.L.ElementCP',  HSZ.HDR.degreeNURBSL); 
                    SLZ.Elements = reshape(Q, [mats length(HSZ.HDR.wavelength)]);
                case 'GMM'
                    [mats, ~]=size(HSZ.L.Elements);
                    %clear HSZ.S.Elements;
                    Q = eval_gaussian_mixture_(HSZ.L.Elements, ...
                        HSZ.L.ElementMean, HSZ.L.ElementStd, HSZ.HDR.wavelength);
                    SLZ.Elements = reshape(Q, [mats length(HSZ.HDR.wavelength)]);
                case 'RAW'
                    SLZ.Elements = HSZ.L.Elements;
                otherwise
                    error('Unknown encoding method for L');
            end
            SLZ.HDR.wavelength = HSZ.HDR.wavelength;
            if isfield(HSZ.HDR, 'wavelength_unit')
                SLZ.HDR.wavelength_unit = HSZ.HDR.wavelength_unit;
            end
            [SLZ.HDR.numEndmembers, ~] = size(HSZ.L.Elements);
            if ~isempty(MaterialList) && length(MaterialList) == SLZ.HDR.numEndmembers
                for i=1:SLZ.HDR.numEndmembers
                    SLZ.HDR.(strcat('MAT', int2str(i))) = MaterialList{i};
                end
            else
                for i=1:SLZ.HDR.numEndmembers
                    SLZ.HDR.(strcat('MAT', int2str(i))) = strcat('Material', int2str(i));
                end
            end
        case 'S'
            if HSZ.HDR.IndexedS ~= 1
                error('The HSZ data structure S is not indexed');
            end
            switch upper(HSZ.HDR.EncodingS)
                case 'NURBS'
                    [mats, cols] = size(HSZ.S.Elements);
                    Q = eval_nurbs(HSZ.S.ElementKnots, HSZ.HDR.wavelength, ...
                        reshape(HSZ.S.Elements, [mats 1 cols]),  HSZ.S.ElementCP',  HSZ.HDR.degreeNURBSS); 
                    SLZ.Elements = reshape(Q, [mats length(HSZ.HDR.wavelength)]);
                case 'GMM'
                    [mats, ~]=size(HSZ.S.Elements);
                    %clear HSZ.S.Elements;
                    Q = eval_gaussian_mixture_(HSZ.S.Elements, ...
                        HSZ.S.ElementMean, HSZ.S.ElementStd, HSZ.HDR.wavelength);
                    SLZ.Elements = reshape(Q, [mats length(HSZ.HDR.wavelength)]);
                case 'RAW'
                    SLZ.Elements = HSZ.S.Elements;
                otherwise
                    error('Unknown encoding method for S');
            end       
            SLZ.HDR.wavelength = HSZ.HDR.wavelength;
            if isfield(HSZ.HDR, 'wavelength_unit')
                SLZ.HDR.wavelength_unit = HSZ.HDR.wavelength_unit;
            end

            [SLZ.HDR.numEndmembers, ~] = size(HSZ.S.Elements);
            if ~isempty(MaterialList) && length(MaterialList) == SLZ.HDR.numEndmembers
                for i=1:SLZ.HDR.numEndmembers
                    SLZ.HDR.(strcat('MAT', int2str(i))) = MaterialList{i};
                end
            else
                for i=1:SLZ.HDR.numEndmembers
                    SLZ.HDR.(strcat('MAT', int2str(i))) = strcat('Material', int2str(i));
                end
            end
        otherwise
            error('The HSZ data structure is not indexed or the field input string specified is not supported');
    end
end

        
    


