%   Syntax:
%       function SLZ = mat2slz(MAT, HDR, normalise)
%
%   Description:
%       this function is used to convert a matrix into a SLZ library 
%
%   Input:
%       MAT: matrix to be converted. size: numEndmember column vectors
%           (each of them has bands elements)
%       Labels: is a string vector containing the names of each library
%               vector. By default, it values will be 'mat1', 'mat2', ...,
%               'matn' where n is the number of library vectors. It should
%               contain the same number of elements as that of library
%               vectors.
%       HDR: header information required. It should contains the following
%            information (but not limited to): wavelength, wavelength_unit
%            and numEndmembers. If HDR is missing, this function will
%            use default values such as wavelength = [number_of_bands, 1]
%            with value starting from 400 to 700 and increase by 10.
%            wavelength unit 'Nanometers'. numEndmembers will be the number
%            of library vectors.
%       normalise: whether normalise input value. Default value is 1 (yes).
%
%   Output
%       SLZ: result in SLZ format

function SLZ = mat2slz_(MAT, Labels, HDR, normalise)
    
    if ~exist('MAT', 'var')
        error('Missing input arguments');
    end
    
    [bands, num] = size(MAT);
    
    if ~exist('normalise', 'var')
        normalise = 1;
    end
    
    if ~exist('HDR', 'var')
        HDR.wavelength_unit = 'Nanometers';
        HDR.wavelength = [400:10:700]';
        HDR.numEndmembers = num;
    else
        if ~isfield(HDR, 'wavelength_unit')
            HDR.wavelength_unit = 'Nanometers';
        end
        
        if ~isfield(HDR, 'wavelength')
            HDR.wavelength = [400:10:700]';
        end
        
        if ~isfield(HDR, 'numEndmembers')
            HDR.numEndmembers = num;
        end
    end
    
    if ~exist('Labels', 'var') || isempty(Labels)
        for i = 1:num
            Labels{i} = sprintf('mat%d', i);
        end
    end
    
    %   now to compose the variable
    if normalise == 1
        MAT = MAT/max(MAT(:));
    end
    
    SLZ.Endmembers = MAT';
    for i = 1:num
        SLZ.Endmembers((i-1)*bands +1:(i-1)*bands+bands) = MAT(:, i);
    end
    SLZ.HDR = HDR;
    SLZ.Labels = char(Labels);
    
end