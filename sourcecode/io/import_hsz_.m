% Usage:
%     HSZ = import_hsz_(filename);
% 
% Writes a compressed NICTA HSZ file. This is an HDF5 file with all
% the variables recovered by the HSZ routine.
%
% Input:
%     filename: The name of the file (including the path) to be read from
%     disk.
%
% Output:
%     HSZ: Structure containing the data delivered by the NICTA pipeline. 
%       This is the same as the output of the NICTAPipeline routine. 
%             
% Example:

%      HSZ = import_hsz_('test.hsz');
%
% See also:
%   NICTAPipeline, export_hsz_
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.7
% Last Update Date: 7 Aug 2014

function HSZ = import_hsz_(filename)

    %   check whether the file exists
    if ~isequal(exist(filename, 'file'), 2)
        error('%s does not exist. Please check.', filename);
    end

    file_fp = fopen(filename);
    if file_fp == -1
        %   when it fails to open the given file, prompt error message
        error('Failed to open file %s', filename);
    end
    fclose(file_fp);

    Info = h5info(filename);
    %Go through the groups
    for i = 1:length(Info.Groups)
        %Start with the header. In reallity, this will search for each of the
        %groups and will allocate structures accordingly
        fieldname = strtok(Info.Groups(i).Name, '/');
        %Do the datasets first
        for j = 1:length(Info.Groups(i).Groups)
            %varname = strtok(Info.Groups(i).Groups(j).Name, strcat('/', fieldname));
            varname = Info.Groups(i).Groups(j).Name(length(strcat('/', fieldname))+2:end);
            DATA    = h5read(filename, strcat(Info.Groups(i).Groups(j).Name, '/DATA'));
            max_     = h5read(filename, strcat(Info.Groups(i).Groups(j).Name, '/MAX'));
            min_     = h5read(filename, strcat(Info.Groups(i).Groups(j).Name, '/MIN'));
            %Do the scaling for data saved as integers
            if isinteger(DATA)
                %Recover the upper bound for the data if it was saved as an integer
                 m = max(DATA);
                 n = ndims(DATA);
                 for k = 2:n
                    m = max(m);
                 end
                %Add a case for those variables that are zero, specially in the header
                if (max_-min_) ~= 0
                    %Do not scale the arrays containing indexes
                    if ~isempty(strfind(varname, 'Indexes'))
                        HSZ.(fieldname).(varname) = double(DATA);
                    else
                        HSZ.(fieldname).(varname) = double(DATA)/double(m)*(max_-min_)+min_;
                    end
                else
                    HSZ.(fieldname).(varname) = double(DATA);
                end
            end
        end
        for j = 1:length(Info.Groups(i).Attributes)
            HSZ.(fieldname).(Info.Groups(i).Attributes(j).Name) = h5readatt(filename, ...
                Info.Groups(i).Name, Info.Groups(i).Attributes(j).Name);
        end
    end

    if isfield(HSZ.HDR, 'wavelength')
        HSZ.HDR.wavelength = uint16(HSZ.HDR.wavelength);
    end

    if isfield(HSZ.HDR, 'wavelength_unit') && ...
       strcmpi(HSZ.HDR.wavelength_unit, 'micrometers')
        HSZ.HDR.wavelength_unit = 'Nanometers';
        HSZ.HDR.wavelength = HSZ.HDR.wavelength*1000;
    end
end
