% Usage:
%     HSPipeline = import_slz_(filename);
% 
% Reads a compressed NICTA HSPipeline library file. This is an HDF5 file.
%
% Input:
%     filename: The name of the file (including the path) to be read from
%     disk.
%
% Output:
%     HSPipeline: Structure containing the library data. 
%             
% Example:
%      HSPipeline = SLZread('test.hsz');
%
% See also:
%   NICTAPipeline, HSZwrite, SLZwrite, HSZread
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2015 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.7
% Date: 05 Feb 2015
% Added support to polygon labels generated by Scyven 
%
% Version: 1.0.6
% Last Update Date: 7 Aug 2014

function HSPipeline = import_slz_(filename)

    %   check file name & path
    if isequal('\', filesep)
        % most systems use \ as separator
        [~, pos] = strtok(filename, '/');
        if ~isempty(pos)
            filename = strrep(filename, '/', filesep);
        end
    else
        [~, pos] = strtok(filename, '\');
        if ~isempty(pos)
            filename = strrep(filename, '\', filesep);
        end
    end

    if ~exist(filename, 'file')
        s = strcat('File does not exist - ', filename);
        error(s);
    end
    
    Info = h5info(filename);
    old_slz_format = 0;
    %   go through groups contained in the file
    for i = 1:length(Info.Groups)
        %   Start with the header. In reallity, this will search for each of the
        %   groups and will allocate structures accordingly
        fieldname = strtok(Info.Groups(i).Name, '/');
        %   Do the datasets first
        for j = 1:length(Info.Groups(i).Groups)
            varname = strtok(Info.Groups(i).Groups(j).Name, strcat('/', fieldname));
            if strcmpi(varname, 'GB')
                varname = 'RGB';
            end
            DATA = h5read(filename, strcat(Info.Groups(i).Groups(j).Name, '/DATA'));
            MAX  = h5read(filename, strcat(Info.Groups(i).Groups(j).Name, '/MAX'));
            MIN  = h5read(filename, strcat(Info.Groups(i).Groups(j).Name, '/MIN'));
            %Do the scaling for data saved as integers
            if isinteger(DATA)
                %Recover the upper bound for the data if it was saved as an integer
                m = max(DATA);
                n = ndims(DATA);
                for k = 2:n
                    m = max(m);
                end
                %Add a case for those variables that are zero, specially in the header
                if MAX-MIN ~= 0
                    %Do not scale the arrays containing indexes
                    if ~isempty(strfind(varname, 'Indexes'))
                        HSPipeline.(fieldname).(varname) = double(DATA);
                    else
                        HSPipeline.(fieldname).(varname) = double(DATA)/double(m)*(MAX-MIN)+MIN;
                    end
                else
                    HSPipeline.(fieldname).(varname) = double(DATA);
                end
            elseif isfloat(DATA)
                HSPipeline.(fieldname).(varname) = double(DATA);
            else
                error('Unsupported data type %s in %s', varname, fieldname);
            end
        end
        %   Read the datasets
        %   get the number of datasets in this Group
        n = length(Info.Groups(i).Datasets);
        if mod(n, 3) == 0 && ~strcmp(Info.Groups(i).Name,'/HDR')
            %   this is a group contains DATA, MAX and MIN datasets
            for j = 1:3:n
                varname = fieldname;
                if strcmpi(varname, 'GB')
                    varname = 'RGB';
                end
                DATA = h5read(filename, strcat(Info.Groups(i).Name, '/DATA'));
                MAX  = h5read(filename, strcat(Info.Groups(i).Name, '/MAX'));
                MIN  = h5read(filename, strcat(Info.Groups(i).Name, '/MIN'));
                %Do the scaling for data saved as integers
                if isinteger(DATA)
                    %Recover the upper bound for the data if it was saved as an integer
                    m = max(DATA);
                    n = ndims(DATA);
                    for k = 2:n
                        m = max(m);
                    end
                    %Add a case for those variables that are zero, specially in the header
                    if (MAX-MIN) ~= 0
                        %Do not scale the arrays containing indexes
                        if ~isempty(strfind(varname, 'Indexes'))
                            HSPipeline.(varname) = double(DATA);
                        else
                            HSPipeline.(varname) = double(DATA)*(MAX-MIN)/double(m)+MIN;
                        end
                    else
                        HSPipeline.(varname) = double(DATA);
                    end
                elseif isfloat(DATA)
                    HSPipeline.(varname) = double(DATA);
                else
                    error('Unsupported data type in %s', varname);
                end
            end
        else
            if isfield(Info.Groups(i),'Labels') %Check for the labels here
                for j = 1:n
                    varname = fieldname;
                    LABELS = char(h5read(filename, strcat(Info.Groups(i).Name, '/DATA')));

                    if ischar(LABELS)
                        HSPipeline.(varname) = LABELS;
                    else
                        error('Unsupported data type in %s', varname);
                    end
                end
            end
        end
        pos = [];
        labelnames = [];
        
        %Do the attributes
        for j = 1:length(Info.Groups(i).Attributes)
            s = Info.Groups(i).Attributes(j).Name;

            if beginswith(s, 'mat')
                if old_slz_format == 0
                    old_slz_format = 1;
                end
                %   then this is a material label stored in old slz format
                %   save the names and values into two vectors, one is
                %   number and the other one is cell array
                sn = s(regexp(s,'\d'));
                pos = [pos;str2num(sn)];
                labelnames = [labelnames;{strtrim(h5readatt(filename, ...
                    Info.Groups(i).Name, Info.Groups(i).Attributes(j).Name))}];
            else
                %   other attributes
                HSPipeline.(fieldname).(s) = h5readatt(filename, ...
                    Info.Groups(i).Name, Info.Groups(i).Attributes(j).Name);
            end
            
        end
        
        if old_slz_format == 1
            %   make sure the material label are in proper order
            [~, idx] = sort(pos);
            labelnames = labelnames(idx);
            %   put it into labels
            
            HSPipeline.Labels = char(labelnames);
        end
    end

    if isfield(HSPipeline.HDR, 'wavelength')
        HSPipeline.HDR.wavelength = uint16(HSPipeline.HDR.wavelength);
    end

    %   make Micrometers into Nanomemters
    if isfield(HSPipeline.HDR, 'wavelength_unit')    
        if strcmpi(HSPipeline.HDR.wavelength_unit, 'micrometers')
            HSPipeline.HDR.wavelength_unit = 'Nanometers';
            HSPipeline.HDR.wavelength = HSPipeline.HDR.wavelength*1000;
        end
    end
    
    if old_slz_format == 1
        HSPipeline.Endmembers = HSPipeline.Endmembers';
    end
    
    %   check to make sure that there are labels. Some old format slz files
    %   don't have labels
    if ~isfield(HSPipeline, 'Labels') && isfield(HSPipeline.HDR,'numEndmembers')
        s = [];
        n = 0;
        for i = HSPipeline.HDR.numEndmembers:-1:1
            ss = sprintf('mat%d', i);
            if isempty(s)
                n = length(ss);
                s = ss;
            else
                while length(ss) < n
                    ss = sprintf([ss,'%s'], ' ');
                end
                s = [ss;s];
            end
        end
        HSPipeline.Labels = s;
        HSPipeline.Endmembers = HSPipeline.Endmembers';
    end
    %Do the header
    if isfield(HSPipeline,'HDR')
        HDRinfo = h5info(filename,'/HDR/');
        if ~isempty(HDRinfo.Datasets)
            Numatt = length(HDRinfo.Datasets);
            for j = 1:Numatt
                HSPipeline.HDR.(HDRinfo.Datasets(j).Attributes.Name) = HDRinfo.Datasets(j).Attributes.Value;
            end
        end
    end
end

    