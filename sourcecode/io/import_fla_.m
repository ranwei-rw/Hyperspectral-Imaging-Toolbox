% Description:
%
%   [I, H] = import_fla_(filename, xstart, xstop, ystart, ystop)
%   Load a hyperspectral FLA (flat) format file
%
% Syntax:
% 1. Load a hyperspectral I I.fla in ENVI compatible format
%       I = import_fla('I.fla');         % where I is 3D I cube (sample by line by bands)
%       [I, H] = import_fla('I.fla');    % where I is the I cube, H contains header information. This equals to 
%
% 2. Load an hyperspectral I with a downsampling rate
%       I = import_fla('I.fla', 'scale', n) % n is the number for downsampling 
%
% 3. Crop a region from given hyperspectral I file (XSTART, YSTART) and (XSTOP, YSTOP) specify
%    the topleft and the bottomright corner of the region
%       I = import_fla('I.fla', XSTART, XSTOP, YSTART, YSTOP);
%
% 4. Extract one band from given hyperspectral I file. Here band = n contains indicates the nth band.
%       I = import_fla('I.fla', 'band', band);
%
% 5. Header only: Load header information only without reading data
%       H = import_fla('I.fla', 'header'); or
%       H = import_fla('I.hdr');
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.1.1
% Last Update Date: 04 Nov 2014 
% Add support for imec header files. Enhanced performance for some unusual
% cases

%%  Logs
% Version: 1.1.0
% Last Update Date: 24 Oct 2014  
% Enabled reading from files with iregular header information like .cue.

% Version: 1.0.8
% Last Update Date: 23 Oct 2014
% Fixed a bug when there is empty description

% Version: 1.0.7
% Last Update Date: 12 Aug 2014

%   Version 1.0.7
%   Changes: 
%       Disabled autoscale; and
%       Cancelled Full Header option and wavelength only option. Modified header only option.
%

% Version 1.0.6
% Bug fixed:
%   Unable to read commenting lines starts with a ;
%   also it reads Sensor type, file type/filetype properly now

% Version: 1.0.5
% Last Update Date: 29 Oct 2013
%   
%   Bug fixed: open fla file twice in reading. This could lead to failure in deleting the same fla file later.

%  Version 1.0.4
%  This is the file for fucntion import_fla from given ENVI Format Hyperspectral Image files
%  Date: 2012.10.29

%  This file is developed by Ran Wei

% Note: the header file of the hyperspectral I has to be present along with the data file in the
% same folder.

%   HEADER DATA TYPE:
%   1   = 8 bit byte;
%   2   = 16-bit signed integer;
%   3   = 32-bit signed long integer;
%   4   = 32-bit floating point;
%   5   = 64-bit double precision floating point;
%   6   = 2x32-bit complex, real-imaginary pair of double precision;
%   9   = 2x64-bit double precision complex, real-imaginary pair of double precision;
%   12  = 16-bit unsigned integer;
%   13  = 32-bit unsigned long integer;
%   14  = 64-bit signed long integer;
%   15  = 64-bit unsigned long integer.
%%
function [I, H] = import_fla_(filename, xstart, xstop, ystart, ystop)

    %% initialise variable values
    crop              = 0; 
    crop_region       = [];
    scale_image       = 0;
    scaleratio_height = 0;
    scaleratio_width  = 0;
    bslice            = 0;
    headeronly        = 0;
    H                 = {};
    
    %% handling input parameters
    switch nargin
        case 5
            %   when all parameters are numeric, it means to get a rect from the I
            if isnumeric([xstart, xstop, ystart, ystop])
                crop = 1;
                crop_region = [xstart, ystart, xstop-xstart+1, ystop-ystart+1];
            else
                error('Unknown input parameter');
            end
        case 3
            %   case: (filename, rows, cols)
            %   or filename, 'scale', 2
            %   or filename, 'band', 3
            if isnumeric(xstop)
                value = xstop;
            else
                error('Input parameter has wrong type');
            end
            if isnumeric(xstart)
                %   case: (filename, rows, cols)
                scale_image = 1;
                scaleratio_height = xstart;
                scaleratio_width  = value;
                disp('resizing image');
            else
                if strcmpi(xstart, 'scale')
                    scale_image = 1;
                    scaleratio_height = value;
                    scaleratio_width  = value;
                elseif strcmpi(xstart, 'band')
                    bslice = 1;                
                    nband  = value;
                    sn = length(nband);
                    vs = string(nband);
                    s = 'Reading band(s) from given file:';
                    for i = 1:sn
                        s = strcat(s, vs(i));
                        if i ~= sn
                            s = strcat(s, ',');
                        end
                    end
                    disp(s);
                else
                    %   could be wrong option inputed
                    error('unknown input argument');
                end
            end
        case 2
            %   4 possibilities. 3 valid and 1 invalid
            %
            %    1     I = import_fla_(filename, scale);
            %    2     I = import_fla_(filename, rect); where rect = [start_y, start_x, height, width];
            %    3     I = import_fla_(filename, 'header');
            %    4     I = import_fla_(filename, 'allelse');
            if ischar(xstart)
                if strcmpi(xstart, 'header') || strcmpi(xstart, 'headeronly')
                    %   option 3
                    headeronly = 1;
                else
                    %   option 4
                    s = sprintf('Unknown parameter: %s', xstart);
                    disp(s);
                end
            elseif isnumeric(xstart)
                if length(xstart) == 4
                    crop = 1;
                    crop_region = xstart;%[xstart(1), xstart(2), xstart(1)+xstart(3)-1, xstart(2)+xstart(4)-1];
                    s = sprintf('Cropping image - topleft corner [%d, %d] dimension: [%d %d]', xstart(1), ...
                                xstart(2), xstart(3), xstart(4));
                    disp(s);
                elseif length(xstart) == 1
                    scale_image = 1;
                    scaleratio_height = xstart;
                    scaleratio_width  = xstart;
                else
                    s = strcat('Unknown parameters: ', sprintf(' %f ', xstart));
                    error(s);
                end
            else
                error('Unknown parameter');
            end
        case 1
            %   simply read the data with a filename with all normal options            
        otherwise
            error('Error in input arguments. Please check.');
    %   end of switch control
    end
    
    %%    
    [pathstr, name, ext] = fileparts(filename);
    %pathstr has no \ at the end, ext has . as its first letter
    
    if ~strcmpi(ext, '.fla') && ...
       ~strcmpi(ext, '.raw') && ...
       ~strcmpi(ext, '.dat') && ...
       ~strcmpi(ext, '.hdr') && ...
       ~strcmpi(ext, '.dem') && ...
       ~strcmpi(ext, '.img') && ...
       ~strcmpi(ext, '.cue')
        error('Invalid file type.');
    end
    
    % Now we have the data and header file name
    if isempty(pathstr)
        %   incase its an empty
        s = what(pathstr);
        pathstr = s.path;
    end
    
    header_file_name = strcat(pathstr, '\', name, '.hdr');
    %   make sure file separator is correct
    if isequal('\', filesep)
        % most systems use \ as separator
        [~, pos] = strtok(header_file_name, '/');
        if ~isempty(pos)
            header_file_name = strrep(header_file_name, '/', filesep);
        end
    else
        [~, pos] = strtok(header_file_name, '\');
        if ~isempty(pos)
            header_file_name = strrep(header_file_name, '\', filesep);
        end
    end
    
    
    if strcmpi(ext, '.hdr')
        %   when a hdr file is given, we also know the user wants a header only option
        headeronly = 1;
        ext = '.fla';
    end
    
    data_file_name = strcat(pathstr, '\', name, ext);
    if isequal('\', filesep)
        % most systems use \ as separator
        [~, pos] = strtok(data_file_name, '/');
        if ~isempty(pos)
            data_file_name = strrep(data_file_name, '/', filesep);
        end
    else
        [~, pos] = strtok(data_file_name, '\');
        if ~isempty(pos)
            data_file_name = strrep(data_file_name, '\', filesep);
        end
    end
    
    if exist(data_file_name, 'file') ~= 2
        s = strcat(data_file_name, {' '}, 'does not exist');
        error(s{1});
    end
    %   check if the data file exist:
    if ~headeronly
        data_file_fp = fopen(data_file_name);
        if data_file_fp == -1
            %   when it fails to open the given file, prompt error message
            error('Failed to open file %s', data_file_name);
        end
        fclose(data_file_fp);
    end    

    % check if header file exists
    data_file_fp = fopen(header_file_name);
    if data_file_fp == -1
        %   when it fails to open the given file, prompt error message
        error('Failed to open file %s', header_file_name);
    end
    % %     Interleave Mode options:
    % %     interleave == 0 <- bsq
    % %     interleave == 1 <- bil
    % %     interleave == 2 <- bip
    
    %define a variable to monitor the input process
    % if continuing_reading is of a value of larger than 0, 
    % any input will be put in corresponding variables;
    % continuing_reading = 1 ==> is reading things for description
    EMPTY_STR               = 'EMPTY_STR';
    acquisition_time        = EMPTY_STR;
    wavelength              = [];
    fwhm                    = [];
    default_bands           = [];
    allothers               = []; 
    byteorder               = 0;
    interleave              = 0;
    reflectance_scaler      = -1;
    description             = {};
    map_info                = EMPTY_STR;
    projection_info         = EMPTY_STR;
    sensor                  = EMPTY_STR;
    header_offset           = 0;
    data_type_word          = 0;
    lines                   = 0;
    samples                 = 0;
    bands                   = 0;
    z_plot_titles           = EMPTY_STR;
    wave_unit               = EMPTY_STR;
    z_plot_range            = [];
    z_plot_average          = [];
    pixel_size              = [];
    default_stretch         = {};
    band_names              = {};
    item                    = '';
    value                   = '';
    complete                = 0;
    rest                    = '';
    buffer                  = '';
    filetype                = EMPTY_STR;
    weights                 = 0;
    encoding                = EMPTY_STR;
    abundances              = 0;
    classes                 = 0;
    classlookup             = [];
    classnames              = {};
    x_start                 = -1;
    y_start                 = -1;
    %%support for imec entries
    maxvalue                = -1;
    fps                     = -1;
    tint                    = -1;
    errors                  = '';
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Reading header file starts
    %
    %   The fist step of this procedure is to make sure the input buffer 
    %   is always in the format:
    %   item name = content or
    %   item name = { content } white space is optional
    %   if the above format is not satisfied, read a new line can put together
    %   the content and re-check
    %
    %   read in the first line
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    line = 0;    
    while (line ~= -1)
        
        if isempty(rest) ~= 1
            line = rest;
            rest = '';
        else
            line = fgets(data_file_fp);
            if (line == -1)
                break;
            end
        end
        buffer = strtrim(line);
        if isempty(buffer)
            continue;
        end
        
        %   check whether this string is started with a ;
        iscomm = strfind(buffer, ';');
        if ~isempty(iscomm)
            if iscomm(1) == 1
                if isempty(allothers)
                    allothers = buffer;
                else
                    allothers = strcat(allothers, ';;', buffer);
                end
                continue;
            end
        end
        [complete, item, value, rest] = compose_pattern_(buffer);
        
        while complete ~= 1
            line = fgets(data_file_fp);
            if line == -1
                error(strcat('Faild to parse: ', buffer));
                break;
            else
                %   here we need to consider special cases. For example first
                %   line is sensor = which is simply an empty value but
                %   the second line is sample = 1000; is already another line
                %   of entry. We need to prevent joining these two line
                %   check whether the new line is a complete sentence. If it is, 
                %   then this previous line has to be set to complete
                [complete, ~, ~, ~] = compose_pattern_(line);
                
                if complete
                    %   put the new line into rest for next circle of
                    %   process
                    
                    %   check the last letter is '=' or not
                    if isequal(buffer(end), '=')
                        buffer = strcat(buffer, '{}');
                    else
                        buffer = strcat(buffer, '={}');
                    end
                    [complete, item, value, ~] = compose_pattern_(buffer);
                    rest = line;
                else
                    buffer = strcat(buffer, line);
                    buffer = strtrim(buffer);
                    [complete, item, value, rest] = compose_pattern_(buffer);
                end
            end
        end
        
        %   now we should have item as entry name and value as the input string for the item
        
        % put it into buffer to see whether it fits the pattern of
        % item name = content or item name = {content}
        standard = '';
        if strcmpi(item, 'envi')
            %   this means the header is designed using ENVI format standard
            standard = 'ENVI';
            continue;
        elseif strcmpi(item, 'description')
            %   create a cell with fields
            description = {'description' value};
            continue;
        elseif strcmpi(item, 'acquisition time') || strcmpi(item, 'acquisition date')
            acquisition_time = sscanf(value, '%s');
            continue;
        elseif strcmpi(item, 'header offset')
            header_offset = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'samples')
            samples = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'lines')
            lines = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'bands')
            bands = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'x start')
            x_start = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'y start')
            y_start = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'data type')
            data_type_word = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'sensor type')
            sensor = value;
            continue;
        elseif strcmpi(item, 'band names')
            band_names = sscanf(value, '%s');
            continue;
        elseif strcmpi(item, 'reflectance scale factor')
            reflectance_scaler = sscanf(value, '%f');
            continue;
        elseif strcmpi(item, 'byte order')
            byteorder = sscanf(value, '%d');
            continue;
        elseif strcmpi(item, 'wavelength units')
            wave_unit = value;
            continue;
        elseif strcmpi(item, 'interleave')
            if strcmpi(value, 'bsq')
                interleave = 0;
                continue;
            elseif strcmpi(value, 'bil')
                interleave = 1;
                continue;
            elseif strcmpi(value, 'bip')
                interleave = 2;
                continue;
            else
                error('unknown data interleave mode: %s', buffer);
            end
        elseif strcmpi(item, 'fwhm')
            %   get the rest of the line
            [range, count] = sscanf(value, '%f, ', inf);
            if count >= 1
                fwhm = range;
            end
            continue;
        elseif strcmpi(item, 'wavelength')
            %   get the rest of the line
            %   remove all spaces from it
            i = isspace(value);
            value = value(~i);
            [range, count] = sscanf(value, '%f,', inf);
            if count >= 1
                wavelength = range;
            end
            continue;
        elseif strcmpi(item, 'default bands')
            [range, count] = sscanf(value, '%f, ', inf);
            if count >= 1
                default_bands = range;
            end
            continue;
        elseif strcmpi(item, 'z plot titles')
            z_plot_titles = value;
            continue;
        elseif strcmpi(item, 'file type') || strcmpi(item, 'filetype')
            filetype = value;
            continue;            
        elseif strcmpi(item, 'z plot range')
            [range, count] = sscanf(value, '%f, ', inf);
            if count >= 1
                z_plot_range = range;
            end
            continue;
        elseif strcmpi(item, 'z plot average')
            [range, count] = sscanf(value, '%f, ', inf);
            if count >= 1
                z_plot_average = range;
            end
            continue;
        elseif strcmpi(item, 'default stretch')
            [range, count] = sscanf(value, '%s%c', inf);
            if count >= 1
                default_stretch = range;
            end
            continue;
        elseif strcmpi(item, 'pixel size')
            [range, count] = sscanf(value, '%f, ', inf);
            if count >= 1
                pixel_size = range;
            end
            continue;
        elseif strcmpi(item, 'projection info')
            projection_info = value;
            continue;
        elseif strcmpi(item, 'map info')
            map_info = value;
            continue;
        elseif strcmpi(item, 'encoding')
            encoding = value;
            continue;
        elseif strcmpi(item, 'weights')
            [range, count] = sscanf(value, '%f', inf);
            if count >= 1
                weights = range;
            end
            continue;
        elseif strcmpi(item, 'abundances')
            [range, count] = sscanf(value, '%f', inf);
            if count >= 1
                abundances = range;
            end
            continue;
        elseif strcmpi(item, 'classes')
            [range, count] = sscanf(value, '%f', inf);
            if count >= 1
                classes = range;
            end
            continue;    
        elseif strcmpi(item, 'class lookup')
            [range, count] = sscanf(value, '%f, ', inf);
            if count >= 1
                classlookup = range;
            end
            continue;    
        elseif strcmpi(item, 'class names')
            [range, count] = recog_string(value, ', ');
            if count >= 1
                classnames = range;
            end
            continue; 
        elseif strcmpi(item, 'max value')
            [range, count] = sscanf(value, '%f', inf);
            if count >= 1
                maxvalue = range;
            end
            continue;
        elseif strcmpi(item, 'fps')
            fps = sscanf(value, '%f', inf);
            continue;
        elseif strcmpi(item, 'tint')
            tint = sscanf(value, '%f', inf);
            continue;
        elseif strcmpi(item, 'errors')
            errors = sscanf(value, '%s');
            continue;
        else
            if isempty(allothers)
                allothers = buffer;
            else
                allothers = strcat(allothers, ';;', buffer);
            end
        end
        
    end
    fclose(data_file_fp); 
    %%
    %Get the interleave code in
    if interleave == 0 
       intleave = 'BSQ';
    elseif interleave == 1
       intleave = 'BIL';
    elseif interleave == 2 
       intleave = 'BIP';
    end

    if isempty(description)
        description = {'description' 'BLANKDATA'};
    end
    
    %Go on with the header
    if ~isempty(classnames) && ...
       ~isempty(classes) && ...
       ~isempty(weights) && ...
       strcmpi(filetype, 'NICTA HSPipeline')
        H = struct('samples', samples, 'lines', lines, 'bands', bands, 'description', description{2}, ...
                 'header_offset', header_offset, 'byteorder', byteorder, 'interleave', intleave, ...
                 'encoding', encoding, 'weights', weights, 'abundances', abundances, ...
                 'classes', classes, 'classlookup', classlookup, 'classnames', classnames);
    elseif strcmpi(filetype, 'NICTA HSPipeline') && ~isequal(encoding, EMPTY_STR)
        H = struct('samples', samples, 'lines', lines, 'bands', bands, 'description', description{2}, ...
                 'header_offset', header_offset, 'byteorder', byteorder, 'interleave', intleave, ...
                 'encoding', encoding, 'allothers', allothers, 'filetype', filetype);
    else        
        H = struct('samples', samples, 'lines', lines, 'bands', bands, 'wavelength', wavelength, ...
                  'datatype', data_type_word, 'interleave', intleave, 'description', description{2}, ...
                  'header_offset', header_offset, 'byteorder', byteorder);
              
        if ~isequal(sensor, EMPTY_STR)
            H.('sensor') = sensor;
        end
        
        if ~isequal(wave_unit, EMPTY_STR)
            H.('wave_unit') = wave_unit;
        end
        if ~isequal(filetype, EMPTY_STR)
            H.('filetype') = filetype;
        end
        if ~isempty(fwhm)
            H.('fwhm') = fwhm;
        end 
        if ~isequal(acquisition_time, EMPTY_STR)
            H.('acquisition_time') = acquisition_time;
        end 
        if ~isempty(default_bands)
            H.('default_bands') = default_bands;
        end
        
        if ~isempty(reflectance_scaler) && reflectance_scaler ~= -1
            H.('reflectance_scaler') = reflectance_scaler;
        end

        if ~isequal(map_info, EMPTY_STR)
            H.('map_info') = map_info;
        end
        if ~isequal(projection_info, EMPTY_STR)
            H.('projection_info') = projection_info;
        end
        if ~isequal(z_plot_titles, EMPTY_STR)
            H.('z_plot_titles') = z_plot_titles;
        end
        if ~isempty(z_plot_range)
            H.('z_plot_range') = z_plot_range;
        end
        if ~isempty(z_plot_average)
            H.('z_plot_average') = z_plot_average;
        end
        if ~isempty(pixel_size)
            H.('pixel_size') = pixel_size;
        end
        if ~isempty(default_stretch)
            H.('default_stretch') = default_stretch;
        end
        if ~isempty(band_names)
            H.('band_names') = band_names;
        end
        if ~isempty(allothers)
            H.('allothers') = allothers;
        end
        %   updates for imec entries
        if maxvalue ~= -1
            H.('maxvalue') = maxvalue;
        end
        
        if x_start ~= -1
            H.('x_start') = x_start;
        end
        
        if y_start ~= -1
            H.('y_start') = x_start;
        end
        if fps ~= -1
            H.('fps') = fps;
        end
        if tint ~= -1
            H.('tint') = tint;
        end
        
        if ~isempty(errors)
            H.('errors') = errors;
        end
    end
    
    %   make wavelength to a row vector
    H.wavelength = reshape(H.wavelength, [1, length(H.wavelength)]);
    
    if headeronly == 1
        I = H;
        disp('Header information is read. Actual image data is skipped');
        return;
    end
 %%       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Reading data file starts
    %
    %   Prepare parameters
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if data_type_word == 4      % float ENVI datatype 4
        data_type = 'float32';
    elseif data_type_word == 12 % short, ENVI datatype 12
        data_type = 'uint16';
    elseif data_type_word == 2 % short, ENVI datatype 2
        data_type = 'int16';
    elseif data_type_word == 1  % byte, ENVI datatype 1
        data_type = 'char';
    elseif data_type_word == 5  % byte, ENVI datatype 5
        data_type = 'double';
    end
   
    height = lines;
    width = samples;
    
    if ~byteorder
        data_file_fp = fopen(data_file_name);
    else
        data_file_fp = fopen(data_file_name, 'r', 'b');
    end

    if data_file_fp == -1
        %   when it fails to open the given file, prompt error message
        error('Failed to open file %s', data_file_name);
    end

    %% bsq mode
    if interleave == 0 % bsq mode
        if crop == 0 && scale_image == 0 && bslice == 0 % read whole I
            I = zeros(height, width, bands);
            for k = 1:bands
                It = fread(data_file_fp, width*height, data_type);
                I(:, :, k) = reshape(It, width, height)';           
            end
        elseif scale_image == 1
            %   here are some differences. if scaleratio_height and scaleratio_width is less than 1
            %   we know we need to shrink image and actual size should be calculated
            %   if scaleratio_height and scaleratio_width is greater than 1, then they are the 
            %   dimensions we want
            I = zeros(height, width, bands);
            for k = 1:bands
                It = fread(data_file_fp, width*height, data_type);
                I(:, :, k) = reshape(It, width, height)';           
            end
            %   then resize it
            if scaleratio_height < 1 && scaleratio_width < 1
                %   this is to get shrink to a percentage
                scaleratio_height = floor(height*scaleratio_height);
                scaleratio_width  = floor(width*scaleratio_width);
            end
            s = sprintf('Resized image: %d by %d', scaleratio_height, scaleratio_width);
            disp(s);
            I = resize_(I, scaleratio_height, scaleratio_width);
            H.lines = scaleratio_height;
            H.samples = scaleratio_width;
        elseif bslice == 1 % read a slice
            sn = length(nband);
            I = zeros(height, width, sn);
            for ni = 1:sn
                current_band = nband(ni);
                
                if current_band > bands
                    fclose(data_file_fp);
                    error('prog:input', 'Wrong band number %d is given', nband);
                end

                if strcmpi(data_type, 'uint16')
                    fseek(data_file_fp, width*height*2*(current_band-1), 'bof');
                elseif strcmpi(data_type, 'char')
                    fseek(data_file_fp, width*height*(current_band-1), 'bof');
                else
                    fseek(data_file_fp, width*height*4*(current_band-1), 'bof');    
                end
                F = fread(data_file_fp, width*height, data_type);
                I(:, :, ni) = reshape(F, width, height)';
                H.wavelength(ni) = H.wavelength(current_band);
            end
            H.bands = reshape(nband, [1, sn]);
            H.wavelength = H.wavelength(1:sn);
        else        % read a specified I region
            if crop_region(1) < 1 || crop_region(2) < 1 
                fclose(data_file_fp);
                error('Given parameters are out of I range');
            end
            if crop_region(1)+crop_region(3)-1 > H.samples
                crop_region(3) = H.samples + 1 - crop_region(1);
                H.samples = crop_region(3); 
            end
            
            if crop_region(2)+crop_region(4)-1>H.lines
                crop_region(4) = H.lines - crop_region(2) + 1;
                H.lines = crop_region(4);
            end
            
            x = crop_region(1);
            y = crop_region(2);
            w = crop_region(3);
            h = crop_region(4);
            if data_type_word == 4      % float ENVI datatype 4
                bitdepth = 4;
            elseif data_type_word == 12 % short, ENVI datatype 12
                bitdepth = 2;
            elseif data_type_word == 2 % short, ENVI datatype 2
                bitdepth = 2;
            elseif data_type_word == 1  % byte, ENVI datatype 1
                bitdepth = 1;
            elseif data_type_word == 5  % byte, ENVI datatype 5
                bitdepth = 4;
            end
%             %   Was like this. But larger files could have very large
%             %   dimensions. So creating image array for whole image is not
%             %   an option 
%             I = zeros(height, width, bands);
            I = zeros(h, w, bands);
            for k = 1:bands
                %   jump to start position for each band
                %   firstly, go to the end of previous line
                fseek(data_file_fp, width*(y-1)*bitdepth, 'cof');
                for m = y:y+h-1
                    %   just to begin position for each row
                    fseek(data_file_fp, (x-1)*bitdepth, 'cof');
                    %   read data
                    I(m-y+1, :, k) = fread(data_file_fp, w, data_type)';
                    %   jump to the end of current line
                    fseek(data_file_fp, (H.samples-(x+w))*bitdepth, 'cof');
                end
                %   jump to the end of current band
                fseek(data_file_fp, (H.lines-(y+h))*bitdepth, 'cof');
            end
%                 It = fread(data_file_fp, width*height, data_type);
%                 I(:, :, k) = reshape(It, width, height)';           
%             end
%             I = crop_(I, crop_region);
%             H.lines = crop_region(4);
%             H.samples = crop_region(3); 
        end
    elseif interleave == 1 
        %% %bil mode
        %   no matter what choise will be use, we read in the whole I first
        I = zeros(height, width, bands);
        for j = 1:height
            for k = 1:bands
                I(j, :, k) = fread(data_file_fp, width, data_type);
            end
        end
            
        if crop == 0 && scale_image == 0 && bslice == 0, % read whole I
            fclose(data_file_fp);
            %disp('Image data is read');
            return;
        elseif scale_image == 1 % read I and scale
            %   we create a new matrix and put the elements into 
            %   remember, matrix element is called by (Y, X, Z) order
            if scaleratio_height < 1 && scaleratio_width < 1
                %   this is to get shrink to a percentage
                scaleratio_height = floor(height*scaleratio_height);
                scaleratio_width  = floor(width*scaleratio_width);
            end
            s = sprintf('Resized image: %d by %d', scaleratio_height, scaleratio_width);
            disp(s);
            I = resize_(I, scaleratio_height, scaleratio_width);
            H.lines = scaleratio_height;
            H.samples = scaleratio_width;            
        elseif bslice == 1 % read a slice
            if nband > bands
                fclose(data_file_fp);
                error('Wrong band number is given');
            end
            I = I(:, :, nband);
            H.wavelength = H.wavelength(nband);
            H.bands = 1;
        else % read a specified I region
            if crop_region(1) < 1 || crop_region(2) < 1 || ...
               crop_region(1)+crop_region(3)-1 > H.lines || crop_region(2)+crop_region(4)-1>H.samples
                fclose(data_file_fp);
                error('Given parameters are out of I range');
            end

            I = zeros(height, width, bands);
            for k = 1:bands
                It = fread(data_file_fp, width*height, data_type);
                I(:, :, k) = reshape(It, width, height)';           
            end
            I = crop_(I, crop_region);
            H.lines = crop_region(4);
            H.samples = crop_region(3);
        end
        fclose(data_file_fp);
        % disp('Image data is read');
        return;
    else
        %% % bip mode
        %   no matter what choise will be use, we read in the whole I first
        I = zeros(height, width, bands);
        
        for j = 1:height
            for i = 1:width
                for k = 1:bands
                    I(j, i, k) = fread(data_file_fp, 1, data_type);
                end
            end
        end
        
        if crop == 0 && scale_image == 0 && bslice == 0 % read whole I
            fclose(data_file_fp);
            % disp('Image data is read');
            return;
        elseif scale_image == 1
            if scaleratio_height < 1 && scaleratio_width < 1
                %   this is to get shrink to a percentage
                scaleratio_height = floor(height*scaleratio_height);
                scaleratio_width  = floor(width*scaleratio_width);
            end
            s = sprintf('Resized image: %d by %d', scaleratio_height, scaleratio_width);
            disp(s);
            I = resize_(I, scaleratio_height, scaleratio_width);
            H.lines = scaleratio_height;
            H.samples = scaleratio_width;
        elseif bslice == 1 % read a slice
            if nband > bands
                fclose(data_file_fp);
                error('Wrong band number is given');
            end

            I = I(:, :, nband);
            H.wavelength = H.wavelength(nband);
            H.bands = 1;
        else % read a specified I region
            if crop_region(1) < 1 || crop_region(2) < 1 || ...
               crop_region(1)+crop_region(3)-1 > H.lines || crop_region(2)+crop_region(4)-1>H.samples
                fclose(data_file_fp);
                error('Given parameters are out of I range');
            end

            I = zeros(height, width, bands);
            for k = 1:bands
                It = fread(data_file_fp, width*height, data_type);
                I(:, :, k) = reshape(It, width, height)';           
            end
            I = crop_(I, crop_region);
            H.lines = crop_region(4);
            H.samples = crop_region(3);
        end
        fclose(data_file_fp);
        % disp('Image data is read');
        return;
    end

    fclose(data_file_fp);
    %disp('Image data is read');
    return
end

%%  function crop_

function Q = crop_(I, rect)

    if ndims(I)==3
        [~, ~, bands] = size(I);    
        for i = 1:bands
             Q(:, :, i) = double(I(round(rect(2)):round(rect(2)+rect(4)-1), ...
                                   round(rect(1)):round(rect(1)+rect(3)-1), i));
        end
    else
        Q(:, :) = double(I(round(rect(2)):round(rect(2)+rect(4)-1), ...
                           round(rect(1)):round(rect(1)+rect(3)-1)));
    end

end

%%  function resize_
function Q = resize_(I, rows, cols)

    [~, ~, bands] = size(I);
    
    for i = 1:bands
        if exist('rows', 'var')
            if ~exist('cols', 'var')
                Q(:, :, i) = imresize(I(:, :, i), rows, 'nearest');
            else
                Q(:, :, i) = imresize(I(:, :, i), [rows cols], 'nearest');
            end
        else
            error('Please provide a scaling value or a valid image size');
        end
    end

end


%%  function compose_pattern_
function [satisfied, item, value, leftover] = compose_pattern_(content)
    %   satisfied = 1 when the input content satisfies the requirement 
    %               of item = content or item = {content}
    %             = 0 when it dosen't
    %   if there is any leftover liek item = {content} this_is_leftover
    %   then it will be put in the variable leftover. Otherwise, it will be empty.
    %   only makes sense when return satisfied == 1
    %   item contains the entry name if this is a good entry
    %   value will be the value for this item. Without {} no matter whether is has in original
    %   content
    
    %   left is the number { found in content
    %   right is the number } function in content
    
    satisfied   = 0;
    leftover    = '';
    left        = 0;
    right       = 0;
    item        = '';
    value       = '';
    brace_left  = '{';
    brace_right = '}';
    
    if isempty(content)
        satisfied = 0;        
        return;
    end
    
    content = strtrim(content);
    
    if isempty(content)
        satisfied = 0;        
        return;
    end
    
    %   search for the mark ='
    app = regexp(content, '=', 'start');
    if isempty(app) == 1
        %   there is no = existing in the content.
        %   return false
        %   the only exception is ENVI
        satisfied = strcmpi(content, 'ENVI');
        item = content;
        return;
    end
    
    %now we have to parts: left part of buffer and right part of buffer
    leng = length(content);
    
    if app(1) == leng
        %   when the right portion is empty (right side of =), it certains fail
        return;
    end
                    
    %   at this stage, we are sure that we have both of two sides of =
    %   check for any existance of {
    item = content(1 : app(1) - 1);
    value = content(app(1)+1 : leng);
    item = strtrim(item);
    %   now we have two parts saved in item and value
    app = regexp(value, brace_left, 'start');
    
    if isempty(app)
        %   no { is found at all, which means it fits in the pattern of item = value
        item = strtrim(item);
        value = strtrim(value);
        satisfied = 1;
        return;
    end
    
    %   now we reach the point that value contains at least one {
    s = size(app);
    left = s(end);

    appl = regexp(value, brace_right, 'start');
    if isempty(appl)
        %   input value string has { but no } => incomplete
        return;
    end
    
    s = size(appl);
    right = s(end);
    if left ~= right
        %   incomplete 
        return;
    else
        %   matched number of { and }
        %   remove them, cut the part at right side of } to leftover
        leng = length(value);
        if appl(end) ~= leng
            %when there is more outside the rightmost }
            leftover = value(appl(end)+1:leng);
            leftover = strtrim(leftover);
        end
        value = value(app(1)+1 : appl(end)-1);
        value = strtrim(value);
        satisfied = 1;
        return;
    end
    
    
%   end of function import_fla_ 
end