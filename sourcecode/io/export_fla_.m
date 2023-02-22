
%   This function is designed to output hyperspectral image data into a standard ENVI file (fla/hdr pair)
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2015 All Rights Reserved.
% Author: Ran Wei
% Version: 1.1.3
% Last Update Date: 17 Apirl 20145
% Added section to handle when users input only a data cube without header
% information

% Version: 1.1.2
% Last Update Date: 17 Nov 2014 
% small bug fixed in filename part

% Version: 1.1.1
% Last Update Date: 04 Nov 2014 
% Add support for imec header files. Enhanced performance for some unusual
% cases

% Version: 1.1.0
% Last Update Date: 24 Oct 2014  
% Enabled output files with different given filetype, like .cue, .raw, etc

% Version: 1.0.9
% Last Update Date: 23 Oct 2014  
% Fixed output problem for bip and bil formats
%
% Version: 1.0.8
% Last Update Date: 23 Oct 2014  
% Fixed empty elements like description, sensor and reflectance scale factor

% Version: 1.0.7
% Last Update Date: 18 Sep 2014

%   Changed done in 1.0.7
%       argument order is swapped, from DATA, filename to filename, DATA
%   Bug fixed in 1.0.6
%       output wavelength properly for float type
%       output sensor type, filetype and allother properly
%   % Last Update Date: 29 May 2014
%
% Version: 1.0.5
% Last Update Date: 29 Oct 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Handle Filename structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_fla_(filename, DATA, Overwrite)
    if isa(filename, 'string')
        filename = char(filename);
    end    
    overwrite_option_is_given = false;

    if exist('Overwrite', 'var')
        overwrite_option_is_given = true;
    end

    %   1. keyword: type
    %   check whether it exists a . in given file name
    app = regexp(filename, '\.', 'start');
    leng = length(filename);
    suffix = '';
    if ~isempty(app)
        if app(1) == leng
            %   if the file name has a . but not any suffix like: `abc.'
            filename = filename(1 : app(1)-1);
        elseif app(1) == 1
            %   situations like ./ .\ .. etc 
            if length(app) == 1
                %   this is the only location. which means no poxtfix for this item.
                if filename(2) == '/' || filename(2) == '\'
                    suffix = '';
                else
                    error('Could not process this filename');
                end
            else
                %   go for the last . detected
                if app(end) == leng
                    filename = filename(1 : app(1)-1);
                else
                    %   use the last . as the separator of postfix
                    suffix = filename(app(end) : leng);
                    filename = filename(1 : app(end)-1);
                end
            end
        else

            suffix = filename(app(1) : leng);
            filename = filename(1 : app(1)-1);

            %   Tell from the data type to see whether a fla or a hsz is required
            if isfield(DATA, 'main')
                %   there is a cell called main in DATA -==> this is a hsz file
                if strcmpi(suffix, '.fla') || strcmpi(suffix, '.hdr')
                    error('Could not save hsz structure as a .fla file. Please use .hsz suffix');
                end
            end
        end
    end


    %   when is a header file extension is given
    if isempty(suffix)
        header_name = strcat(filename, '.hdr');
        data_name   = strcat(filename, '.fla');
    else
        if isequal(suffix, '.hdr')
            header_name = strcat(filename, '.hdr');
            data_name   = strcat(filename, '.fla');
        else
            header_name = strcat(filename, '.hdr');
            data_name   = strcat(filename, suffix);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Prepare for writing out: Check existance of files
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist(header_name, 'file') || exist(data_name, 'file')
        if overwrite_option_is_given == true 
            if Overwrite == false
                disp('File exists. Exiting');
                return;
            else
                disp('File exists. Overwriting');
            end
        else
            text = sprintf('File already exists. Replace it? [Y/N] (Default to No)');
            % make sure \ is replaced by \\
            text = strrep(text, '\', '\\');
            choice = input(text, 's');
            if isempty(choice)
                %   default value is N, overwrite it
                choice = 'N';
            end
            if strcmpi(choice, 'Y') ~= 1
                %   when user doesn't want to overwrite
                text = sprintf('User does not want to replace it. Exiting\n');
                % make sure \ is replaced by \\
                text = strrep(text, '\', '\\');
                fprintf(text);
                return;
            else
                disp('Overwriting');
            end
        end
    end

    if isstruct(DATA)
        if ~isfield(DATA, 'HDR') || ~isfield(DATA, 'I')
            error('The data that is provided is not a valid fla/hdr structure');
        end
    else
        %   it means that the input is only a matrix
        if isnumeric(DATA)
            %   this is a data matrix, without valid HDR information
            disp('No header information found. Generating from data given');
            I.I = DATA;
            [height, width, bands] = size(DATA);
            s = sprintf('Lines: %d', height);
            disp(s);
            I.HDR.lines   = height;
            s = sprintf('samples: %d', width);
            disp(s);
            I.HDR.samples = width;
            s = sprintf('bands: %d', bands);
            disp(s);
            I.HDR.bands   = bands;
            disp('Default interleave: BSQ');
            I.HDR.interleave = 'bsq';
            I.HDR.description = 'Information deduced from data given';
            disp('Default data type: double');
            I.HDR.datatype = 5;
            disp('Wavelength info: None');
            DATA = I;
        else
            error('The data that is provided is not a valid fla/hdr structure');
        end
        
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Before saving to a file, check the data structure for certain compulsory items in the header 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(DATA.HDR, 'samples')
        error('hstoolbox:io:flaread:noSamplesFound', 'No information about Samples is found');
    end

    if ~isfield(DATA.HDR, 'lines')
        error('hstoolbox:io:flaread:noLinesFound', 'No information about lines is found');
    end

    if ~isfield(DATA.HDR, 'bands')
        error('hstoolbox:io:flaread:noBandsFound', 'No information about bands is found');
    end

    if ~isfield(DATA.HDR, 'interleave')
        error('hstoolbox:io:flaread:noInterleaveFound', 'No information about interleave is found');
    end

    if ~isfield(DATA.HDR, 'datatype')
        error('hstoolbox:io:flaread:noDatatypeFound', 'No information about datatype is found');
    end

    %   From this point, the data structure inside of DATA should be valid for fla/hdr

    %   check whether the data has a valid structure.

    

    %   output file names are verified to be correct

    %   two options: output as a HSZ or a pair of fla/hdr

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % open and write header file
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fp = fopen(header_name, 'w');
    if fp == -1
        error('Failed to open file %s', header_name);
    end
    %   HEADER DATA TYPE:
    %   1   = 8 bit byte;                                                                   uint8
    %   2   = 16-bit signed integer;                                                        int16
    %   3   = 32-bit signed long integer;                                                   int32
    %   4   = 32-bit floating point;                                                        single
    %   5   = 64-bit double precision floating point;                                       double
    %   6   = 2x32-bit complex, real-imaginary pair of double precision;
    %   9   = 2x64-bit double precision complex, real-imaginary pair of double precision;
    %   12  = 16-bit unsigned integer;                                                      uint16
    %   13  = 32-bit unsigned long integer;                                                 uint32
    %   14  = 64-bit signed long integer;                                                   int64
    %   15  = 64-bit unsigned long integer.                                                 uing64

    %   define main variables:
    height  = 0;
    width   = 0;
    bandnum = 0;
    % height = lines, width = samples
    fprintf(fp, 'ENVI\r\n');
    if isfield(DATA.HDR, 'description')
        if ~strcmpi(DATA.HDR.description, 'BLANKDATA')
            text = 'Description = {';
            text = [text DATA.HDR.description '}'];
            fprintf(fp, '%s\r\n', text);
        end
    end
    
    if isfield(DATA.HDR, 'sensor')
        fprintf(fp, 'sensor type = %s\r\n', DATA.HDR.sensor);
    end

    if isfield(DATA.HDR, 'filetype')
        fprintf(fp, 'file type = %s\r\n', DATA.HDR.filetype);
    else
        fprintf(fp, 'file type = ENVI Standard\r\n');
    end

    fprintf(fp, 'samples = %d\r\n', DATA.HDR.samples);
    fprintf(fp, 'lines = %d\r\n', DATA.HDR.lines);
    fprintf(fp, 'bands = %d\r\n', DATA.HDR.bands);
    fprintf(fp, 'interleave = %s\r\n', DATA.HDR.interleave);
    fprintf(fp, 'data type = %d\r\n', DATA.HDR.datatype);

    if isfield(DATA.HDR, 'wavelength')
        fprintf(fp, 'wavelength = {\r\n');
        num = max(size(DATA.HDR.wavelength));
        wa = reshape(DATA.HDR.wavelength, num, 1);
        if isfloat(DATA.HDR.wavelength)
            for i = 1:num-1
                fprintf(fp, '              %f,\r\n', wa(i));
            end
            fprintf(fp, '              %f\r\n', wa(num));
        else
            for i = 1:num-1
                fprintf(fp, '              %d,\r\n', wa(i));
            end
            fprintf(fp, '              %d\r\n', wa(num));
        end
        fprintf(fp,'}\r\n');
    end

    if isfield(DATA.HDR, 'header_offset')
        fprintf(fp, 'header offset = %d\r\n', DATA.HDR.header_offset);
    end

    if isfield(DATA.HDR, 'byteorder')
        fprintf(fp, 'byte order = %d\r\n', DATA.HDR.byteorder);
    end

    if isfield(DATA.HDR, 'reflectance_scaler')
        fprintf(fp, 'reflectance scale factor = %f\r\n', DATA.HDR.reflectance_scaler);
    end

    if isfield(DATA.HDR, 'wave_unit')
        fprintf(fp, 'wavelength units = %s\r\n', DATA.HDR.wave_unit);
    end
    if isfield(DATA.HDR, 'acquisition_time')
        fprintf(fp, 'acquisition time = %s\r\n', DATA.HDR.acquisition_time);
    end
    if isfield(DATA.HDR, 'z_plot_titles')
        fprintf(fp, 'z plot titles = %s\r\n', DATA.HDR.z_plot_titles);
    end
    
    if isfield(DATA.HDR, 'default_bands')
        m = max(size(DATA.HDR.default_bands));
        if m ~= 0
            fprintf(fp, 'default bands = {');
            for i = 1:m
                if i == m
                    fprintf(fp, '%d}\r\n', DATA.HDR.default_bands(i));
                else
                    fprintf(fp, '%d, ', DATA.HDR.default_bands(i));
                end
            end
        end
    end

    %   entries for imecs
    if isfield(DATA.HDR, 'maxvalue')
        fprintf(fp, 'max value = %f\r\n', DATA.HDR.maxvalue);
    end
    if isfield(DATA.HDR, 'x_start')
        fprintf(fp, 'x start = %d\r\n', DATA.HDR.x_start);
    end
    if isfield(DATA.HDR, 'y_start')
        fprintf(fp, 'y start = %d\r\n', DATA.HDR.y_start);
    end
    if isfield(DATA.HDR, 'fps')
        fprintf(fp, 'fps = %f\r\n', DATA.HDR.fps);
    end
    if isfield(DATA.HDR, 'tint')
        fprintf(fp, 'tint = %f\r\n', DATA.HDR.tint);
    end
    if isfield(DATA.HDR, 'errors')
        fprintf(fp, 'errors = {%s}\r\n', DATA.HDR.errors);
    end
    
    if isfield(DATA.HDR, 'allothers')
        fprintf(fp, '%s\r\n', DATA.HDR.allothers);
    end

    fclose(fp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     Finished writing out header file
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % write data file
    fp = fopen(data_name, 'w');
    if fp == -1
        error('Failed to open file %s', data_name);
    end
    word_name = '';

    switch DATA.HDR.datatype
        case 12
            word_name = 'uint16';
        case 2
            word_name = 'int16';
        case 1
            word_name = 'char';
        case 4
            word_name = 'float';    
        case 5
            word_name = 'double';    
        otherwise
            word_name = 'char';
    end

    %   it is assumed that actual image data is always arranged in BSQ mode
    %   if the output interleave format is different, care should be taken.
    %   In summary, all data (if they are not in right format, convert it
    %   to the right mode and output them in one action
    switch upper(DATA.HDR.interleave)
        case 'BIP'
            % Band-interleaved-by-pixel
            % need to convert from BSQ to BIP and then output in one batch
            d = zeros(DATA.HDR.bands * DATA.HDR.samples * DATA.HDR.lines, 1);
            counter = 1;
            for j = 1:DATA.HDR.lines
                for i = 1:DATA.HDR.samples
                    d(((counter-1)*DATA.HDR.bands+1):counter*DATA.HDR.bands) =  reshape(DATA.I(j, i, :),  [DATA.HDR.bands, 1]);
                    counter = counter+1;
                end
            end
            fwrite(fp, d, word_name);
        case 'BIL'
            % Band-interleaved-by-line
            d = zeros(DATA.HDR.bands * DATA.HDR.samples * DATA.HDR.lines, 1);
            counter = 1;
            for j = 1:DATA.HDR.lines
                for i = 1:DATA.HDR.bands
                    d(((counter-1)*DATA.HDR.samples+1):counter*DATA.HDR.samples) = DATA.I(j, :, i)';
                    counter = counter + 1;
                end
            end
            fwrite(fp, d, word_name);
        otherwise
            %   default mode: BSQ
            d = zeros(DATA.HDR.samples, DATA.HDR.lines, DATA.HDR.bands);
            for i = 1:DATA.HDR.bands
                d(:, :, i) = DATA.I(:, :, i)';
            end
            fwrite(fp, d, word_name);
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     Finished writing out data file
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fclose(fp);
    
%   end of file.
end
