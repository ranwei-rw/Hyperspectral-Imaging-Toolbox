%   This function is designed to generate release code for Scylluras Toolbox. 
%
%   What this function needs to do:
%       1. Copy all existing folders under src (except release and tests) if they are not existing, into Release folder
%       2. Copy all .m files into their corresponding folder under Release
%       3. Pcode all _.m files and put the generated .p file into corresponding folder under release
%       4. Check all other files (except itself and all sav files). If any
%          of them is missing from Release, copy it over from src folders. 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei 
% Version: 1.0.2
% Date: 26 Aug 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Version 1.0.2
%       Added modifying date check.
%
%   Version 1.0.1
%   Basic version

function generate_release(source, dest, level)
    
    if ~exist('level', 'var')
        level = 0;
    end
    
    if ~exist('source', 'var')
        source = './';
    end
    
    if ~exist('dest', 'var')
        dest = './Release/';
    end
    
    newlevel = level + 1;
    %   get the content of source
    contents = dir(source);
    %   analyse them.
    for i = 1:length(contents)
        m = contents(i);
        if strcmpi(m.name, '.') || ...
           strcmpi(m.name, '..') || ...
           strcmpi(m.name, 'html') || ...
           strcmpi(m.name, 'Release') || ...
           strcmpi(m.name, 'generate_release.m') || ...
           strcmpi(m.name, 'function_list.m') || ...
           strcmpi(m.name, 'tan') || ...
           strcmpi(m.name, 'SuiteSparse') || ...
           strcmpi(m.name, 'tests') || ...
           strcmpi(m.name, 'Publishing_Readme.txt')
            continue;
        end
        if m.isdir == 1
            % only pcode file on level 0 and 1
            %   iteratively go through
            nest_source = strcat(source, m.name, '/');
            nest_dest   = strcat(dest, m.name, '/');
            %   need to check whether this folder exists at target
            %   location. If not create it
            if ~exist(nest_dest, 'dir')
                mkdir(dest, m.name);
            end
            generate_release(nest_source, nest_dest, newlevel);
        else
            %   now deal with files
            [~, name, ext] = fileparts(m.name);
            if strcmpi(ext, '.asv')
                %   omit auto backup file
                continue;
            end
            if strcmpi(ext, '.m') && level <= 1
                if name(end) == '_'
                    %   this is a source file to be pcode-ed
                    pcode(strcat(source, m.name));
                    %   pfile will be generated under the folder where this
                    %   function is called first time
                    pfile = strcat(name, '.p');
                    [status, message] = movefile(pfile, dest);
                    if status == 0
                        error(message);
                    end
                else
                    % the interface m file, copy it over 
                    %   check modifying date
                    dest_name = strcat(dest, name, ext);
                    d = dir(dest_name);
                    if ~isempty(d)
                        if isequal(d.datenum, m.datenum)
                            continue;
                        end
                    end
                    [status, message] = copyfile(strcat(source, m.name), dest);
                    if status == 0
                        error(message);
                    end
                end
            else
                %   other cases, simply copy it over if it's been there
                %   already
                dest_name = strcat(dest, name, ext);
                if ~exist(dest_name, 'file')
                    %   a normal file, simply copy it over when the
                    %   destination is not there
                    copyfile(strcat(source, m.name), dest);
                end
            end
        end
    end
end
