% Usage:
%     export_hsz_(filename, HSZ);
%     export_hsz_(filename, HSZ, options);
% 
% Writes a compressed NICTA HSPipeline file. This is an HDF5 file with all
% the variables recovered by the HSPipeline routine.
% Inputs:
%     filename: The name of the file (including the path) to be written to disk
%     HSZ: The data structure delivered by the HSPipeline routine
%     options: Write to disk options these are
%         'compression': Level of compression for the HDF5 data (0-9). 
%             The default is 9 (maximum compression).
%         'datatype': Type for the data written to disk on the HDF5 datasets. 
%             The default is 'uint16', but 'uint8' can also be used. 
%             
% Example:
%      options.compression = 5;
%      options.datatype = 'uint8';
%      export_hsz_('test.hsz', HSPipeline);
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.6
% Last Update Date: 18 Sep 2014

function export_hsz_(filename, HSZ, options)

    if ~exist('options', 'var') || ...
       ~isfield(options, 'compression') || ...
        options.compression > 9 || ...
        options.compression < 1
            options.compression = 9;            %This is the default, i.e. max compression
    end
    if ~isfield(options, 'datatype') || ...
       (~strcmpi(options.datatype, 'double') && ...
       ~strcmpi(options.datatype, 'uint8') &&...
       ~strcmpi(options.datatype, 'uint16')) 
            options.datatype = 'uint16';       %Encode the materials using double measurments
    end

    %Create a file with standard options

    %   note: here we use low level H5 functions to create a HDF5 file and
    %   close it.
    H5FileID = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5F.close(H5FileID);

    %   in the rest part of this function, high level functions are used to
    %   write datasets and attributes to complete
    names = fieldnames(HSZ);
    %Start with all the numeric values (data)
    for i=1:length(names)
        if isstruct(HSZ.(names{i}))
              subnames = fieldnames(HSZ.(names{i}));
              for j=1:length(subnames)
                  if isnumeric(HSZ.(names{i}).(subnames{j}))      
                      %Create the dataset
                      %Find the number of dimensions for the array
                      n = ndims(HSZ.(names{i}).(subnames{j}));
                      %Go on to find the scaling factors
                      m = max(HSZ.(names{i}).(subnames{j}));
                      q = min(HSZ.(names{i}).(subnames{j}));
                      for k=2:n
                          m = max(m);
                          q = min(q);
                      end
                      %Do the actual write-to-disk operation
                      if strcmpi(options.datatype, 'uint8')
                          if  ~isempty(strfind(subnames{j}, 'Indexes'))  || (m-q)==0
                              if floor(HSZ.(names{i}).(subnames{j}))>2^8
                                  t=uint16(floor(HSZ.(names{i}).(subnames{j})));
                                  h5create(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'DATA'), ...
                                    size(t), 'Datatype', 'uint16', 'ChunkSize', size(t), 'Deflate', ...
                                    options.compression);
                              else
                                  t=uint8(floor(HSZ.(names{i}).(subnames{j})));
                                  h5create(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'DATA'), ...
                                    size(t), 'Datatype', 'uint8', 'ChunkSize', size(t), 'Deflate', ...
                                    options.compression);
                              end
                          else
                              t = uint8(floor((HSZ.(names{i}).(subnames{j})-q)/(m-q)*2^8));
                              h5create(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'DATA'), ...
                               size(t), 'Datatype', 'uint8', 'ChunkSize', size(t), 'Deflate', ...
                               options.compression);
                          end
                      end
                      if strcmpi(options.datatype, 'uint16')
                          if  ~isempty(strfind(subnames{j}, 'Indexes')) || (m-q)==0
                              t = uint16(floor(HSZ.(names{i}).(subnames{j})));
                          else
                              t = uint16(floor(double(HSZ.(names{i}).(subnames{j})-q)/double(m-q)*(2^16-1)));
                          end
                          h5create(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'DATA'), ...
                               size(t), 'Datatype', 'uint16', 'ChunkSize', size(t), 'Deflate', ...
                               options.compression);
                      end
                      %Write it to the file
                      h5write(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'DATA'), t);
                      %Write the scaling of the data
                      h5create(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'MAX'), ...
                               size(m), 'Datatype', 'double', 'ChunkSize', size(m), 'Deflate', ...
                               options.compression);
                      h5write(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'MAX'), m);
                      h5create(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'MIN'), ...
                               size(q), 'Datatype', 'double', 'ChunkSize', size(q), 'Deflate', ...
                               options.compression);
                      h5write(filename, strcat('/', names{i}, '/', subnames{j}, '/', 'MIN'), q);    
                  end
              end
        end
    end
    %Do the header strings (attributes)
    subnames = fieldnames(HSZ.HDR);
    for j=1:length(subnames)
        if ~isnumeric(HSZ.HDR.(subnames{j}))
            h5writeatt(filename, '/HDR', subnames{j}, HSZ.HDR.(subnames{j}));
        end
    end
    %Do the Endmemember header strings (attributes)
    if isfield(HSZ, 'HDREndmembersS')
        subnames = fieldnames(HSZ.HDREndmembersS);
        for j=1:length(subnames)
            if ~isnumeric(HSZ.HDREndmembersS.(subnames{j}))
                h5writeatt(filename, '/HDREndmembersS', subnames{j}, HSZ.HDREndmembersS.(subnames{j}));
            end
        end
    end
    %Do the Illuminant header strings (attributes)
    if isfield(HSZ, 'HDREndmembersL')
        subnames = fieldnames(HSZ.HDREndmembersL);
        for j=1:length(subnames)
            if ~isnumeric(HSZ.HDREndmembersL.(subnames{j}))
                h5writeatt(filename, '/HDREndmembersL', subnames{j}, HSZ.HDREndmembersL.(subnames{j}));
            end
        end
    end
%   end of function
end
