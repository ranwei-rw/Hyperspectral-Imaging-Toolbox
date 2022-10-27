%   function [MAT, MATIND, MATAB] = recover_materials_KM(S, max_clu_num, debug, drate) 
%   
%   Clustering a given set of feature vectors using k-means. 
%    
% Input:
%
%   S:              3D matrix of source image reflectance in vector format
%                    (height by width by bands) 
%   debug:           The level of debugging information to be shown. It
%                    ranges from 1 to 5, 1 (default) for the minimal information
%                    while 5 means the most detailed description. 
%   max_clu_num:     The maximum number of clusters. (default: 20)
%   drate:  the rate that the input image is down sampled at. Default is 1 (no downsampling)
%                    (default: 1)
%   
% Output:
%
%   MAT:    materials found by this function. (materials x bands) 
%   MATAB:  The material abundancy matrix for each pixel (height x width x materials number) 
%   MATIND: The material index image for each pixel (height x width)
%   
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.6
% Last Update Date: 28 Aug 2014


function [MAT, MATIND, MATAB] = recover_materials_KM_(S, max_clu_num, debug, drate) 
    
    %   handle and verify input parameters.
    disp('Start material recovery using k-means.');
    
    if ~exist('max_clu_num', 'var') || numel(max_clu_num) ~= 1 || max_clu_num(1) < 0 
       max_clu_num = 20;
    end
    
    if ~exist('debug', 'var') || numel(debug) ~= 1 || debug(1) < 0 
       debug = 1;
    end
    
    if ~exist('drate', 'var') || numel(drate) ~= 1 || drate(1) < 0
       drate = 5;
    end
   
    if debug >= 2
        disp('Reshaping input reflectance.');
    end
        
    [h, w, bands] = size(S);
    if debug >= 2
        s = sprintf('Image size: height: %d, width: %d, bands: %d', h, w, bands);
        disp(s);
    end
    temp_data = S;    
    if drate > 1
        if debug >= 2
            s = sprintf('Image is down sampled by %d.', drate);
            disp(s);
        end
        % Down-sample to avoid out of memory error.
        DS = zeros(ceil(h/drate), ceil(w/drate), bands);
        for b = 1:bands
            temp = temp_data(:, :, b);
            DS(:, :, b) = downsample(downsample(temp, drate)', drate)';
        end
        % clear I;
        clear temp;
        temp_data = DS;
        clear DS;
        [h, w, bands] = size(temp_data);
        if debug >= 2
            s = sprintf('New Image size: height: %d, width: %d, bands: %d', h, w, bands);
            disp(s);
        end
    end
    
    %Reshape the data for the k-means
    S2D = reshape(temp_data, [h * w, bands]);
    clear temp_data;
    
    if debug >= 2
        disp('Reshaping reflectance finished.');
    end    
    
    %Perform the k-means
    [~, MAT, ~, ~] = kmeans(S2D, max_clu_num);
    
    %Build the matrix abundances by using the distances. This is a form of
    %expectation over a Gaussian distribution about the centers.
    if debug >= 2
        disp('Computing the material abundance matrix....');
    end
    
    [h, w, bands] = size(S);
    S2D           = reshape(S, [w*h , bands]);
    D             = dist(S2D, MAT');
    xm            = min(min(D(find(D>0)))); %Avoid divisions by zero
    D(find(D==0)) = xm;
    Scale         = sum(ones(w*h, max_clu_num)./D, 2);
    MATAB         = (ones(w*h, max_clu_num)./D)./Scale(:, ones(max_clu_num, 1));
    Scale         = 1:max_clu_num;
    MATIND        = Scale(ones(w*h, 1), :);
    MATIND        = reshape(MATIND, h, w, max_clu_num); %  Reshape the matrices and finish off
    MATAB         = reshape(MATAB, h, w, max_clu_num);
    disp('material recovering using K-Means is finished');
end
    