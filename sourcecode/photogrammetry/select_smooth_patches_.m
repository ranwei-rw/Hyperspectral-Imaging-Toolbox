% function [PATCHES, ...
%           LIGTH_EST_MASK, ...
%           CONTRAST_MASK, ...
%           HIGHLIGHT_MASK, ...
%           PATCH_MAP] = select_smooth_patches(I, ...
%                                              pheight, ...
%                                              pwidth, ...
%                                              mean_thresh, ...
%                                              selected_num, ...
%                                              debug)                                                 
%
% Select smooth patches with the most contrast (averaged over all bands) due
% to specularity (not because the patch spans different materials).
%
% Input:
%
%   I:              height x width x bands - the dimensions of input data
%   pheight:        height of each patch, default to 20.
%   pwidth:         width of each patch, default to 20.
%   selected_num:   number of selected patches which is default to 50
%   mean_thresh:    only choose patches with mean radiance across pixels and bands greater than this 
%                   threshold. Default to 100
%   debug:          debug information will shown; values range from 0 to 2
%
% Output:
%
%   LIGTH_EST_MASK: a mask showing 1 at the pixels selected for initial light
%                   estimation, 0 otherwise.
%   CONTRAST_MASK:  1 where the most contrast patches are selected, 0 otherwise (same size as
%                   image)
%   HIGHLIGHT_MASK: the mask corresponding to patches where specularities
%                   need to be removed.
%   PATCH_MAP:      if a patch with coordinates i,j is selected, then patchMap(i, j) = 1.
%   PATCHES:        selected patches.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.7
% Last Update Date: 26 May 2015

%   increased speed

% Version: 1.0.6
% Last Update Date: 25 Aug 2014

function [PATCHES, ...
          LIGTH_EST_MASK, ...
          CONTRAST_MASK, ...
          HIGHLIGHT_MASK, ...
          PATCH_MAP] = select_smooth_patches_(I, ...
                                              pheight, ...
                                              pwidth, ...
                                              mean_thresh, ...
                                              selected_num, ...
                                              debug)

    %   get size information of the input I
    [height, width, bands] = size(I);
    
    IMEAN = mean(I, 3);
    
    % Only select unsaturated patches.
    
    if ~exist('debug', 'var')
        debug = 0;
    end
    
    if ~exist('selected_num', 'var')
        selected_num = 50;
    end
    
    if ~exist('mean_thresh', 'var')
        mean_thresh = 100;
    end
    
    if ~exist('pwidth', 'var')
        pwidth = 20;
    end
    
    
    if ~exist('pheight', 'var')
        pheight = 20;
    end
    
    %   Ge the max value of data
    maxdata = max(I(:));
    depth = ceil(log2(maxdata));
    saturating = 2^depth*0.95;
    
    if debug >= 2
        s = sprintf('Patch height is %d, Patch width is %d, Mean Threshold is %d, Data depth is %d', ...
            pheight, pwidth,  depth);
        disp(s);
    end
    
    HIGHLIGHT_MASK = zeros(height, width);
    
    homo_count = 0;          % count the number of homogeneous patches (with uniform spectral reflectance)
    
    contrast_count = 0;      % count the number of homogeneous patches with high contrast.
    
    cosine_angle_thresh = cos(25/180 * pi);
    proportion_thresh  = 0.70;
    
	pnum_h = floor(height/pheight);
    pnum_w = floor(width/pwidth);
    miny = zeros(pnum_h, pnum_w);
    minx = zeros(pnum_h, pnum_w);
    
    for i = 1 : pheight : pnum_h*pheight
        for j = 1 : pwidth : pnum_w*pwidth
            ph = min(i+pheight-1, height);
            pw = min(j+pwidth-1, width);
            IPATCH3D = I(i : ph, j : pw, :);
            
            % Don't select patches containing saturated pixels.
            if max(IPATCH3D(:)) >= saturating
                continue;
            end

            patch_pixel_num = (ph - i + 1) * (pw  - j + 1);
            IMEAN_PATCH_1D  = reshape(IMEAN(i:ph, j:pw), [patch_pixel_num, 1]);
            
            % If this is a homogeneous surface, it is a candidate for
            % highlight removal. Combine the edge map with the deviation
            % angle from the dichromatic plane to determine if a patch is homogeneous.
            if mean(IMEAN_PATCH_1D) > mean_thresh
                I2D = zeros(bands, patch_pixel_num);
                for k = 1:bands
                    I2D(k, :) = reshape(IPATCH3D(:, :, k), [patch_pixel_num, 1]);
                end
                [U, ~, ~] = svd(I2D);
                A         = U(:, 1:2);
                P         = A /(A' * A)*A';
                I2DProj   = P * I2D;
                
                angle_cosine = sum(I2DProj .* I2D, 1) ./ (sqrt(sum(I2D .* I2D, 1)) .* sqrt(sum(I2DProj .* I2DProj, 1)));
                pixel_count  = sum(angle_cosine > cosine_angle_thresh);    % number of pixel forming an angle of less
                % than devAngleThresh with the dichromatic plane.
                
                if pixel_count >= proportion_thresh * patch_pixel_num
                    HIGHLIGHT_MASK(i:min(i+pheight-1, height), j:min(j+pwidth-1, width)) = 1;
                    if debug ~= 0
                        homo_count = homo_count + 1;
                    end
                    contrast_count = contrast_count + 1;
                    % Store the standard deviation and the coordinates of the left top corner
                    % of each patch.
                    miny(contrast_count, 1) = (i-1)/pheight+1;
                    minx(contrast_count, 1) = (j-1)/pwidth+1;
                end
                
            end
        end
    end
    

    if debug ~= 0
        s = sprintf('%d smooth patches are detected', homo_count);
        disp(s);
    end
    
    % Select patches for light estimation.
    LIGTH_EST_MASK = zeros(height, width);
    for i = 1:contrast_count
        r = miny(i);
        c = minx(i);
        LIGTH_EST_MASK((r-1)*pheight+1:min(r*pheight, height), (c-1)*pwidth+1:min(c*pwidth, width)) = 1;
    end
    
    if contrast_count < selected_num
        selected_num = contrast_count;  
    end
    
    if selected_num == 0
        warning('Scylluras:photogrammetry:select_smooth_patches', ...
                'Could not find valid patches');
    end
    
    % Compute the contrast mask
    CONTRAST_MASK = zeros(height, width);
    PATCH_MAP = zeros(ceil(height/pheight), ceil(width/pwidth));
    
    X = zeros(selected_num, 1);
    Y = zeros(selected_num, 1);
    
    for i = 1:selected_num
        r = miny(i);
        c = minx(i);
        PATCH_MAP(r, c) = 1;
        Y(i, 1) = (r-1)*pheight+1;
        X(i, 1) = (c-1)*pwidth+1;
        CONTRAST_MASK((r-1)*pheight+1:min(r*pheight, height), (c-1)*pwidth+1:min(c*pwidth, width)) = 1;
    end
  
    PATCHES = ones(selected_num, 4);
    PATCHES(:, 1) = Y;
    PATCHES(:, 2) = X;
    PATCHES(:, 3) = PATCHES(:, 3)*pheight;
    PATCHES(:, 4) = PATCHES(:, 4)*pwidth;
    clear I;
end