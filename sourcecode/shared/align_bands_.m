%
% Function: 
%   I_aligned = align_bands(I, max_shift, window_size, crop_orig)
%
% Align all the band images of a hyperspectral image if there are misalignments. This function was originally developed
% for images taken by Fluxdata 3CCD 7Channel Cameras. These cameras can take hyperspectral images of 7 bands. Among these
% bands, there is certain random misalignment which may be resulted from CCD alignment error during manufacturing.
% 
% Input: 
%   I:           an input hyperspectral image.
%   max_shift:   maximal shift in pixel in both directions. Default to 20
%   window_size: the size of the neighborhood used to compute the focus value of every pixel. If window_size = 0, a
%                single focus measure is computed for the whole image. Default to 3.
%   crop_orig:   Whether to crop output image to cut off misalignment edges. Default to 1 meaning: yes
% Output: 
%   I_aligned:   a new hyperspectral image whose pixels in different band are aligned. If crop_orig == 1, depend on the
%                misalignment size, I_aligned will be cropped to smaller than original dimension in height and width
% 
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Cong Phuoc Huynh and Ran Wei

function I_aligned = align_bands_(I, max_shift, window_size, crop_orig)

    if ~exist('max_shift', 'var')
        max_shift = 20;
    end
    
    if ~exist('window_size', 'var')
        window_size = 3;
    end
    
    if ~exist('crop_orig', 'var')
        crop_orig = 1;
    end

    [height, width, band] = size(I);
    
    fm = zeros(band, 1);
    
    %   use the central part of the image as mask
    mask_rows = [round(3/8*height):round(5/8*height)]'; % rows of the mask
    mask_cols = [round(3/8*width) :round(5/8*width)]';  % cols of the mask
    
    n_mask_rows = numel(mask_rows);
    n_mask_cols = numel(mask_cols);
    grad_mask   = zeros(n_mask_rows, n_mask_cols, band);    
    
    for b = 1:band
        band_img = I(:, :, b);
        
        % Energy of laplacian (Subbarao92a)
        f_img = imfilter(band_img, fspecial('laplacian'), 'replicate', 'conv');
        if window_size == 0
            f_img = mean2(f_img.^2);
        else
            f_img = imfilter(f_img.^2, fspecial('average', [window_size window_size]), 'replicate');
        end
        
        fm(b) = mean(f_img(:));

        % compute the gradient mask
        [dx, dy] = gradient(band_img(mask_rows, mask_cols));
        grad_mag = dx.^2 + dy.^2;
        
        % retain the larger gradients
        grad_mask(:, :, b) = logical(grad_mag > median(grad_mag(:)));
    end

    % Select the band with the highest contrast as the reference band.
    [~, ref_band_idx] = max(fm);
    ref_grad_mask = grad_mask(:, :, ref_band_idx);
    
    % Initialise the aligned image
    I_aligned = zeros(height, width, band);
    I_aligned(:, :, ref_band_idx) = I(:, :, ref_band_idx);

    % The row and column indices to be searched for
    search_rows = n_mask_rows + [-max_shift:max_shift]';
    search_cols = n_mask_cols + [-max_shift:max_shift]';
    
    max_offset_x = 0;
    max_offset_y = 0;
    for b = 1:band        
        if (b == ref_band_idx)
            continue;
        end

        % Find the pixel shift distance by cross correlation
        [cc] = normxcorr2(ref_grad_mask, grad_mask(:, :, b));
                
        % Concentrate the search in the region in the middle of the correlation map assuming that the bands are only
        % shifted by a small number of pixels wrt each other.
        cc_search_region     = cc(search_rows, search_cols);        
        [~, max_idx]         = max(abs(cc_search_region(:)));
        [peak_row, peak_col] = ind2sub(size(cc_search_region), max_idx(1));
        
        % The shift (in the y and x directions) of pixels in the current band image (at band b) relative to the
        % reference band. 
        offset_y = peak_row - (max_shift + 1);
        offset_x = peak_col - (max_shift + 1);
        
        % Translating the current band in the opposite direction of (offset_y, offset_x) so that it matches the
        % reference band 
        x_start = 1 + max(0, -offset_x); 
        x_end   = width + min(0, -offset_x); 
        y_start = 1 + max(0, -offset_y); 
        y_end   = height + min(0, -offset_y);
        
        if offset_x > max_offset_x
            max_offset_x = offset_x;
        end
        
        if offset_y > max_offset_y
            max_offset_y = offset_y;
        end
        % Translate the current band image by (-offset_x, -offset_y) and crop out the irrelevant section
        I_aligned(y_start:y_end, x_start:x_end, b) = I([y_start:y_end] + offset_y, [x_start:x_end] + offset_x, b);
    end
    
    if crop_orig == 1 && (max_offset_x ~= 0 || max_offset_y ~= 0)
        I_aligned = I_aligned(max_offset_y+1:height-max_offset_y, max_offset_x+1:width-max_offset_x, :);
    end
    
%   end of function align_bands
end