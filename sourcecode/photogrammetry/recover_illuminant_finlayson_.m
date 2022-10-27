%% Syntax:
%        L = recover_illuminant_finlayson(I, PATCHES, MASK, debug)
%
% Estimate the light spectrum using the dichromatic plane + optimisation
% method proposed in the following paper
% 
% Graham D. Finlayson, Gerald Schaefer: 
% Convex and Non-convex Illuminant Constraints for Dichromatic Colour Constancy. 
% CVPR (1) 2001.
% 
% Input:
%
%   I:       height x width x bands, the input spectral radiance image.
%   PATCHES: selected patches over the given data cube
%   MASK:    masking matrix, if not given, it equals a no masking
%   debug:   Whether to show debug information during processing. 1 for
%            minimal and 3 for maximal
% 
% Output:
%   
%   L: the estimate of the light spectrum.
%
%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei, Antonio Robles-Kelly and Cong Phuoc Huynh
% Version: 1.0.8
% Last Update Date: 17 Nov 2014

% Version: 1.0.8
% Default illuminant (all 1) values will be returned with warning message if no
% valid illuminant is found.

% Version:1.0.7: Removed possible patches whose edges are outside the given
% image.

function L = recover_illuminant_finlayson_(I, PATCHES, MASK, debug)
    
    if ~exist('debug', 'var')
        debug = 0;
    end
    
    use_mask = 1;
    if ~exist('MASK', 'var')
        use_mask = 0;
    end
    
    % Image dimension
    [image_height, image_width, bands] = size(I);
    
    %   warning message if given data is of invalid size
    if image_height == 0 || image_width == 0 || bands == 0
        %   stop and exit
        error('Illuminate Recovery using Finlayson method: size of parameter I is incorrect.');
    end

    if debug > 1
        s = sprintf('image information: height = %d, width = %d, bands = %d', image_height, image_width, bands);
        disp(s);
    end
        
    %Mask the image accordingly
    if use_mask == 1
        for i = 1:bands
            I(:, :, i) = I(:, :, i).*MASK;
        end
    end

    %   chech whether input patches is pre-selected
    [~, p_w] = size(PATCHES);
    select_num = 0;
    pre_selected_patched = 1;
    
    if p_w == 4
        %   which means the selected patches are ready. do not call select_smooth_patches_ function
        if debug == 3
            disp('Pre-set Patches are ready.');
        end
    elseif p_w == 1
        if PATCHES == 0
            if debug == 3
                disp('No patches are given, default 50 patches will be chosen.');
            end
            select_num = 50;
            pre_selected_patched = 0;
        else
            if debug == 3
                s = sprintf('%d patches will be chosen', select_num);
                disp(s);
            end
            select_num = PATCHES;
            pre_selected_patched = 0;
        end
    else
        warning('Scyllarus:Photogrammetry:recover_illuminant_finlayson', ...
                'Error in geometries of selected smooth patches. Could not proceed. A all-ones illuminant estimation is returned');
        L = ones([bands, 1]);
        return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               %
    %          STEP 1               %
    % Smooth the image in each band %
    % using the Wiener filter.      %
    %                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for b = 1:bands
        I(:, :, b) = wiener2(I(:, :, b), [5 5]);    
    end

    if debug > 1
        disp('Input image is smoothed by Wiener filter');
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                       %
    %          STEP 2                       %
    %   Generate Patches required which     %
    %   selects smooth homogeneous patches  %
    %                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if pre_selected_patched == 0
        mean_threshold = 100;
        patchheight = 20;
        patchwidth  = 20;
        if debug > 2
            disp('start to choose smooth patches.');
        end
        
        PATCHES = select_smooth_patches_(I, ...
                                         patchheight, ...
                                         patchwidth, ...
                                         mean_threshold, ...
                                         select_num, ...
                                         debug);

        [h, w] = size(PATCHES);

        if h==0 || w==0
            warning('Scyllarus:Photogrammetry:recover_illuminant_finlayson', ...
                    'No valid patches found - from Function: recover_illuminant. Exiting with a all-one illuminant');
            L = ones([bands, 1]);
            return;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                               %
        %                               %
        %   Variable declaration area   %
        %                               %
        %                               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if debug >= 2
            disp('patches are chosen.');
        end
    else
        if debug >= 2
            disp('smooth patches are given. Skipping patch selection.');
        end
        %   generate lightMask, Select patches for light estimation.
        L_EST_MASK = zeros(image_height, image_width);
        for i = 1 : select_num
            L_EST_MASK(PATCHES(i, 1)+1 : PATCHES(i, 1)+ PATCHES(i, 3), PATCHES(i, 2)+1 : PATCHES(i, 2)+ PATCHES(i, 4)) = 1;
        end
    end
    
   
    %   patchheight: the height of each patch.
    %   patchwidth:  the width of each patch.
    %   patch_Y: (N x 1) the y coordinates of the top left corners of the patches.
    %   patch_X: (N x 1) the x coordinates of the top left corners of the patches.
    %       where N is the number of patches.
    patchheight = PATCHES(1, 3);
    patchwidth  = PATCHES(1, 4);
    patch_Y     = PATCHES(:, 1);
    patch_X     = PATCHES(:, 2);
    
    % check how many patches are given (counting the number of topleft corner
    % y coordinate)
    patch_num = size(patch_Y, 1);
    
    if debug >= 2
        s = sprintf('%d patches are selected\n', patch_num);
        disp(s);
    end
    
    patch_pixels_ = patchheight * patchwidth; % number of pixels per patch.

    PATCHES = zeros(patch_num, bands, patch_pixels_);
    valid_patches = 0;
    for i = 1:patch_num
        y1 = patch_Y(i);
        y2 = min(y1+patchheight-1, image_height);
        x1 = patch_X(i);
        x2 = min(x1+patchwidth-1, image_width);
        if y2 <= image_height || x2 <= image_width
            valid_patches = valid_patches + 1;
            for b = 1:bands
                PATCHES(valid_patches, b, :) = reshape(I(y1: y2, x1: x2, b), [patch_pixels_ 1]);
            end
        end
    end

    if valid_patches < patch_num
        PATCHES = PATCHES(1:valid_patches, bands, patch_pixels_);
        if debug > 1
            s = sprintf('%d patches were removed. %d patches are selected', patch_num - valid_patches, valid_patches);
            disp(s);
        end
    end
    
    % This method is attributed to Finlayson and Schaefer
    I1 = zeros(bands, patch_pixels_);
    M = zeros(bands, bands); % sum of projection vector
    for i = 1:valid_patches
        I1(:, :) = PATCHES(i, :, :);
        [U1, ~, ~] = svd(I1);
        A = U1(:, 1:2);  
        M = M + A * inv(A' * A) * A';
    end

    [V, ~] = eig(M);

    % Choose an eigen vector to be the light spectrum.
    % Do this by choosing the vector with all positive or negative elements
    nCols = size(V, 2);
    L = ones(bands, 1)*-1;
    for i = 1:nCols
        eigV = real(V(:, i));
        if (abs(sum(sign(eigV))) == bands)
            L = abs(eigV);
            break;
        end
    end
    
    if isequal(L, ones(bands, 1)*-1)
        %   now valid illuminant estimation was found. Warning users
        L = ones(bands, 1);
        warning('Scyllarus:Photogrammetry:recover_illuminant_finlayson', ...
                'No valid illuminant was found. A all-ones illuminant estimation is returned');
    end
    

end