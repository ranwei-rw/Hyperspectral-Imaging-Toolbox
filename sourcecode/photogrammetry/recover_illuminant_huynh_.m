% Syntax:
%
%       L = recover_illuminant_huynh(I, alpha, PATCHES, debug)
%
% Input:
%
%       I:        hyperspectral image stored as a 3D array.
%       alpha:    the variable used to test the result. default to 50
%       PATCHES:  number of PATCHES used. default to 20;
%       debug:    Level of debug information to be shown during processing. 0(by default) for minimum and 3 is maximum
%
% Output:
%
%       L: the illuminate computed, which is a 2D array of size (bands x 1), where nBands is the number of bands,
%          storinglight spectral power. 
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014-2015 All Rights Reserved.
% Author: Ran Wei and Cong Phuoc Huynh
% Version: 1.0.9
% Last Update Date: 19 Jun 2015
%   speed up by about 6%

% Version: 1.0.8
% Last Update Date: 26 Aug 2014

%   version 1.0.7
%   Cancelled input argument bitdepth which is not used in the code
%
%   version 1.0.6
%   Cancelled input argument MASK. It is less useful to put a masking matrix here. If peopel need, they can do masking
%   explicitly for I. Also, cancelled temp variable in calling wiener2 function to save some time
%
function L = recover_illuminant_huynh_(I, alpha, PATCHES, debug)

    if ~exist('debug', 'var')
        debug = 1;
    end
    
    if ~exist('PATCHES', 'var')
        PATCHES = 0;
    end
    
    if ~exist('alpha', 'var')
        alpha = 50;
    end
    
    if isinteger(I)
        I = double(I);
    end
    % Image dimension
    [height, width, bands] = size(I);

    %   warning message if given data is of invalid size
    if height == 0 || width == 0 || bands == 0
        %   stop and exit
        error('Illuminate Recovery: Data size is incorrect.');
    end

    if debug > 1
        disp('image information:');
        s = sprintf('height = %d, width = %d, bands = %d.', height, width, bands);
        disp(s);
    end
    
    %   check whether input PATCHES is pre-selected
    [~, p_w] = size(PATCHES);
    max_patch_num = 0;
    pre_selected_patched = 1;
    if p_w == 4
        %   which means the selected patches are ready. do not call select_smooth_patches_ function
        if debug > 1
            disp('Pre-set Patches are ready.');
        end
    elseif p_w == 1
        if PATCHES == 0
            max_patch_num = 50;
        else
            max_patch_num = PATCHES;
        end
        
        if debug > 1
            s = sprintf('%d Patches will be automatically chosen.', max_patch_num);
            disp(s);
        end
        pre_selected_patched = 0;
    else
        error('Error in geometries of selected smooth PATCHES');
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
                    %          STEP 1                       %
                    %   Generate Patches required which     %
                    %   selects smooth homogeneous patches  %
                    %                                       %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if debug > 1
        disp('start to choose smooth patches.');
    end
    
    if pre_selected_patched == 0
        ImeanThresh = 100;
        patchheight = 20;
        patchwidth  = 20;
        
        [PATCHES, L_EST_MASK] = select_smooth_patches_(I, ...
                                                       patchheight, ...
                                                       patchwidth, ...
                                                       ImeanThresh, ...
                                                       max_patch_num, ...
                                                       debug);

        [h, w] = size(PATCHES);

        if h==0 || w==0
            warning('Scyllarus:photogrammetry:recover_illuminant_huynh', ...
                    'No valid patches found - from Function: recover_illuminant. Exiting');
            L = 0;
            return;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                               %
        %                               %
        %   Variable declaration area   %
        %                               %
        %                               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if debug > 1
            disp('patches are chosen.');
        end
    else
        if debug > 1
            disp('smooth patches are given. Skipping patch selection.');
        end
        %   generate lightMask Select patches for light estimation.
        L_EST_MASK = zeros(height, width);
        for i = 1:max_patch_num
            L_EST_MASK(PATCHES(i, 1)+1 : PATCHES(i, 1)+ PATCHES(i, 3), PATCHES(i, 2)+1 : PATCHES(i, 2)+ PATCHES(i, 4)) = 1;
        end
    end

    LMAX = zeros(bands, 1);
    for b = 1:bands
        LMAX(b) = max(max(I(:, :, b) .* L_EST_MASK));
    end

    %   allocate memory for illuminates
    L = ones(bands, 1) * mean(LMAX);

    % initial value of the objective function, set to be infinity. 
    previous_average_cost = 1e+20; 
    PATCH_Y                = PATCHES(:, 1);
    PATCH_X                = PATCHES(:, 2);
    patchheight            = PATCHES(1, 3);
    patchwidth             = PATCHES(1, 4);
    patch_number           = size(PATCHES, 1);%   check how many patches are given (counting the number of topleft corner y coordinate)

    if debug > 1
        s = sprintf('%d patches are selected', patch_number);
        disp(s);
    end

    %   find how many pixels is in a patch. Assume that all patches are in the same size
    patch_pixels      = patchwidth * patchheight;
    PATCH_G           = zeros(patch_number, patchheight, patchwidth);
    PATCH_K           = zeros(patch_number, patchheight, patchwidth);
    PATCH_G_SQUARE    = zeros(patch_number, patchheight, patchwidth);
    PATCH_G_SUM       = zeros(patch_number, 1);
    PATCH_K_SUM       = zeros(patch_number, 1);
    PATCH_GK_SUM      = zeros(patch_number, 1);
    DICHROMATIC       = zeros(patch_number, 1);% == 1 if the corresponding p^{th} patch is used to find the optimal L, 0 otherwise.
    current_iteration = 1;
    max_iterations    = 100;      %   maximum number of iterations of the big while loop.
    TEMP_L            = zeros(bands, 1);
    G_OUT             = zeros(patchheight, patchwidth);
    K_OUT             = zeros(patchheight, patchwidth);
    best_cost_value   = 1e+20;
    BEST_L            = TEMP_L;
    D_OUT             = zeros(bands, 1); % temporary variable for diffuse radiance component.
    PATCH_D           = zeros(patch_number, bands); % diffuse radiance component of the patches.

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                               %
                    %                               %
                    %   Actual computation starts   %
                    %                               %
                    %                               %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if debug > 1 
        disp('solving iteration begins:');
    end

    % store the L in each iteration.
    while current_iteration < max_iterations
        if debug > 1
            s = sprintf('Iteration: %d', current_iteration);
            disp(s);
        end
        current_iteration = current_iteration + 1;
        %   ignore the iteration information

        % Recover G_OUT, K_OUT, S_out and adjust the radiance for each patch to conform to the DICHROMATIC model.    
        total_dich_error   = 0;
        total_smooth_error = 0;
        
        for p = 1 : patch_number
            %   compute the geometry of current patch
            y1 = PATCH_Y(p);
            y2 = min(y1 + patchheight - 1, height);
            x1 = PATCH_X(p);
            x2 = min(x1 + patchwidth - 1, width);

            if debug > 2
                if current_iteration == 1
                    s = sprintf('Patch geometry: (%d, %d, %d, %d)\n', y1, x1, y2, x2);
                    disp(s)
                end
            end
            %   get subset data of this patch, this patch is of n bands 
            current_patch_ = I(y1 : y2, x1 : x2, :);

            %   get DICHROMATIC model for this pathch
            [G_OUT, ...
             K_OUT, ...
             ~, ...
             p_dich_error, ...
             p_smooth_error, ...
             is_dichromatic] = patch_dichromatic_decompose_(current_patch_, ...
                                                             L, ...
                                                             alpha, ...
                                                             debug);

            total_dich_error   = total_dich_error + p_dich_error * is_dichromatic;
            total_smooth_error = total_smooth_error + p_smooth_error* is_dichromatic;
            PATCH_G(p, :, :)   = G_OUT;
            PATCH_K(p, :, :)   = K_OUT;
            PATCH_G_SQUARE(p, :, :) = G_OUT.*G_OUT;
            PATCH_G_SUM(p)     = sum(sum(PATCH_G_SQUARE(p, :, :)));
            PATCH_K_SUM(p)     = sum(sum(K_OUT.*K_OUT));
            PATCH_GK_SUM(p)    = sum(sum(G_OUT.*K_OUT));
            DICHROMATIC(p, 1)  = is_dichromatic;
        end   

        % Check if the totalCost has actually decreased for the choice of G_OUT, K_OUT, 
        % S_out, L. We need to do this because we are not sure whether constraining
        % G_OUT, K_OUT, S_out, L to the DICHROMATIC plane (obtained by SVD) will limit the
        % solution space.
        % totalCost = total_dich_error + total_smooth_error; % take the average DICHROMATIC error per patch

        PS = PATCH_G_SUM(~isnan(PATCH_G_SUM(p))) ~= 0;
        temp1 = sum(PATCH_K_SUM);
        temp2 = sum(PATCH_GK_SUM(PS).*PATCH_GK_SUM(PS)./PATCH_G_SUM(PS));

        for b = 1 : bands
            temp3 = 0;
            temp4 = 0;

            for p = 1:patch_number
                y1 = PATCH_Y(p);
                y2 = min(y1 + patchheight-1, height);
                x1 = PATCH_X(p);
                x2 = min(x1 + patchwidth-1, width);

                temp3 = temp3 + sum(sum(I(y1 : y2, x1 : x2, b) .* reshape(PATCH_K(p, :, :), y2-y1+1, x2-x1+1)));
                
                if PATCH_G_SUM(p) ~=0 && ~isnan(PATCH_G_SUM(p))
                    temp4 = temp4 + PATCH_GK_SUM(p) * sum(sum(I(y1 : y2, x1 : x2, b) .* reshape(PATCH_G(p, :, :), y2-y1+1, x2-x1+1))) / PATCH_G_SUM(p);
                end
            end
            
            TEMP_L(b, 1) = (temp3 - temp4)/(temp1 - temp2);
        end

        for p = 1 : patch_number
            y1 = PATCH_Y(p);
            y2 = min(y1 + patchheight-1, height);
            x1 = PATCH_X(p);
            x2 = min(x1 + patchwidth-1, width);
           

            for b = 1:bands
                D_OUT(b, 1) = sum(sum(I(y1:y2, x1:x2, b) .* reshape(PATCH_G(p, :, :), y2-y1+1, x2-x1+1) - reshape(PATCH_G_SQUARE(p, :, :), y2-y1+1, x2-x1+1) * TEMP_L(b, 1)));
            end
            D_OUT = D_OUT / PATCH_G_SUM(p);
            
            PATCH_D(p, :) = D_OUT;

        end

        % Compute the new DICHROMATIC error
        I_STACK = zeros(bands, patch_number * patch_pixels); 
        I_DICH_STACK = zeros(bands, patch_number * patch_pixels);
        
        for p = 1:patch_number
            y1 = PATCH_Y(p);
            y2 = min(y1 + patchheight-1, height);
            x1 = PATCH_X(p);
            x2 = min(x1 + patchwidth-1, width);
            for b = 1:bands
                I_STACK(b, (p-1)*patch_pixels + 1:p*patch_pixels) = reshape(I(y1:y2, x1:x2, b) * DICHROMATIC(p, 1), [patch_pixels 1]);
                I_DICH_STACK(b, (p-1)*patch_pixels + 1:p*patch_pixels) = ...
                    DICHROMATIC(p, 1) * reshape(PATCH_D(p, b) * PATCH_G(p, :, :) + PATCH_K(p, :, :) * TEMP_L(b), [patch_pixels 1]);
            end
        end

        temp_dich_error = 0;
        for b = 1:bands
            I_DIFF = reshape(I_STACK(b, :) - I_DICH_STACK(b, :), [patch_number*patch_pixels  1]);
            temp_dich_error = temp_dich_error + sum(sum(I_DIFF .* I_DIFF));
        end


        % The total cost at this step of each iteration is always lower than 
        % that after the previous step, where we optimise the cost with respect 
        % to PATCH_G, patch_s_, and PATCH_K, keeping L fixed.
        average_cost = (temp_dich_error + total_smooth_error) / sum(DICHROMATIC);  % take the average DICHROMATIC error per patch

        if best_cost_value > average_cost
             BEST_L = TEMP_L;
             best_cost_value = average_cost;
        end

        %    if angle(L, TEMP_L) < 0.01 || (previous_average_cost - average_cost) / previous_average_cost < 0.01
        %   this is what Cong was using. But, when the right side is negative, it's definitely < 0.01. Which
        %   may not be a convergence condition.
        %   instead, I use the abs value of difference. Or we can even use the angle difference alone
        % Check the angle difference between the L in the previous and the
        % current iteration.
        % If the angle difference decreases => sign of convergence, 
        % otherwise: sign of divergence.

        if angle(L, TEMP_L) < 0.01 || abs((previous_average_cost - average_cost) / previous_average_cost) < 0.01
            break;
        else
            % update L
            L = TEMP_L;
            previous_average_cost = average_cost;
        end

        L = abs(L/max(abs(L))) * max(LMAX);
    end

    if debug > 1
        disp('solving iteration ends');
    end

    % If this process doesn't converge after a maximum allowable number of
    % iterations (max_iterations), we pick the one with lowest cost.
    % if ((current_iteration == max_iterations) && (LAngleDiff >= 0.05))
    if current_iteration == max_iterations
        L = BEST_L;
    end
    
%   end of function
end