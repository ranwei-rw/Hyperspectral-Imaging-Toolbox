%%  Syntax:
%   K = recover_specularity_KL(I, L, cluster_num, debug)
%   
%%  Description:
%   Recover the specularity map from a single hyperspectral image. The method used
%   in this function was described in DICTA Paper `Specularity Removal from Imaging Spectroscopy
%   Data via Entropy Minimisation' by Lin Gu and Antonio Robles-Kelly
%      
%%  Input: 
%
%       I:           Input image cube, can be in either multi-spectrum or RGB image
%       L:           Light spectrum vector.
%       cluster_num: The number of clusters used for the K-means algorithm
%       debug:       Debugging information display level (1 - 3). 
%
%%  Output: 
%
%       K:           the map of specular coefficients at the image pixels,
%                    stored as a 2D array of size height x width.  
%
%
%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2014 All Rights Reserved.
% Author: Ran Wei, Lin Gu and Antonio Robles-Kelly. 
% Version: 1.0.6
% Last Update Date: 1 Sep 2014

function K = recover_specularity_KL_(I, L, cluster_num, debug)

    if ~exist('debug', 'var')
        debug = 1;
    end
    
    if ~exist('cluster_num', 'var')
        cluster_num = 20;
    end
    
    if debug >= 3
        disp(strcat(num2str(now), ' - start to recover specularity'));
    end

    %   firstly, get size information of input data matrix I
    [height, width, bands] = size(I);
    if debug >= 2
        s = sprintf('Image dimensions are: Height: %d, Width: %d, Bands: %d', height, width, bands);
        disp(s);
    end

    %   Define R as the output of I/L, firstly, define R as a 2D matrix of height*width by bands
    %   also, put the reshaped original data into it
    R = reshape(I, height*width, bands);

    %   Get the value of R = I/L (L is LIGHT), which is done band by band
    for i = 1 : bands
        R(:, i) = R(:, i) / L(i);    
    end

   
    D      = R;                             %   define a variable for illuminant
    R_bar  = mean(R, 2);                    %   define R bar
    R      = R - R_bar(:, ones(1, bands));  %   get R - R_bar, stored in R
    R_std  = std(R, 1, 2);                  %   get the standard deviation of R-Rbar
    R_mask = ones(height*width, 1);         %   define a mask to remove certain noise
    R_mask(R_std < (mean(R_std)/100.0)) = 0;

    %   get the noise pixels removed from both R_bar and R_standard
    R_log = logical(R_mask);
    R     = R(R_log, :);
    R_std = R_std(R_log);
    R     = R ./ R_std(:, ones(bands, 1));
    D     = D(R_log, :);

    %   define K in vector form of (height*width by 1)
    K_1D = zeros(height * width, 1);
    K_S  = zeros(length(R_std), 1);
    %   get the kmeans for the input points where CLUSTERS is a height*width by 1 vector with values
    %   ranging from 1 to cluster_num (or less, depends on how many clusters were actually categorised)
    if debug >= 2
        disp('Computing K Means for I/L');
    end

    %   doing some downsampling or rescaling to get a proper input parameter for kmeans
    downsample = 1;
    if height <= cluster_num
        cluster_num = height;
        downsample = 0;
    end
    samples = ones(cluster_num, bands);
    if downsample == 1
        %   get the R downsampled to cluster_num height
        down_num = fix(height/cluster_num);
        for i = 1:cluster_num
            samples(i, :) = R(1+i*down_num, :);
        end
    else
        samples = R(1:cluster_num, :);
    end
    [CLUSTERS, ~] = kmeans(R, cluster_num, 'Start', samples, 'distanc', 'correlation', 'emptyaction', 'singleton');

    %   Go through each cluster to get 
    cluster_num = max(CLUSTERS);

    for i = 1:cluster_num
        mask_clusters = CLUSTERS;
        mask_clusters(CLUSTERS == i) = 1;
        mask_clusters(CLUSTERS ~= i) = 0;
        diffuse_patch      = D(logical(mask_clusters), :);
        diffuse_patch_mean = mean(diffuse_patch, 2);
        diffuse_patch_AC   = diffuse_patch - diffuse_patch_mean(:, ones(bands, 1));
        %   then use PCA to get the main axis
        [coeff, ~, ~] = princomp(diffuse_patch_AC);
 
        p = coeff(:, 1);
        p = (p - mean(p))/std(p);
        x = [p ones(bands, 1)];
        diffuse_num = size(diffuse_patch, 1);
        patch_para = zeros(2, diffuse_num);
        for j = 1:diffuse_num
            patch_para(:, j) = x \ diffuse_patch(j, :)';
        end

        % use Finalyson's entropy minimisation is applied
        %   define entropy vector (E) for 180 degrees [1 180]
        E     = zeros(180, 1);
        theta = 0.001;
        for k = 1:180     
            theta       = theta + pi/180.0;                                 %   get the value of entropy for each angle
            triang      = [sin(theta) cos(theta)];                          %   get its cos and sin
            d           = triang * patch_para;                              %   d = a*cos() + b*sin()
            unit_length = (max(d) - min(d))/64;                             %   get its distance
            bins        =  min(d) : unit_length : max(d);                   %   get the hist bin ranges
            pos         = (histc(d, bins) + 1) / (length(d) + length(bins));%	get the possibility of each p(i|theta)
            E(k)        = - sum(pos .* log(pos));
        end
        %   find the index of the value which has the minimum entropy among 180 angles
        %   since E is a vector, minE is then a one element vector, or, a scalar
        [~, minE] = min(E);
        %   then we can get the angle which leads to the minimal entropy
        ratio     = abs(tan(pi * minE / 180.0));
        if debug == 3
            s = sprintf('Best ratio of angle found: %g pi for cluster %d', minE/180.0, i);
            disp(s);
        end

        patch_temp = patch_para;
        for k = 1:100
            pca_ratio = mean(patch_temp(2, :) ./ patch_temp(1, :));
            k_patch   = patch_temp(2, :) - pca_ratio * patch_para(1, :);
            k_patch(k_patch < 0) = 0;
            patch_temp(2, :)     = patch_temp(2, :) - k_patch;
            if sum(k_patch)/diffuse_num < 0.001
                break;
            end
        end

        K_S(logical(mask_clusters)) = patch_para(2, :) - ratio * patch_para(1, :);

    end
    K_1D(R_log) = K_S;
    %Use the absolute value for the intersect
    K = (reshape(K_1D, height, width));

    if debug >= 2
        disp('Leaving recover specularity function using entropy minimisation');
    end
    if debug == 3
        disp(strcat(num2str(now), ' - recovering specularity is finished'));
    end

%   end of function
end