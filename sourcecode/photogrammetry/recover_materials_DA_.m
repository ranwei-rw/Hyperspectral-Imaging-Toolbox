% Syntax
%     [MAT, MATIND, MATAB] = recover_materials_DA_(S, debug, drate, max_clu_num, tmax, tmin, cool_rate, split_threshold) 
%   
%   Clustering a given set of feature vectors using the deterministic
%   annealing approach introduced by Kenneth Rose, Deterministic Annealing
%   for Clustering, Compression, Classification, Regression and Related
%   Optimization Problems, Proceeding IEEE, 1998. 
%    
% Input:
%   S:               3D matrix of source image reflectance in vector format(height by width by bands) 
%   debug:           The level of debugging information to be shown. It
%                    ranges from 1 to 3, 1 (default) for the minimal information
%                    while 5 means the most detailed description. 
%   tmax:            The maximum temperature of the deterministic annealing process. (default: 0.02)
%   tmin:            The minimum temperature of the deterministic annealing process (default 0.00025)
%   cool_rate:       The cooling factor in each iteration of the outer loop. (default: 0.8)
%   max_clu_num:     The maximum number of clusters. (default: 20)
%   split_threshold: The (dot product) threshold below which a cluster should be split. When the dot
%                    product between the real and the surrogate centroid vectors fall below this
%                    threshold, then the cluster needs to be split. (default: cos(5*pi/180))
%   drate:           the rate that the input image is down sampled at. Default is 1 (no down-sampling)
%   
%   Distance metric used in this implementation is the spectral angle
%   (or the dot product) between any two normalised reflectance spectra.
%
% Output:
%   MAT:    materials found (materials x bands) 
%   MATAB:  The material abundancy matrix for each pixel (height x width x materials number) 
%   MATIND: The material index image for each pixel (height x width)
%   
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.6
% Last Update Date: 28 Aug 2014

function [MAT, MATIND, MATAB] = recover_materials_DA_(S, debug, drate, max_clu_num, tmax, tmin, cool_rate, split_threshold) 
    
    %   handle and verify input parameters.    
    if ~exist('split_threshold', 'var') || numel(split_threshold) ~= 1 || split_threshold(1) < 0 
       split_threshold = cos(5*pi/180);
    end

    
    if ~exist('cool_rate','var') || numel(cool_rate) ~= 1 || cool_rate(1) < 0 
       cool_rate = 0.8;
    end
    
    if ~exist('tmin','var') || numel(tmin) ~= 1 || tmin(1) < 0 
       tmin = 0.00025;
    end
    
    if ~exist('tmax','var') || numel(tmax) ~= 1 || tmax(1) < 0 
       tmax = 0.02;
    end
    
        
    if ~exist('max_clu_num','var') || numel(max_clu_num) ~= 1 || max_clu_num(1) < 0 
       max_clu_num = 20;
    end
    
    if ~exist('drate', 'var') || numel(drate) ~= 1 || drate(1) < 0 
       drate = 1;
    end
    
    if ~exist('debug', 'var') || numel(debug) ~= 1 || debug(1) < 0 
       debug = 0;
    end
    
    if debug > 1
        disp('Start material unmixing using Deterministic Annealing.');
    end
    %   firstly, check whether the input reflectance has already been set to column format:
    %   height*width by bands
        
    [height, width, bands] = size(S);
    if debug >= 2
        s = sprintf('Input image size: height: %d, width: %d, bands: %d', height, width, bands);
        disp(s);
    end
    
    if drate > 1
        if debug >= 2
            disp('Image is down sampled by %d.', drate);
        end
        % Down-sample to avoid out of memory error.
        DS = zeros(ceil(height/drate), ceil(width/drate), bands);
        for b = 1:bands
            DS(:, :, b) = downsample(downsample(S(:, :, b), drate)', drate)';
        end

        S = DS;

        [height, width, bands] = size(S);
        if debug >= 2
            s = sprintf('New Image size: height: %d, width: %d, bands: %d', height, width, bands);
            disp(s);
        end
    end
    
    S2D = reshape(S, [height*width, bands])';

    %Clear unused variables
    clear DS S;
    %Continue...
    [feature_num, spectra_num] = size(S2D);
    if debug >= 2
        %   output parameters used by this function call
        s = sprintf('Maximum Temperature: %f\nMinimum (Stopping) Temperature: %f\Maximum Material Number: %d\n', ...
                     tmax, tmin, max_clu_num);
        disp(s);
        s = sprintf('Cooling Down Rate: %f\nMaterial clusters split threshold: %f', cool_rate, split_threshold);
        disp(s);
    end
    t = tmax;
    k = 2;
    % number of clusters = 1, We associate each cluster with a pair of centroid vectors, one is the
    % real center, the other is a perturbed version of it. The distance measure between them is a
    % measure of how cohesive the cluster is. When the distance between these two vectors exceeds a
    % threshold, then the cluster needs to be split. Note that the real cluster centroids and
    % association probabilities need to be averaged between each of these vector pairs.
    
    DELTA_Y = 0.1 * std(S2D, 0, 2);
    Y = zeros(feature_num, k);

    Y(:, 1) = median(S2D, 2);

    Y(:, 2) = Y(:, 1) + DELTA_Y;
    
    % normalize each column vector in Y
    CMAGNITUDE = sqrt(sum(Y.*Y, 1));
    CMAGNITUDE(CMAGNITUDE == 0) = 1;    % avoid divisions by zero
    Y                 = Y ./ CMAGNITUDE(ones(feature_num, 1), :);    
    CWEIGHT           = ones(k, 1)/k; % the weight p(y_{j}) (or prior probability) of each cluster.
    ASSOC             = ones(spectra_num, k)/k; % initialise p(y|x)
    assoc_thresh      = 0.005; % half a percent.
    y_change_thresh   = 0.001; % no more than 0.001 in each feature dimension
    break_flag        = 0; % if this flag is set to 1, then the outer loop will
    max_inner_iter    = 20; % maximum allowable number of iterations of the inner loop.
    
    if debug >= 1
        disp('Deterministic Annealing Starts');
    end
    
    while (1)
        if debug >= 2
            s = sprintf('Current temperature %g, Number of clusters %g.', t, k/2);
            disp(s);
        end
        if k >= 2*max_clu_num
            if debug >= 2
                disp('Maximum number of clusters is reached.');
            end
            break_flag = 1;
        end
        if t * cool_rate < tmin
            % Now optimise just the distortion term
            % sum_{x, y} p(x, y) d(x, y) by setting t = 0
            if debug >= 2
                disp('Current temperature below min threshold.');
            end
            break_flag = 1;
        end
        inner_iter_num = 0;
        while inner_iter_num < max_inner_iter   % Inner loop
            
            inner_iter_num = inner_iter_num + 1;
            if inner_iter_num > 1 && break_flag ~= 1
                PRE_ASSOC = sparse(ASSOC);
                PRE_Y = Y;
            end
            
            % subtract by the smallest distance to avoid very small values of ASSOC leading to divide by zero error while normalizing
            % ASSOC across clusters. 
            MIN_DIS    = min(1 - (abs(S2D' * Y)), [], 2); 
            % unnormalised association
            ASSOC = CWEIGHT(:, ones(spectra_num, 1))' .* exp(-(1 - (abs(S2D' * Y)) - MIN_DIS(:, ones(k, 1)))/t);
            m = sum(ASSOC == 0);
            if m ~= 0
                ASSOC = sparse(ASSOC);
            end
            
            ASSOC_SUM  = sum(ASSOC, 2);
            % normalisation constants
            if m ~= 0
                ASSOC = sparse(ASSOC ./ ASSOC_SUM(:, ones(k, 1)));    % normalised association..
            else
                ASSOC = ASSOC ./ ASSOC_SUM(:, ones(k, 1));
            end               
            
            CWEIGHT    = sum(ASSOC, 1)'/spectra_num;                      % marginal probability or weight of each cluster.
            Y          = S2D * ASSOC/spectra_num ./ CWEIGHT(:, ones(feature_num, 1))';
            CMAGNITUDE = sqrt(sum(Y .^ 2, 1));                            % centroid magnitudes
            Y          = Y ./ CMAGNITUDE(ones(feature_num, 1), :);        % normalise each column vector in Y
            
            % Termination condition of the inner loop.
            if inner_iter_num > 1
                assocChange = mean(abs(PRE_ASSOC - ASSOC));
                YChange     = mean(abs(PRE_Y - Y));
                if mean(assocChange) < assoc_thresh && mean(YChange) < y_change_thresh
                     break;
                end
            end
            if break_flag == 1
                break;
            end
        end
        
        if break_flag == 1
            break;
        end
        
        % Check for cluster split.
        pre_k = k;
        for i = 1:2:pre_k
            if (k < 2 * max_clu_num && Y(:, i)' * Y(:, i + 1) < split_threshold) % if use dot produce distance metric.
                if debug >= 2
                    s = sprintf('New cluster splits from %d.', floor(i/2) + 1);
                    disp(s);
                end
                Y(:, k + 1) = Y(:, i + 1);      %   Create a new pair of cluster centroids
                Y(:, k + 2) = Y(:, i + 1) + DELTA_Y;
                Y(:, i+ 1)  = Y(:, i) + DELTA_Y; %   Reassign the perturbed version to the real center of the current cluster.
                Y(:, k+ 2)  = Y(:, k+ 2) / norm(Y(:, k+ 2), 2);     % normalise new centroid vectors
                Y(:, i+ 1)  = Y(:, i+ 1) / norm(Y(:, i+ 1), 2);
                CWEIGHT(k + 1, 1) = CWEIGHT(i + 1, 1)/2;            % Update the cluster weights.
                CWEIGHT(k + 2, 1) = CWEIGHT(i + 1, 1)/2;
                CWEIGHT(i, 1)     = CWEIGHT(i, 1)/2;
                CWEIGHT(i + 1, 1) = CWEIGHT(i, 1);
                k = k + 2;
            end
        end
        t = cool_rate * t;
    end
    
    % Now assign the "real" centroids.
    MAT         = Y(:, 1:2:k)';
    MATAB       = ASSOC(:, 1:2:k) + ASSOC(:, 2:2:k);
    [~, MATIND] = max(MATAB, [], 2);
    MATIND      = reshape(MATIND, [height, width]);
    mat_num = size(MAT, 1);
    TEMP_AB = zeros(height, width, mat_num);
    for i = 1:mat_num
        TEMP_AB(:, :, i) = reshape(MATAB(:, i), [height width]);
    end
    MATAB = TEMP_AB;
    clear TEMP_AB; 
    disp('Material unmixing using Deterministic Annealing is finished.');
    if debug >= 2
        s = sprintf('%d materials are found.', mat_num);
        disp(s);
    end
end