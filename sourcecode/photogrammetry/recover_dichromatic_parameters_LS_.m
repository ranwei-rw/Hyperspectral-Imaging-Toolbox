% Syntax:
%   [K, G, S] = recover_dichromatic_parameters_LS(I, L, debug, neighb_size, gray_threshold)
%   
%   Description:
%   This function removes specularities from a hyperspectral image. All the parameters below
%   are defined in the context of the dichromatic reflection model.
% 
% Input: 
% 
%   I: the radiance image (stored as a 3D array with size height x width x bands).
%   L: the illuminant power spectrum. It could be either a 1D vector of
%      size bands x 1, or a 2D or 3D pixel-wise matrix.
%
%   neighb_size (optional): The size of the neighbourhood used for the recovery of the
%       reflectance. The default value is neighb_size = 5;
%   gray_threshold (optional): a threshold used to determine whether a material is a shade of gray. If
%       the reflectance spectra within a cluster does not deviate from a uniform spectrum (a flat
%       line) beyond this threshold, then we will assume that it is a shade of gray and purely
%       diffuse. (default value is 2) gray_threshold = 2;
%   
%   debug: Defines the level of displaying debugging information. Default is 1, the least
%       information will be given
% 
% Output: 
%
%   K: the map of specular coefficients at the image pixels, stored as a 2D array of size height x
%      width.  
%   S: the normalised reflectance cube per pixel and wavelength, stored as a 3D array of size height
%      x width x bands. The reflectance spectrum at each pixel is normalised to a unit L2-norm. That
%      is the vector S(i, j, :) is normalised.
%       
%   G: the shading factor g per pixel. G is stored as a 2D array with a size of height x width. 
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly.
% Version: 1.0.9
% Last Update Date: 28 May 2015
%   Update: speed improvement for about 1/57%

% Version: 1.0.8
% Last Update Date: 13 Oct 2014

function [K, G, S] = recover_dichromatic_parameters_LS_(I, L, debug, neighb_size, gray_threshold)
    
    %   Verify input parameters.
    if ~exist('debug', 'var') || numel(debug) ~= 1 || debug(1) < 0 
       debug = 1;
    end
    
    if ~exist('neighb_size', 'var') || numel(neighb_size) ~= 1 || neighb_size(1) < 0 
       neighb_size = 5;
    end
    
    if ~exist('gray_threshold', 'var') || numel(gray_threshold) ~= 1 || gray_threshold(1) < 0 
       gray_threshold = 2;
    end

    %   Remove the first few bands of a visible spectral image (where the signal is weak and can
    %   induce noisy segmentation results). Do this only when the band number is greatert than 15
    if debug >= 1
        disp('Start to compute specularities.');
    end
    
    if debug >= 3
        s = sprintf('Parameters input: ignore Threshold: %d, gray Threshold: %d.', neighb_size, gray_threshold);
        disp(s);
    end
    
    [height, width, bands] = size(I);
    
    if debug >= 3
        s = sprintf('Image dimensions: height = %d, width = %d, bands = %d.', height, width, bands);
        disp(s);
    end

    %   get the matrix R = I/L
    R = zeros(height, width, bands);
    
    % Here the illuminant estimate could be either a column vector or a
    % pixel-wise 2D or 3D matrix
    [lh, lw, ln] = size(L);
    
    %   to check the dimensions of L
    
    if ln ~= 1
        %   3d
        diml = 3;
        if (lh ~= 1 && lh ~= height) || (lw ~= 1 && lw ~= width)
            %   L has to be a 3D matrix with either 1 value per band or a
            %   full dimension pixel-wise matrix
            error('Dimensions of L is not correct');
        end
    elseif lh ~= 1 && lw ~= 1
        %   2D matrix, same for all bands
        if lh ~= height || lw ~= width
            error('Dimensions of L is not correct.');
        end
        diml = 2;
    else
        %   vecotr
        if lh == 1 && lw == 1
            error('Dimensions of L is not correct. A vector or matrix is required');
        end
        diml = 1;
    end
    
    switch diml 
        case 3
            for b = 1:bands  
                R(:, :, b) = I(:, :, b) ./ L(:, :, b);
            end
        case 2
            for b = 1:bands  
                R(:, :, b) = I(:, :, b) ./ L;
            end
        case 1
            for b = 1:bands  
                if L(b) > 0
                    R(:, :, b) = I(:, :, b) ./ L(b);
                else
                    R(:, :, b) = I(:, :, b);
                end
            end
        otherwise
            error('Dimensions of L is not correct');
    end

 
    NormLambda = sqrt(sum(sqrt(sum(I, 1).^2), 2));
    bandWeight = reshape(NormLambda, 1, bands)/sqrt(sum(NormLambda, 3).^2);
    bandWeight = bandWeight';
             
    %   Normalise R before clustering.
    if debug >= 2
        disp('Normalising R.');
    end
    
    Rnorm = sqrt(sum(R.^2, 3));
    normalisedR = zeros(height, width, bands);
    
    for b = 1:bands        
        normalisedR(:, :, b) = R(:, :, b) ./ Rnorm;
    end
%     
    %   Convert the 3D reflectance cube(size height x width x bands) into a 2D reflectance array    
    %   (size height X width X bands) as input for DA clustering.    
    R2D = zeros(bands, height * width);
    for b = 1:bands
        R2D(b, :) = reshape(normalisedR(:, :, b), [height * width 1]);
    end
        
    if debug >= 2
        disp('Start to compute K. Please be patient.');
    end
    
    G = zeros(height, width);
    K = zeros(height, width);
    S = zeros(height, width, bands);

    dotnum = 9;
    if debug >= 3
        disp('Computing');
    end
    
    %Set optimisation to echo-off
    options = optimset('Display', 'off');
    for i = neighb_size+1:height-neighb_size-1
        for j = neighb_size+1:width-neighb_size-1
            %Display debug information
            if debug >= 3
                if dotnum == 80
                    dotnum = 0;
                    fprintf('\n');
                end
                fprintf('.');
                dotnum = dotnum + 1;
            end

            % The mean reflectance of the material in the neighbourhood
            clusterS = reshape(mean(mean(R(i-neighb_size:i+neighb_size, j-neighb_size:j+neighb_size, :), 1), 2), [bands 1]);            
         
            % Normalise the neighbourhood reflectance.
            clusterS = clusterS / norm(clusterS);
            
            % Do the optimisation        
            A = [clusterS .* bandWeight bandWeight]; % the left hand side of the problem.        
            RIJ = reshape(R(i, j, :), [bands 1]);
            x = lsqnonneg(A, RIJ.* bandWeight);
            minr = min(RIJ);
            maxr = max(RIJ./clusterS);
            if x(2) > minr

                if maxr ~= 0 && minr ~= 0
                    x = lsqlin(A, RIJ.* bandWeight, [], [], [], [], [0, 0], [maxr, minr], x, options);
                else
                    x = [0, 0];
                end
            end

            G(i, j) = x(1);
            K(i, j) = x(2);
        end
    end
    
    if debug >= 3
        disp('\nFinished computing K');
    end   
    
    if debug >= 2
        disp('Computing S.');
    end
     
    % Note that R2D has been normalised so that each column has a unit L2-norm.
    % Angle (in degrees) with the uniform spectrum ones(bands, 1)
    angle_R       = acos(sum(R2D, 1) / sqrt(bands))/pi * 180; 
    grayMask1D    = logical(angle_R < gray_threshold); 
    grayMask2D    = reshape(grayMask1D, [height, width]); % obtain a 2D map of gray pixels.
    K(grayMask2D) = 0;

    % G is the L2-norm of the original reflectance spectra when K = 0
    G(grayMask2D) = Rnorm(grayMask2D);   
    
    % Assign reflectance to each pixel once G and K are known.
    for b = 1:bands
        S(:, :, b) = (R(:, :, b) - K) ./ G;
    end
    S(isnan(S)) = 0; % fix divide by zero error when G = 0.
    S(isinf(S)) = 0; % fix divide by zero error when G = 0.
    if debug >= 2
        disp('Finished computing S.');
    end

    if debug >= 2
        disp('Normalising the reflectance and recoverying the appropriate factors.');
    end
    
    Rnorm = sqrt(sum(R.^2, 3));
    S = S./Rnorm(:, :, ones(bands, 1));
    S(isnan(S)) = 0; % fix divide by zero error when G = 0.
    S(isinf(S)) = 0; % fix divide by zero error when G = 0.
    %Scale the shading factor accordingly
    G = G.*Rnorm;
    if debug >= 2
        disp('Dichromatic parameters are successully recovered. Exiting.');
    end
    
%   end of function: recover_dichromatic_parameters_LS_
end
