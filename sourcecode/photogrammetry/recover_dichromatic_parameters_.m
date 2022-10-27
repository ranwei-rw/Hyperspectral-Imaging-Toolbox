%% Dichromatic parameter recovery
%
%%   Syntax
%       [k, g, K, S] = recover_dichromatic_parameters(I, L)
%       [k, g, K, S] = recover_dichromatic_parameters(I, L, options)
%   
%% Description:
%   This function recovers the photometric parameters in the context of the 
%   dichromatic reflection model whereby
%
%       I = L (g S+k K)
%      
%   where, the illuminant power spectrum is given on L, g is the shading
%   factor, S is the surface reflectance, k is the specular coefficient and
%   K is the specular highlight.
% 
%% Input: 
% 
%   I:       the radiance image (stored as a 3D array with size height x width x bands).
%   L:       the illuminant power spectrum, stored as a 1D vector of size (bands x 1). 
%
%   options: Structure with the following fields
%            method:     String specifying the method to be used. This can be 
%                        'LS' for the linear least squares solver, 
%                        'KL' for the energy minimization method, or 
%                        'TI' for the method of Tan and Ikeuchi neighbourhood size (optional, used by the LS option): a
%                             threshold used to ignore the processing of regions with a small number of pixels. Those
%                             with a number of pixels smaller than or equal to this are considered diffuse. (default
%                             value is 5) ignoreThresh = 5;  
%           gray_threshold:  (optional, used by the LS option): a threshold used to determine whether a material is a shade of
%                        gray. If the reflectance spectra within a cluster does not deviate from a uniform spectrum (a flat
%                        line) beyond this threshold, then we will assume that it is a shade of gray and purely diffuse.
%                        (default value is 2) gray_threshold = 2; 
%           numclusters: (option, used by the KL option): the number of clusters used for the K-means, default: 20.
%           debug:       Defines the level of displaying debugging information. Default is 1, the least information will
%                        be given 
% 
%% Output: 
%
%   K: The wavelength dependant specularity per pixel, stored as a 3D array of size height x
%      width x bands.  
%   S: The normalised reflectance cube per pixel and wavelength, stored as a 3D array of size height
%      x width x bands. The reflectance spectrum at each pixel is normalised to a unit L2-norm. That
%      is the vector S(i, j, :) is normalised.
%   k: The specular coefficient per pixel. k is stored as a 2D array with a size of height X width.   
%   g: The shading factor per pixel. g is stored as a 2D array with a size of height X width. 
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei, Antonio Robles-Kelly and Cong Phuoc Huynh
% Version: 1.0.9
% Last Update Date: 14 Oct 2014


function [k, g, K, S] = recover_dichromatic_parameters_(I, L, options)

%   Verify input parameters.
    if ~exist('options', 'var') || ~isfield(options, 'method')
       options.method = 'LS';        %Set the least squares method by default
    end
    
    if ~isfield(options, 'debug') || numel(options.debug) ~= 1 || options.debug(1) < 0 || options.debug(1) > 5
       options.debug = 1;
    end
    
    %   Setup the variable for the number of bands
    [rows, cols, bands] = size(I);
    
    switch upper(options.method)
    %   Do the least squares method if selected (default)    
        case 'LS'
            disp('Using least-squares method for the specularity recovery.');
        %   Verify input parameters.
            if ~exist('options.neigbourhoodsize', 'var') || ...
               numel(options.neigbourhoodsize) ~= 1 || ...
               options.neigbourhoodsize < 0 || ...
               options.neigbourhoodsize > rows || ...
               options.neigbourhoodsize > cols
                options.neigbourhoodsize = min([5, rows, cols]);
            end

            if ~exist('options.gray_threshold', 'var') || ...
               numel(options.gray_threshold) ~= 1 || ...
               options.gray_threshold(1) < 0 
                options.gray_threshold = 2;
            end
            
            [k, g, St] = recover_dichromatic_parameters_LS_(I, L, options.debug, options.neigbourhoodsize, options.gray_threshold);

            norm        = sqrt(sum(St.^2, 3));
            S           = St./norm(:, :, ones(bands, 1));
            S(isnan(S)) = 0;                            %   fix divide by zero error when G = 0.
            S(isinf(S)) = 0;                            %   fix divide by zero error when G = 0.
            g           = g.*norm;                      %   Recover the shading factor
            g(isnan(g)) = 0;                            %   fix divide by zero error when G = 0.
            g(isinf(g)) = 0;                            %   fix divide by zero error when G = 0.
            K           = ones(size(I))/sqrt(bands);    %   Recover the wavelength dependent specularity
            K(isnan(K)) = 0;                            %   fix divide by zero error when G = 0.
            K(isinf(K)) = 0;                            %   fix divide by zero error when G = 0.
            
        case 'KL'
            disp('Using the energy minimisation method for the specularity recovery.');
            if ~exist('options.numclusters', 'var') || numel(options.numclusters) ~= 1 || options.numclusters < 1
                options.numclusters = 20;   %Fix the default number of clusters
                disp('Using default number of clusters: 20');
            end
            %   Recover the specularity
            k           = recover_specularity_KL_(I, L, options.numclusters, options.debug);
            k(k<0)      = 0;%   Avoid negative values         
            T           = min(I, [], 3);%   Avoid invalid pixel   
            indx        = find(T-k < 0);
            k(indx)     = T(indx);            
            L           = to3D(L, [rows, cols, bands]);
            St          = I./L - k(:, :, ones(bands, 1));
            norm        = sqrt(sum(St.^2, 3));
            S           = St./norm(:, :, ones(bands, 1));
            S(isnan(S)) = 0;                            %   fix divide by zero error when G = 0.
            S(isinf(S)) = 0;                            %   fix divide by zero error when G = 0.
            g           = norm;                         %   Recover the shading factor
            g(isnan(g)) = 0;                            %   fix divide by zero error when G = 0.
            g(isinf(g)) = 0;                            %   fix divide by zero error when G = 0.
            K           = ones(size(I))/sqrt(bands);    %   Recover the wavelength dependent specularity
            K           = K./L;
            norm        = sqrt(sum(K.^2, 3)); 
            K           = K./norm(:, :, ones(bands, 1));
            K(isnan(K)) = 0;                            %   fix divide by zero error when G = 0.
            K(isinf(K)) = 0;                            %   fix divide by zero error when G = 0.
            
        case 'TI'
            %   Prepare for the directory checks
            pold = pwd;
            p = mfilename('fullpath');
            p = strtok(p, '\recover_dichromatic_parameters.m');
            cd(p);
            
            disp('Using the method of Tan and Ikeuchi for the specularity recovery.');
            
            if ~isequal(exist('TanRGBDespec.mexw64', 'file'), 3) && ~isequal(exist('TanRGBDespec.mexw32', 'file'), 3)
                %     tell users file name is incorrect
                disp('TanRGBDespec.mexw64 or TanRGBDespec.mexw32 does not exist. Now compiling it');
                cd('.\tan')
                mex -g TanRGBDespec.cpp zHighlightRemoval.cpp;
                if exist('.\tan\TanRGBDespec.mexw64', 'file') == 0 && exist('.\tan\TanRGBDespec.mexw32', 'file') == 0
                    disp('Could not compile TanRGBDespec. Exitng');
                    return;
                end
                if exist('.\tan\TanRGBDespec.mexw32', 'file') == 0
                    movefile('.\tan\TanRGBDespec.mexw64');
                    disp('TanRGBDespec is ready.');
                elseif exist('.\tan\TanRGBDespec.mexw64', 'file') == 0
                    movefile('.\tan\TanRGBDespec.mexw32');
                    disp('TanRGBDespec is ready.');
                end
                cd(pold);
            end
            
            K           = recover_specularity_tan(I);
            k           = sqrt(sum(K.^2, 3));
            L           = to3D(L, [rows, cols, bands]);
            St          = I./L - k(:, :, ones(bands, 1)).*K;    %   Recover the reflectance
            norm        = sqrt(sum(St.^2, 3));
            S           = St./norm(:, :, ones(bands, 1));
            S(isnan(S)) = 0;    %   fix divide by zero error when G = 0.
            S(isinf(S)) = 0;    %   fix divide by zero error when G = 0.            
            g           = norm; %   Recover the shading factor
        otherwise
            s = sprintf('Unknown option method: %s', options.method);
            disp(s);
    end
end

function L3D = to3D(L, size3D)

    if length(size3D) ~= 3
        error('Error in size information');
    end
    height = size3D(1);
    width  = size3D(2);
    bands  = size3D(3);
    
    [lh, lw, ln] = size(L);
    
    diml = 0;
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

    L3D = zeros(size3D);
    switch diml 
        case 3
            L3D = L;
        case 2
            for b = 1:bands  
                L3D = repmat(L, 1, 1, bands);
            end
        case 1
            L = reshape(L, [1, 1, bands]);
            L3D = repmat(L, height, width);
        otherwise
            error('Dimensions of L is not correct');
    end
    
%   end of function
end