% Function to estimates the light source of an input_image. Depending on
% the parameters the estimation is equal to Grey-Wolrd, Max-RGB, general
% Grey-World, Shades-of-Gray or Grey-Edge algorithm. 
%
% SYNOPSIS
%   [L, C] = estimate_lightsource_(I, order, norm, MASK);
%
% INPUT :
%   I:      Hyperspectral image data cube
%	order:  the order of differentiation (ranges from 0 to 2). Default to
%           0. Users can also use method names, ie 'grey_world', 'maxrgb',
%           'shades_of_grey' or 'grey_edge' instead of actual number for
%           order. 
%	norm:   minkowski norm used (if norm ==  -1 means infinity). Default to 2
%   sigma:  value of sigma, default to 0
%   MASK:   Binary images with zeros on image borders should be ignored. If
%           MASK is missing, the default value will be generated according
%           to the value of order. 
% OUTPUT: 
%   L:      illuminant estimation in the form of a row vector
%   C:      color corrected image
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2015 All Rights Reserved.
% Author: Ran Wei


%   this function was based on general_cc.m with the following
%   explaination.
% general_cc: estimates the light source of an input_image. 
%
% Depending on the parameters the estimation is equal to Grey-Wolrd, Max-RGB, general Grey-World,
% Shades-of-Gray or Grey-Edge algorithm.
%
% SYNOPSIS:
%    [white_R ,white_G ,white_B,output_data] = general_cc(input_data,njet,mink_norm,sigma,mask_im)
%    
% INPUT :
%   input_data    : color input image (NxMx3)
%	njet          : the order of differentiation (range from 0-2). 
%	mink_norm     : minkowski norm used (if mink_norm==-1 then the max
%                   operation is applied which is equal to minkowski_norm=infinity).
%   mask_im       : binary images with zeros on image positions which
%                   should be ignored.
% OUTPUT: 
%   [white_R,white_G,white_B]           : illuminant color estimation
%   output_data                         : color corrected image

% LITERATURE :
%
% J. van de Weijer, Th. Gevers, A. Gijsenij
% "Edge-Based Color Constancy"
% IEEE Trans. Image Processing, accepted 2007.
%
% The paper includes references to other Color Constancy algorithms
% included in general_cc.m such as Grey-World, and max-RGB, and
% Shades-of-Gray.

function [L, C] = estimate_lightsource_(I, order, norm, sigma, MASK)

    if ~exist('sigma', 'var')
        sigma = 1;
    end

    if ~exist('norm', 'var')
        norm = 2;
    end

    if ~exist('order', 'var')
        order = 2;
    end

    if ischar(order)
        switch order
            case 'grey_world'
            case 'maxrgb'
            case 'shades_of_grey'
                order = 1;
            case 'grey_edge'
                order = 2;
            otherwise
                error('unknown method name');
        end
    end
    
    [height, width, bands] = size(I);

    if ~exist('MASK', 'var')
        MASK = zeros(height, width);
    else
        if height ~= size(MASK, 1) || width ~= size(MASK, 2)
            error('Incorrect dimensions of MASK');
        end
    end

    %   finished the replacement of function norm_derivative.m in the
    %   following section of about 15 lines below
    if order == 0
        for i = 1:bands
            R = I(:, :, i);
            I(:, :, i) = get_2d_gaussian_derivative_(R, sigma, 0, 0);
        end
    elseif order == 1
        for i = 1:bands
            R = I(:, :, i);
            I(:, :, i) = sqrt(get_2d_gaussian_derivative_(R, sigma, 1, 0).^2 + get_2d_gaussian_derivative_(R, sigma, 0, 1).^2);
        end
    elseif order < 0
        error('incorrect value of order');
    else
        %computes frobius norm
        for i = 1:bands
            R = I(:, :, i);
            I(:, :, i) = sqrt(get_2d_gaussian_derivative_(R, sigma, 2, 0).^2 + get_2d_gaussian_derivative_(R, sigma, 0, 2).^2 + get_2d_gaussian_derivative_(R, sigma, 1, 1).^2);
        end
    end

    C = I;
    I = abs(I);
    
    %   make a border according to the value of the fourth argument sigma

    if order == 0
        MASK(2:height-1, 2:width-1) = 1;
    else
        MASK(4:height-3, 4:width-3) = 1;
    end
    % %   ignore saturated pixel handling because it's probably not
    % %   realistic for hyperspectral images which normally have 12bits or
    % %   more depth
    L = zeros(bands, 1);

    % minkowski norm = (1, infinity]
    if norm ==  -1
        %minkowski-norm is infinit: Max-algorithm
        for i = 1:bands
            M    = I(:, :, i);
            L(i) = max(M(:).*MASK(:));
        end
    else
        K = power(I, norm);
        for i = 1:bands
            M    = K(:, :, i).*MASK;
            L(i) = power(sum(M(:)), 1/norm);
        end
    end
    if sum(L) ~= 0
        L = L./sqrt(sum(L.^2));
        for i = 1:bands
            C(:, :, i) = C(:, :, i)/(L(i)*sqrt(bands));
        end
        L = L(~isnan(sum(L, 2)), :)';
    end
    
end    




