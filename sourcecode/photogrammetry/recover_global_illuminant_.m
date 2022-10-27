%% Recover a single illuminant from a hyperspectral image
%
%% Syntax
%   L = recover_global_illuminant(I)
%   L = recover_global_illuminant(I, options)
% 
%% Description:
%   Recover a signle illuminant from a single hyperspectral image. 
%
%% Input:
%
%   I: hyperspectral image stored as a 3D array.
%   options: Structure with the following fields
%           bitdepth:      Is the data type for the spectral cube, i.e. number of bits per spectral measurement. By
%                          default this is 16. 
%           method:        Selects between the following methods 
%                          'HRK':    Employs the method of Huynh and Robles-Kelly (A Solution of the Dichromatic Model
%                                    for Multispectral Photometric Invariance, International Journal of Computer Vision
%                                    2010). This is the default choice.
%                          'FS':     Uses the method of Finlayson and Schaefer (Convex and Non-convex Illuminant 
%                                    Constraints for Dichromatic Colour Constancy, CVPR 2001).
%                          'GW':     Uses the Grey World method.
%                          'SG':     Uses the Shade of Grey method.
%                          'WP':     Uses the White Patch method.
%                          '1stOGE': Uses the 1st order Grey Edge method.
%                          '2ndOGE': Uses the 2nd order Grey Edge method.
%           drate:         Image downsampling rate for the Grey World, Shade of Grey, White Patch and Grey Edge methods.
%                          The default is 1, i.e. no downsampling.
%           order_of_grey: The order of the L^p mean used for the Shade of Grey method. The default is 1.
%           alpha:         The value for the regularisation term used for the HRK (Huynh and Robles-Kelly) method. The
%                          default for this is 50.  
%           patches:       Pre-selected patches. This could be a set of geometry data of patches with a format of
%                          (Top_left_y, top_left_x, height, width). This can be left empty.
%           debug:         Defines the level of debugging information shown at execusion time (DEBUG < 3). the default
%                          is 0 (least).  
%
%% Output:
%
%   L: a 2D array of size (1 x bands), where bands is the number of wavelength
%           indexed bands in the input image.
%% See also:
%   recover_multi_illuminant
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly. 
% Version: 1.0.7
% Last Update Date: 21 Aug 2014

%%
function L = recover_global_illuminant_(I, options)


    %Commence by checkin those variables applicable to both methods
    if ~exist('options', 'var') || ...
       ~isfield(options, 'method') || ...
       (~strcmpi(options.method, 'FS') && ~strcmpi(options.method, 'GW') && ...
        ~strcmpi(options.method, 'SG') && ~strcmpi(options.method, 'WP') && ...
        ~strcmpi(options.method, '1stOGE') && ~strcmpi(options.method, '2ndOGE'))
        %   default method is HRK
        options.method = 'HRK';
    end
    
    if ~isfield(options, 'debug') || isempty(options.debug) || options.debug > 3 || options.debug < 0
        options.debug = 0;
    end
    
    if options.debug > 0
        disp('No valid options.method is found. Using default method: HRK');
    end
    
    if ~isfield(options, 'bitdepth') ||...
       (options.bitdepth ~= 8 && options.bitdepth ~= 12) || ...
       numel(options.bitdepth) ~= 1
        options.bitdepth = 16;    
    end
    
    if ~isfield(options, 'patches') || isempty(options.patches)
        options.patches = 20;    
    end
    
    %Commence with the processing
    switch upper(options.method)
        case 'HRK'
            if options.debug > 0
                disp('Using the method of Huynh and Robles-Kelly for recoverying the illuminant power spectrum.');
            end
            %Use the method of Huynh and Robles-Kelly
            if ~isfield(options, 'alpha') || options.alpha<0 || numel(options.alpha) ~= 1
                options.alpha = 50;
            end
            L = recover_illuminant_huynh_(I, options.alpha, options.patches);
            if std(L) == 0 || isnan(sum(L))
                options.method = 'GW';
                if options.debug > 0
                    disp('The method of Huynh and Robles-kelly failed to converge... using Grey World instead.');
                end
                L = recover_illuminant_other_methods(I, options);
            end
        case 'FS'
            if options.debug > 0
                disp('Using the method of Finlayson and and Schaefer for recoverying the illuminant power spectrum.');
            end
            %Use the method of Finlayson and Schaefer
            L = recover_illuminant_finlayson_(I, options.patches);
        case {'GW', 'WP', 'SG', '1STOGE', '2NDOGE'}
            L = recover_illuminant_other_methods(I, options);
        otherwise
            error('Unknown method in options');        
    end
    
    %normalise the illuminant vector
    L = (L/norm(L))';
end