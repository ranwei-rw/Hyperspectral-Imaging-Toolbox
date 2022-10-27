%   Description
%       [L, MAPPING, KPMASK, SPMATRIX, CNRPOINTS] = recover_multi_illuminant_LM(I, method, theta, pitch_height, pitch_width, order)
%   recovers multiple illuminant from given images.
%
%   Input:  
%       I:            Image data cube
%       method:       a string specifying which method will be used.
%                     Available options are: 'grey_edge', 'grey_world',
%                     'maxrgb' and 'shades_of_grey'. Default value is
%                     'grey_edge'. 
%       theta:        values of theta. It could be a vector or a single
%                     value variable, defaults are [0.001, 0.03, 0.1, 3, 10].
%       pitch_height: height of pitches used in this function, default to 20.
%       pitch_width:  width of pitches used in this function, default to 20
%       order:        differential order default to 1
%   Ouptut:
%       L:            illuminants estimated, in row format
%       MAPPING:      matrix indicates each pitches is assigned to which
%                     illuminant
%       SPMATRIX:     Specularity mask
%       CNRPOINTS:    matrix holds coordinates of top left corners of
%                     pitches
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2015 All Rights Reserved.
% Author: Ran Wei. This code is based on the code of Lawrence Mutimbu of NICTA
%   version 1.0.0

function [L, MAPPING, SPMATRIX, CNRPOINTS] = recover_multi_illuminant_LM_(I, method, theta, pitch_height, pitch_width, order)
    
    % was like [lambda_p, q, energies] = recover_multi_illuminants_(dataDir, files, filesidx, PH, pW, theta_p, algorithm, diff_order)
    %   dataDir => I, files => none, filesinx=>none, PH => pitch_hight,
    %   pW => pitch_width, theta_p => theta, algorithm => method,
    %   diff_order => order
    
    if ~exist('order', 'var')
        order = 1;
    end

    if ~exist('pitch_width', 'var')
        pitch_width = 20;
    end
    
    if ~exist('pitch_height', 'var')
        pitch_height = 20;
    end
    
    if ~exist('theta', 'var')
        theta = [0.001, 0.03, 0.1, 3, 10];
    end
    
    if ~exist('method', 'var')
        method = 'grey_edge';
    else
        method = lower(method);
        if ~strcmp(method, 'grey_edge') && ~strcmp(method, 'grey_world') && ...
           ~strcmp(method, 'maxrgb') && ~strcmp(method, 'shades_of_grey')
            error('Unknown colour constancy algorithm.');
        end
    end
    
	%   Define Parameters for training
    %	lam_pvec = [0, 0.1, 0.3, 0.5, 0.8];
	lam = [0.1, 0.4, 0.7];
	qvec = [0, 0.5, 1];
    
    if isinteger(I)
        I = double(I);
    end

    maxi = max(I(:));
    if maxi == 0 
        mini = min(I(:));
        if mini == 0
            error('Could not handle all zero matrix');
        else
            %   negative value. flip them
            disp('Fliping pixel values');
            I = -I;
            maxi = -mini;
        end
    end

    if maxi ~= 1 
        %   normalise input 
        I = I/maxi;
    end

    %   replacing illuminantlabels by get_illuminant_labels_

    [LABELS, KPMASK, SPMATRIX, CNRPOINTS] = get_illuminant_labels_(I, pitch_height, pitch_width, method, order);


    if size(LABELS, 1) > 1
        [L, MAPPING] = graph_cut_optimise_(I, LABELS, KPMASK, SPMATRIX, CNRPOINTS, pitch_height, pitch_width, theta, lam, qvec);
    else
        L = LABELS;
        MAPPING = ones(size(SPMATRIX));
    end

    
end