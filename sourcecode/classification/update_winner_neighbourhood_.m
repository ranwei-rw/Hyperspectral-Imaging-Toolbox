%   This function is used to update the neighbourhood in a 3D neuron network around a winner neuron 
% Syntax
%    UPD_N = update_winner_neighbourhood_(NEURONS, INPUT, pos, options, COR, debug)
%
% Description
%
%   Update the area around a winner neuron in a som 3D neuron matrix.
%    
% Input:
%
%   NEURONS: the neuron network in 2D form which is of size (bands, length)
%
%   INPUT:   the vector of input which is of size (bands, 1)
%
%   pos:     the position of the winner neuron in matrix, around which
%            other neurons will be updated 
%
%   options: options includes height, width, sigma, learn_rate, radius and threshold. 
%            height, width, depth: dimensions of NEURON
%            sigma and learn_rate: options used to update neurons. by
%                                  default, sigma == 4, learn_rate == 0.75
%            radius:               the radius of updating neurons around the winner
%            threshold:            the threshold of weight below which the
%                                  particular neuron will not be updated.   
%            debug:                Level of debuging information, default
%                                  to 1, max 3
%   COR:     Coordinates of neurons in 2D format. If COR is not given, it
%            will be generated automatically.
%
%   debug:   level of debugging information to be shown.
%
% Output:
%   
%   UPD_N:   Updated neuron network
%   COR:     Coordinates of neurons in 2D format. This matrix can be reused
%            later to save time.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 - 2015 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.3
% Last Update Date: 16 June 2015
% added an option to reuse pre-calculated COR matrix to speed up.


% Version: 1.0.2
% Last Update Date: 22 Sep 2014

function [UPD_N, COR] = update_winner_neighbourhood_(NEURONS, INPUT, pos, options, COR, debug)

    %   make sure the options is of proper structure
    
    if ~isfield(options, 'width')
        error('Argument options does not have width information. Please check');
    end
    
    if ~isfield(options, 'height')
        error('Argument options does not have height information. Please check');
    end
    
    if ~isfield(options, 'depth')
        error('Argument options does not have depth information. Please check');
    end
    
    if ~exist('debug', 'var')
        debug = 1;
    end
    
    if ~exist('COR', 'var')
        nocor = 1;
    else
        nocor = 0;
    end
    
    if ~isfield(options, 'sigma')
        if debug > 1
            disp('Argument options does not have sigma information. Using default value: 4.0');
        end
        options.sigma = 4;
    end

    if ~isfield(options, 'learn_rate')
        if debug > 1
            disp('Argument options does not have learn_rate information. Using default value: 0.75');
        end
        options.learn_rate = 0.75;
    end

    if ~isfield(options, 'radius')
        if debug > 1
            disp('Argument options does not have radius information. Using default value: 4');
        end
        options.radius = 4;
    end

    if ~isfield(options, 'threshold')
        if debug > 1
            disp('Argument options does not have threshold information. Using default value: 0.01');
        end
        options.threshold = 0.01;
    end

    height     = options.height;
    width      = options.width;
    depth      = options.depth;
    sigma      = options.sigma^2;
    learn_rate = options.learn_rate;
    neuron_num = height*width*depth;
    
    if nocor
        %   here we need to calculate the distance between this neuron and the
        %   winner neuron based on their 3D Euclidean distance
        COR = zeros(3, neuron_num);

        P = 1:neuron_num;
        Z = ceil(P/(height*width));
        R = mod(P, height*width);
        X = ceil(R/height);
        X(X==0) = width;
        Y = mod(R, height);
        Y(Y==0) = height;
        COR(1, :) = Y;
        COR(2, :) = X;
        COR(3, :) = Z;
    end
    
    %  firstly, we need to know which layer the winner neuron is on (its z coordinate)
    CO = COR(:, pos);
    
    %   need to work out the range of updating. using the radius infomation
    yp = min((CO(1) + options.radius), height);
    yn = max(1, (CO(1) - options.radius));
    xp = min((CO(2) + options.radius), width);
    xn = max(1, (CO(2) - options.radius));
    zp = min((CO(3) + options.radius), depth);
    zn = max(1, (CO(3) - options.radius));
    %   now we have the neurons that to be updated, we need then to get
    %   them into a new sub matrix and calculate the new coordinates of CO
    %   in this new matrix
    
    TBUPD = zeros(height, width, depth); %TBUPD == to be updated
    TBUPD(yn:yp, xn:xp, zn:zp) = 1;
    INDEX = find(TBUPD == 1); %now INDEX contains the index of neurons to be updated
    
    CORD = COR;
    
    for i = 1:3
        CORD(i, :) = COR(i, :) - CO(i);
    end
    
    D = sum(CORD.*CORD);
    WEIGHTS = learn_rate * exp(-D ./ (2.0*sigma));
    WEIGHTS(WEIGHTS < options.threshold) = 0;
    WEIGHTS = WEIGHTS(INDEX);
    
    %   get the neurons out into a new sub matrix
    N = NEURONS(:, INDEX);
    bands = size(N, 1);
    
    for i = 1:bands
        N(i, :) = N(i, :) + (INPUT(i) - N(i, :)).*WEIGHTS;
    end
    
    %   put updated neurons back into original positions
    NEURONS(:, INDEX) = N;

    UPD_N = NEURONS;
    clear NEURONS;
    clear N;

    %end of function
end