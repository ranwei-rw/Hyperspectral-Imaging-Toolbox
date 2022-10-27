%   This function is designed to train a SOM neuron network using given hyperspectral image cube data
%
% Syntax
%    [TRAINED_N2D, UPR] = som3d_training(N2D, I, options, loops, seed)
%
% Description
%
%   Train a 3D SOM neuron network using given hyperspectral image cube data
%    
% Input:
%
%   N2D: the neuron network matrix in 2D form which is of size (bands, length)
%
%   I2D:  2D training data converted from 3D Hyperspectral data cube - size (bands, height by width)
%
%   options: options includes height, width, sigma and learn_rate. height and width are first two dimensions of NEURON size.
%            sigma and learn_rate are options used to update neurons
%
%   loops: the number of total training steps. Default value is 10,000.
%   seed:  the seed used by randi function. Default is MATLAB default value
%
% Output:
%   
%   TRAINED_N2D: the neuron network trained
%   UPR    : map of update number on every neurons in 1D format. length of
%            UPR is the same as neuron number. values in UPR indicates
%            updates happened on that neuron.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.1
% Last Update Date: 23 Sep 2014

%   version changes:
%   1.0.4. Added seed function for random number generator

function [TRAINED_N2D, UPR]= som3d_training_(N2D, I2D, options, loopnum, seed)

    % Check number of inputs
    if nargin < 3
        error('Incorrect number of input arguments. Please check. Or see help using help som3d_training_');
    end
    
    if ~exist('loopnum', 'var')
        loopnum = 100000;
    end
    % check the input arguments to make sure that they are of legal size
    [bands, length] = size(N2D);
    
    [ibands, data_frame_size] = size(I2D);
    if bands ~= ibands
        error('the size of neuron vectors are not the same as hyperspectral bang number. Please check.');
    end
    
    %   define a map of training to record the number of winner neuron    
    UPR = zeros(length, 1);
    
    if exist('seed', 'var')
        rng(seed);
    end
    V = I2D(:, randi(data_frame_size, loopnum, 1));
    for i = 1:loopnum % define number of training loops 
        VP = V(:, i);
        N = N2D;
        %   Find winner neuron vector for each of the input from som
        %   network
        % for optimal reasons, using for loop and substract a scalar for each row or column along its short edge
        % here is to along its column direction
        for j = 1:bands
            N(j, :) = N(j, :) - VP(j);
        end

        % get its square
        [~, pos] = min(sum(abs(N), 1));
        
        %   update neuron matrix
        UPR(pos) = UPR(pos) + 1;
        %   update neuron network using this winner neuron
        update_rate = power(1.0/options.sigma, i/loopnum);
        options.sigma = options.sigma*update_rate;
        options.learn_rate = options.learn_rate*update_rate;
        if i == 1
            [N2D, COR] = update_winner_neighbourhood_(N2D, VP, pos, options);
        else
            N2D = update_winner_neighbourhood_(N2D, VP, pos, options, COR);
        end
    end
    
    TRAINED_N2D = N2D;
end
    