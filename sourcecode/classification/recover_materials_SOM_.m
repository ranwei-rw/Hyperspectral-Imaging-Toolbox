%   Recover materials from given hyperspectral image cube by using Self Organising Map (SOM) of neural network.
%   [WINNER, MSSOM, SOM, COUNT] = recover_materials_SOM_(I, som3d_height, som3d_width, som3d_depth, SOM, loopnum, seed)
%   
%   Input
%
%       I           : A 3D hyperspectral image cube
%       som3d_height: Height of som neuron network matrix, Default to 10
%       som3d_width : Width of som neuron network matrix, Default to 10
%       som3d_depth : depth of som neuron network matrix, Default to 10
%       SOM         : Optional pre-set neuron network matrix, by default
%                     it is a 3D random matrix of size mentioned above
%       loopnum     : is the number of loops to be conducted in training,
%                     default to 100,000
%       seed  :  the seed used by randi function. Default is MATLAB default value
%
%   Output
%
%       WINNER: material index categorised. A 2D matrix which has the same
%               frame size as the original image.
%       MSSOM : the SOM neurons after meanshift.
%       COUNT : pixel number for each neurons (material recognised)
%       SOM   : The neuron network used in this computation. It's in 2D
%               form (bands by neuron number). This SOM is not clusterred by
%               meanshift.
%       seed  :  the seed used by randi function. Default is MATLAB default value
%
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.4
% Last Update Date: 1 Dec 2014
%

%   version changes:
%   1.0.4. Added seed function for random number generator

function [WINNER, MSSOM, SOM, COUNT] = recover_materials_SOM_(I, som3d_height, som3d_width, som3d_depth, SOM, loopnum, seed)

    %   I is the input hyperspectral image cube, it should be of the size
    %   height by width by bands. To make the process faster, in this function we need to transform
    %   it into 2D form which is height*width by bands
    if isstruct(I)
        if isfield(I, 'I')
            I = I.I;
        else
            error('Please check input data');
        end
    end
    [height, width, bands] = size(I);
    
    if bands == 1
        error('I should be a 3D matrix. Please check it');
    end
    
    given_som = 0;
    if exist('SOM', 'var') && ~isscalar(SOM)
        %   a som matrix is given, use it instead of creating a new one
        given_som = 1;
    end
    
    if ~exist('som3d_height', 'var')
        som3d_height = 10;
    end
    
    if ~exist('som3d_width', 'var')
        som3d_width = 10;
    end
    
    if ~exist('som3d_depth', 'var')
        som3d_depth = 10;
    end
    
    if ~exist('loopnum', 'var')
        loopnum = 100000;
    end
    
    %options includes height, width, sigma and learn_rate. height and width are first two dimensions of NEURON size.
    %            sigma and learn_rate are options used to update neurons
    %   then we need to train the neuron network
    options.height     = som3d_height;
    options.width      = som3d_width;
    options.depth      = som3d_depth;
    options.radius     = 8;
    options.threshold  = 0.01;
    options.sigma      = 4;
    options.learn_rate = 0.75;
    length = height*width;
    
    I2D = reshape(I, [length, bands])';
    
    if ~given_som
        if exist('seed', 'var')
            rng(seed);
        end
        %   then we need to create the som2d weight matrix
        SOM = rand([bands, som3d_height * som3d_width * som3d_depth]);
        %   if SOM is not given, we need to train it. Otherwise, we'll jump
        %   over training section
        SOM = som3d_training_(SOM, I2D, options, loopnum);
    end

    %   Training finished. Once we get the SOM network trained. We need to
    %   feed this SOM with actual image pixels for testing/classifying. 
    %  Create a list for all image vectors for their winner neuron location in SOM2D and cout for them
    
    
    WINNER = zeros(length, 1);
    COUNT  = zeros(som3d_height * som3d_width * som3d_depth, 1);
    
    for i = 1:length
        V = I2D(:, i);
        N = SOM;
        %   Find winner neuron vector for each of the input from som
        %   network
        % for optimal reasons, using for loop and substract a scalar for each row or column along its short edge
        % here is to along its column direction
        for j = 1:bands
            N(j, :) = N(j, :) - V(j);
        end

        [~, pos] = min(sum(abs(N), 1));

        COUNT(pos) = COUNT(pos) + 1;
        WINNER(i) = pos;
        clear N;
    end
    
    %   Here: WINNER(i) records the index of the neuron for pixel i.
    %   Count(pos) records how many pixels that neurons pos has. This count
    %   is also used as weights and spsizes (line 84-85 in
    %   meanshift_som.cpp). The difference between weights and spsizes is
    %   that weights is then normalised against to its mean later (line
    %   91-93 in meanshift_som.cpp). vector spsizes is later copied into
    %   meanshift data structure where it keeps the same name.
    
    
    
    %   USE meanshift to cluster all neurons into different clusters. In
    %   meanshift algorithm, SOM neurons are called modes.
    [MSSOM, MAPPING] = meanshift_(SOM, COUNT);
    
    %   now we have the final mapping m
    for i = 1:length
        WINNER(i) = MAPPING(WINNER(i));
    end
    
    %reshape WINNER matrix
    WINNER = reshape(WINNER, [height, width]);

    %end of function
end