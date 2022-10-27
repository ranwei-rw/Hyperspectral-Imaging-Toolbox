%   Syntax
%       G = I2Graph(I);
%       G = I2Graph(I, options);
%       
%  Description:
%       Computes the weighted adjacency matrix for the image cube I.
% 
%  Input: 
%       I:       the radiance image (stored as a 3D array with size height x width x bands).
%       options: a structure containing the following fields
%                nn:    Neigbhourhood to be used to build the adjacency matrix. For instance, if nn = 2, a neighbourhood
%                       of 2 pixels about. 
%                sigma: constant used to weight the color/spectra term on the weight matrix. Note the weights denote
%                       similarity, not distance. This is, a high value of the weight between two nodes implies these
%                       are "close" to each other. The default value is 1. 
%                gamma: Constant used to weight the spatial term. The default value is 1.
% 
%  Output: 
%       G: Structure containing the graph model. It contains the following fields
%           width, height: width and height of the input image.
%           ii, jj, vv:    Triplet for the graph edge term. These are such that the weighted adjacency matrix of the
%                          graph is given by W = sparse(G.ii, G.jj, G.vv, N, N), where N = G.width*G.height;
%
%  Example:
%
%       imdata = imread('./shared/samples/ngc6543a.jpg');
%       G = I2graph(imdata);
%
% This computer code is subject to copyright: (bands) National ICT Australia Limited (NICTA) 2014-2015 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.3
% Last Update Date: 18 June 2015
% Speed up

% Version: 1.0.2
% Last Update Date: 29 July 2014

function G = I2Graph_(I, options)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup the parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('options', 'var')
        nn  = 2;
        sigma = 0.1;
        gamma = 0.2;
    else
        if ~isfield(options, 'nn')
            nn = 2;
        else
            nn = options.nn;
        end

        if ~isfield(options, 'sigma')
            sigma = 0.1;
        else
            sigma = options.sigma;
        end

        if ~isfield(options, 'gamma')
            gamma = 0.2;
        else
            gamma = options.gamma;
        end      
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	normalise input image and get its dimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = double(I);
    I = I/max(I(:));
    [height, width, bands] = size(I);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Compute the adjacency matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx = 1;
    Q = zeros(2*nn+1, 1);
    for i = 1:nn
        Q(i)    = i;
        Q(i+nn) = -i;
    end
    Q(2*nn+1) = 0;
    %   now Q has 2*nn+1 elements
    for i = 1:height
        for j = 1:width
            idxi = height*(j-1)+i;
            for k = 1:2*nn+1
                for m = 1:2*nn+1
                    if Q(k)+Q(m) == 0
                        continue;
                    end
                    
                    if (i-Q(k))>0 && (j-Q(m))>0 && (i-Q(k))<height+1 && (j-Q(m))<width+1
                        idxj = height*(j-Q(m)-1)+i-Q(k);
                        G.vv(idx) = exp(-norm(reshape(I(i, j, :) - I(i-Q(k), j-Q(m), :), 1, bands))^2*sigma)*exp(-(Q(m)*Q(m)+Q(k)*Q(k))*gamma);
                        G.ii(idx) = idxi;
                        G.jj(idx) = idxj;
                        idx = idx+1;
                    end
                end
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Complete the structure and finish
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G.height = height;
    G.width = width;
    
end