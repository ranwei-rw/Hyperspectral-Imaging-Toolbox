%% Syntax
%       G = data2graph(Data);
%       G = data2graph(Data, options);
%       
%% Description:
%       Computes the weighted adjacency matrix for the 2D image matrix I.
% 
%% Input: 
% 
%   Data: the data matrix (stored as a 2D array with size n x m, where n is the number of data points and m is their
%         dimensionality). 
%   options: a structure containing the following fields
%            nn: Neigbhourhood to be used to build the adjacency matrix. For instance, if nn = 2, a neighbourhood of 2
%                pixels about. 
%            sigma: constant used to weight the data term on the weight matrix. Note the weights denote similarity, not
%                   distance. This is, a high value of the weight between two nodes implies these are "close" to each
%                   other. The default value is 10.
% 
%   Output: 
%
%   G: Structure containing the graph model. It contains the following fields: 
%      w, h: width and height of the input image. 
%      ii, jj, vv: Triplet for the graph edge term. These is such that the weighted adjacency matrix of the graph is
%                  given by W = sparse(G.ii, G.jj, G.vv, N, N),  where N = G.w*G.h;
%
%   Example:
%
%       load fisheriris
%       xdata = meas(51:end, 3:4);
%       G = data2graph(xdata);
%
%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.2
% Last Update Date: 28 July 2014

function G = data2graph_(I, options)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup the parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('options', 'var')
        options = [];   
    end;

    if ~isfield(options, 'nn')
        nn = 2; 
    else
        nn = options.nn;
    end;

    if ~isfield(options, 'sigma')
        sigma = 10;
    else
        sigma = options.sigma;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup the variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [N, w] = size(I);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the adjacency matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx  = 1;
    for i = 1:N
        [k, d] = knnsearch(I, I(i, :), 'K', nn+1);
        for j = 1:nn
            idxi      = i;
            idxj      = k(j+1);
            G.vv(idx) = exp(-d(j+1)^2*sigma);
            G.ii(idx) = idxi;
            G.jj(idx) = idxj;
            idx       = idx+1;
            G.vv(idx) = exp(-d(j+1)^2*sigma);
            G.ii(idx) = idxj;
            G.jj(idx) = idxi;
            idx = idx+1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Complete the structure and finish
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G.h = N;
    G.w = w;

end