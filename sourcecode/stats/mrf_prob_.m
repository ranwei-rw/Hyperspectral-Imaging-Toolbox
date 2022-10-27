% Syntax
%       [PLX_MAP, PLX_P] = mrf_prob(PLX_MAP);
%       
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.2
% Last Update Date: 1 Aug 2014

function [PLX_MAP, PLX_P] = mrf_prob_(PLX_MAP)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check the varaibles and compute the graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   output of I2Graph records the dimensions of PLX_MAP
    [height, width, bands] = size(PLX_MAP);
    
    if bands <= 1
        error('input argument should be a 3-dimensional array');
    end
    
    G      = I2Graph_(PLX_MAP);
    G.nCls = bands;
    tau    = 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the data term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G.mCls = diag(ones(G.nCls, 1)); 
    for i = 1:height
        for j = 1:width
            idxi = (j-1)*height+i;
            for k = 1:G.nCls
                G.e(idxi, k) = exp(-norm(reshape(PLX_MAP(i, j, :), 1, G.nCls)-G.mCls(k, :))^2*tau);
            end    
            t = sum(G.e(idxi, :));
            if t~=0 && ~isnan(t)
                G.e(idxi, :) = G.e(idxi, :)/sum(G.e(idxi, :));
            else
                %   this condition may never be met
                G.e(idxi, :) = 1/G.nCls;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the random walk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [LABELS, PROB] = randwalk_(G);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Transfer the labels and probabilities
    %to PLX_MAP and PLX_P
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    LABELS  = reshape(LABELS, height, width);
    PLX_P   = reshape(PROB, height, width, G.nCls);
    PLX_MAP = zeros(height, width, G.nCls);
    for i = 1 : G.nCls                
        PLX_MAP(:, :, i) = LABELS == i;
    end

%   end of function mrf_prob_
end