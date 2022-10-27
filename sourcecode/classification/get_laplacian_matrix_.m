%   Usage
%       L = get_laplacian_matrix(W);
%
%   Description:
%
%     Computes the Laplacian matrix for the weight matrix W.
%
%   Input:
%       
%       W:  weight matrix, a square matrix.
%
%   Output:
%       
%       L: Laplacian matrix computed from W.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = get_laplacian_matrix_(W)

    %   Recover the cliques from the x and y coords
    [n, ~] = size(W);

    %   normalise input weight matrix
    minw = min(W(:));
    maxw = max(W(:));
    
    if maxw - minw ~= 0
        W = (W - minw)/ maxw - minw;
    else
        if minw ~= 0
            W = W/minw;
        end
    end

    D_INV = zeros(n, n);
    %   Compute the Laplacian
    D = sum(W, 2);
    D = diag(D);
    IND = D~=0;
    D_INV(IND) = 1./sqrt(D(IND));
    
    L = D_INV*(D - W)*D_INV;


end