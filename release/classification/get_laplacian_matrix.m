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

function L = get_laplacian_matrix(W)

    L = get_laplacian_matrix_(W);
    
end