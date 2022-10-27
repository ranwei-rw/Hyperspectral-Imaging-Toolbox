%   Usage
%       bias = bias = get_linearsvm_bias_(LABELS, POINTS, alpha, computegraph, h)
%
%   Description:
%       Computes the bias for linear SVM test results
%
%   Inputs:
%       LABELS: Labels for the two training clases. Labels are 1 and 2, 
%               respectively.
%       POINTS: Input data points
%       alpha: Vector of alpha weights
%       h:     bandwidth used. Default to 2.
% 
%   Output
%       bias: Bias
%
% See also 
%   spectral_svmtrain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (C) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bias = get_linearsvm_bias_(LABELS, POINTS, alpha, h)

    if ~exist('h', 'var')
        h = 2;
    end

    TRI = delaunayn(POINTS);
    W = get_delaunay_weight_(TRI, POINTS, h);
    L = get_laplacian_matrix_(W);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the bias
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y = 2*(LABELS-1.5);
    bias = 1/sum(heaviside(alpha))*sum((Y.*alpha)*L-sum(Y));
    
end
