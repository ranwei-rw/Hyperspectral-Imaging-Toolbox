%   Euclidean distance weight function
%  
%   Syntax
%
%       Z = dist(W, P); or If you want find the distance from each neuron
%                       to the other with, use
%       Z = dist(W);
% 
%   Description
% 
%       This function computes weights to an input to get weighted inputs.
% 
%   Inputs: 
%       W:	S-by-R weight matrix
%       P:  R-by-Q matrix of Q input (column) vectors
%   Output:
%       Z:  S-by-Q weight matrix
%   
%   Usage: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
%   Author: Ran Wei and Antonio Robles-Kelly
%   Version: 1.0.1
%   Last Update Date: 16 July 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = dist_(W, P)
    if nargin < 1 || nargin > 2
        error('Error found in input arguments. Please check.');
    end
    
    if nargin == 1
        P = W;
        W = W';
    end
    %   check whether dimensions meet    
    [s r1] = size(W);
    [r2 q] = size(P);
    if r1 ~= r2
        Error('Dimensions of input arguments do not match. Please check.');
    end
    
    %   Then we create an output matrix
    Z = zeros(s, q);
    
    %   compute values of it
    if s <= q
        %   when q is short, do calculation along 
        for i = 1:s
            %   compute Z(i, :) for each row of W
            V = P;
            for j = 1:r1
                V(j, :) = (P(j, :) - W(i, j)) .^ 2;
            end
            %   get sum for each 
            Z(i, :) = sqrt(sum(V, 1));
        end
    else
        %   when q is short, do calculation along 
        for i = 1:q
            %   compute Z(i, :) for each row of W
            V = W;
            for j = 1:r1
                V(:, j) = (W(:, j) - P(j, i)) .^ 2;
            end
            %   get sum for each 
            Z(:, i) = sqrt(sum(V, 2));
        end
    end

end