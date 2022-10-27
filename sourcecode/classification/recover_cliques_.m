
%   Usage: 
%       [CLIQUE, ORDER] = recover_cliques_(FACETS, POINTS);
%
%   Recovers the cliques and returns them on CLIQUEs and the order for each
%   clique
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CLIQUE, ORDER] = recover_cliques_(FACETS, POINTS)

    %   POINTS is a n by 2 matrix
    [n, ~] = size(POINTS);
    %   n should equal to max value in FACETS
    N = 30;
    CLIQUE = zeros(n, N);
    ORDER = zeros(n, 1);
    [n, m] = size(FACETS);


    for i = 1:n
        for j = 1:m
            %   get the value of current facets pixel as the centre point
            center = FACETS(i, j);
            for k = 1:m
                allocated = 0;
                for p = 1:ORDER(center)
                    if CLIQUE(center, p) == FACETS(i, k)
                        allocated = 1;
                        break;
                    end
                end
                if allocated == 0 && FACETS(i, k)~= center && ORDER(center)<N+1
                    ORDER(center) = ORDER(center)+1;
                    CLIQUE(center, ORDER(center)) = FACETS(i, k);
                end
            end
        end
    end
                
end