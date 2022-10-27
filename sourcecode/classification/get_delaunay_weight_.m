% Usage
%   W = get_delaunay_weight_(TRI, POINTS, h);
%
% Description:
%
%   Computes the weight matrix for the Delaunay triangulation with
%   simplexes TRI and node coordinates POINTS. It uses the constant h for
%   the graph weights computed using an exponential decay over the distance
%   for adjacent nodes.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = get_delaunay_weight_(TRI, POINTS, h)

    %   Recover the cliques from the x and y coords
    [n, m] = size(POINTS);
    MAXP = max(POINTS);
    MINP = min(POINTS);
    %   Don't take into consideration distances below threshold
    threshold = 0.00;

    %   Normalise POINTS
    for i = 1:m
        if MAXP(i) - MINP(i) ~= 0
            POINTS(:, i) = (POINTS(:, i) - MINP(i))/(MAXP(i) - MINP(i)); %4 is the scaling factor
        end
    end

    %Recover the cliques
    [clique, clique_order] = recover_cliques_(TRI, POINTS);


    %   Compute the matrix W using the cliques
    W = zeros(n, n);
    
    for i = 1:n
        for j = 1:clique_order(i)
            W(i, clique(i, j)) = norm(POINTS(i, :) - POINTS(clique(i, j), :));
        end
    end

    %   Compute the final weight matrix W
    for i = 1:n
        for j = i:n
            if W(j, i) > threshold   %Avoid elements of W which are out of the value for 3 sigma
                W(j, i) = exp(-h*W(i, j));%Apple
            else
                W(j, i) = 0;
            end
            W(i, j) = W(j, i);
        end
    end

end