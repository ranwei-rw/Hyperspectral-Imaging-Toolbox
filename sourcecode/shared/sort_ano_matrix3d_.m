% Function 
%
%           RESULT = sort_ano_matrix3d-(SOURCE, TARGET, dim)
%
% this function is designed to sort matrix TARGET according to the result
% of sort(SOURCE, dim). This version works only on 3D matrices. Also,
% SOURCE and TARGET should have identical size. dim is the dimension along
% which the sorting operation goes.
%
%  For 2D matrices, please see sort_ano_matrix2d
%
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.6
% Last Update Date: 11 NOV 2013


function RESULT = sort_ano_matrix3d_(SOURCE, TARGET, dim, MODE)

    %   by default, sorting is carried on ascending
    if ~exist('MODE', 'var')
        MODE = 'ascend';
    end

    if ~exist('dim', 'var')
        dim = 1;
    end


    if ~isequal(size(SOURCE), size(TARGET))
        error('Dimensions of SOURCE and TARGET do not match');
    end

    [~, IND] = sort(SOURCE, dim, MODE);

    if (dim == 3)

        ix1 = repmat(repmat([1:size(SOURCE, 1)]', [1, size(SOURCE, 2)]), [1, 1, size(SOURCE, 3)]);
        
        ix2 = repmat(repmat([1:size(SOURCE, 2)], [size(SOURCE, 1)], 1), [1, 1, size(SOURCE, 3)]);
        
        %   get the index
        ix = sub2ind(size(SOURCE), ix1, ix2, IND);
        
    elseif (dim == 2)
        
        ix1 = repmat(repmat([1:size(SOURCE, 1)]', [1, size(SOURCE, 2)]), [1, 1, size(SOURCE, 3)]);
        
        ix2 = ones(size(SOURCE));
        for i = 2:size(SOURCE, 3)
            ix2(:, :, i) = i;
        end
        
        %   get the index
        ix = sub2ind(size(SOURCE), ix1, IND, ix2);
        
    else
        %   sort along the height
        ix2 = repmat(repmat([1:size(SOURCE, 2)], [size(SOURCE, 1)], 1), [1, 1, size(SOURCE, 3)]);
        ix3 = ones(size(SOURCE));
        for i = 2:size(SOURCE, 3)
            ix3(:, :, i) = i;
        end
        ix = sub2ind(size(SOURCE), IND, ix2, ix3);
    end

    RESULT = TARGET(ix);

end