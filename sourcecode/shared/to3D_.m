%   this function was initially designed to convert a 1D illuminant vector
%   L into a 3D form so that other functions may use it easier.
%   Description:
%       function L3D = to3D(L, size3D)
%   
%   Input:
%       L: the 1D vector or 2D matrix
%       size3D: the size of 3D matrix that L is going to be converted to
%
%   Output:
%       L3D: 3D version of L with the 3D size specified in size. L's
%            elements will be repeated on band direction (3rd one in size3D).
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.0
% Last Update Date: 14 Oct 2014

function L3D = to3D_(L, size3D)

    if length(size3D) ~= 3
        error('Error in size information');
    end
    height = size3D(1);
    width  = size3D(2);
    bands  = size3D(3);
    
    [lh, lw, ln] = size(L);
    
    diml = 0;
    if ln ~= 1
        %   3d
        diml = 3;
        if (lh ~= 1 && lh ~= height) || (lw ~= 1 && lw ~= width)
            %   L has to be a 3D matrix with either 1 value per band or a
            %   full dimension pixel-wise matrix
            error('Dimensions of L is not correct');
        end
    elseif lh ~= 1 && lw ~= 1
        %   2D matrix, same for all bands
        if lh ~= height || lw ~= width
            error('Dimensions of L is not correct.');
        end
        diml = 2;
    else
        %   vecotr
        if lh == 1 && lw == 1
            error('Dimensions of L is not correct. A vector or matrix is required');
        end
        diml = 1;
    end

    L3D = zeros(size3D);
    switch diml 
        case 3
            L3D = L;
        case 2
            for b = 1:bands  
                L3D = repmat(L, 1, 1, bands);
            end
        case 1
            L = reshape(L, [1, 1, bands]);
            L3D = repmat(L, height, width);
        otherwise
            error('Dimensions of L is not correct');
    end
    
%   end of function

end