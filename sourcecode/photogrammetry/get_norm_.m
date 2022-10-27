%  Calculate norm (distance) between two input matrices
%
% Syntax
%   NORM = get_norm(M, N, form);
%
% Description
%
%   Calculate the norm between two input matrices M and N in the form of L1 norm or L2 norm. M and N can be vectors, 
%   2D and 3D matrix as long as they are of same dimensions.
%    
% Input:
%   M and N: vectors, 2D or 3D matrix of source image
%   form:    'L1' or 'L2' or 'L1norm' or 'L2norm' or 'angle', for L1, L2
%            norm and angle distance accordingly, case insensitive. When
%            form is angle, the values of output NORM is cos(theta)
%
% Output:
%   NORM:    norm value(s) between two input matrices. If M is a scalar or a vector, NORM will be a scalar. If M 
%            is a 2D matrix, NORM will be a row vector. If M is a 3D matrix of h by w by b, NORM will be a matrix of 
%            h by w.
%   
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014-2015 All Rights Reserved.

% Author: Ran Wei
% Version: 1.0.1
% Last Update Date: 5 Nov 2015

function NORM = get_norm_(M, N, form)

    %   check the number of inputs;
    if nargin < 2
        error('Not enough input arguments');
    end
    
    %   check condition
    if ~exist('form', 'var') 
        form = 'L2';
    end
    
    form = lower(form);
    
    if strcmp(form, 'l1') || strcmpi(form, 'l1norm')
        f = 1;
    elseif strcmpi(form, 'l2') || strcmpi(form, 'l2norm')
        f = 2;
    elseif strcmpi(form, 'angle')
        f = 3;
    else
        error('Unknown option for form');
    end

    %   check whether inputs M and N are of same size
    samesize = 1;
    if ~isequal(size(M), size(N))
        samesize = 0;
    end
    
    %   now inputs are verified, proceed to calculate 
    %   special case: a value
    if samesize
        %   when inputs have same size
        [h, w, b] = size(M);
        if h*w*b == 1
            %   one value matrix
            NORM = abs(M - N);
            return;
        end

        %   when M and N are vectors
        if isrow(M) || iscolumn(M)
            DIFF = M - N;
            NORM = DIFF;
            if f == 1
                NORM = sum(abs(DIFF));            
            else
                DIFF = DIFF .^ 2;
                NORM = sqrt(sum(DIFF));
            end
        end

        % 2D and 3D cases, support they are column vectors of interest
        if b == 1
            % 2D cases
            DIFF = M - N;
            NORM = DIFF;
            if f == 1
                NORM = sum(abs(DIFF), 1);  
            else
                DIFF = DIFF .^ 2;
                NORM = sqrt(sum(DIFF, 1));
            end
        else
            % 3D cases, suppose that the norm is taken along 3rd dimension.

            if f == 1
                NORM = sum(abs(M - N), 3);  
            elseif f == 2
                NORM = sqrt(sum((M - N) .* (M - N), 3));
            else
                P = sqrt(sum((M .* M), 3)) .* sqrt(sum((N .* N), 3));
                P(P == 0) = 1e-10;
                NORM = sum(M.*N, 3)/P;
            end
        end
    else
        %   when inputs have different sizes
        [mh, mw, mb] = size(M);
        [nh, nw, nb] = size(N);
        
        if mb == nb && nb ~= 1
            %   3D matrix
            if mh ~= nh || mw ~= nw
                %   this case happens when two 3D matrices are passed in, they
                %   are of different frame but have same depth (Bands). It
                %   means people want to know the norm between them for each of
                %   the pixels. It will result in a 2D matrix NORM while the
                %   height of NORM is mh*mw and width of NORM is nh*nw
                %   1 reshape them into a new 3D
                M = reshape(M, [mh*mw, 1, mb]);
                N = reshape(N, [1, nh*nw, nb]);
                
                M3D = repmat(M, [1, nh*nw, 1]);
                N3D = repmat(N, [mh*mw, 1, 1]);
                
                %   now these two matrix has the same dimension. compute
                %   them
                if f == 1
                    NORM = sum(abs(M3D - N3D), 3);
                elseif f == 2
                    NORM = sqrt(sum((M3D - N3D) .* (M3D - N3D), 3));
                else
                    P = sqrt(sum((M3D .* M3D), 3)) .* sqrt(sum((N3D .* N3D), 3));
                    P(P == 0) = 1e-10;
                    NORM = sum(M3D.*N3D, 3)/P;
                end
            else
                %   two identical 3D matrix which will already handled by
                %   above code.
            end
            
        elseif mb == 1 && nb == 1;
            %   two 2D matrix
            if mw == nw && mh == nh
                %   identical 2D matrx
                if f == 1
                    NORM = abs(M - N);
                elseif f == 2
                    NORM = sqrt((M-N).*(M-N));
                else
                    P = M.*N;
                    P(P == 0) = 1e-10;
                    NORM = M.*N./P;
                end
            elseif mw ~= nw && mh == nh
                %   two matrix with same height, actually it compares
                %   distance between columns
                M = reshape(M, [mw, 1, mh]);
                N = reshape(N, [1, nw, nh]);
                M3D = repmat(M, [1, nw, 1]);
                N3D = repmat(N, [mw, 1, 1]);
                
                %   now these two matrix has the same dimension. compute
                %   them
                if f == 1
                    NORM = sum(abs(M3D - N3D), 3);
                elseif f == 2
                    NORM = sqrt(sum((M3D - N3D) .* (M3D - N3D), 3));
                else
                    P = sqrt(sum((M3D .* M3D), 3)) .* sqrt(sum((N3D .* N3D), 3));
                    P(P == 0) = 1e-10;
                    NORM = sum(M3D.*N3D, 3)/P;
                end
            elseif mw == nw && mh ~= nh
                M = reshape(M, [mh, 1, mw]);
                N = reshape(N, [1, nh, nw]);
                M3D = repmat(M, [1, nh, 1]);
                N3D = repmat(N, [mh, 1, 1]);
                
                %   now these two matrix has the same dimension. compute
                %   them
                if f == 1
                    NORM = sum(abs(M3D - N3D), 3);
                elseif f == 2
                    NORM = sqrt(sum((M3D - N3D) .* (M3D - N3D), 3));
                else
                    P = sqrt(sum((M3D .* M3D), 3)) .* sqrt(sum((N3D .* N3D), 3));
                    P(P == 0) = 1e-10;
                    NORM = sum(M3D.*N3D, 3)/P;
                end
            else
                error('Please check the size of matrix');
            end
            
        elseif mb ~= nb
            %   3D matrix but with differnt depth.
            error('Please check the size of matrix');
        end
        
    end
end