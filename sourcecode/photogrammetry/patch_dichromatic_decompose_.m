%
%   function [G, ...
%             K, ...
%             S, ...
%             dich_error, ...
%             smooth_error, ...
%             is_valid] = patch_dichromatic_decompose(I, ...
%                                                     L, ...
%                                                     alpha, ...
%                                                     debug)
%  This function is used to estimate g, k and S for a patch. It also
%   returns whether this patch is recognised as a valid dichromatic patch.
%   Given the light spectrum (i.e. the optimised g, k and S have positive
%   values). If so this patch will be used to optimise L in the subsequent
%   step in each iteration of the algorithm in LightEst4. Otherwise the
%   patch will be rejected from the computation. Note if I has more than 50
%   bands, then only 50 bands will be used.
%
%   Input:
%       I:     height x width x bands, the input patch.
%       L:     bands x 1, the known light spectrum.
%              NOTE: L must have all positive components.
%       alpha: a variable controlling the error threshold
%       debug: whether to show debug information, default value is 0
%
%   Output:
%
%       G: the shading map - shading factor
%       K: the specularity map. - specular coefficient
%       S: the reflectance of the patch. - spectral reflectance
%       is_valid: whether the select patch has the dichromatic
%       effect (is_valid == 1)
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.7
% Date: 4 June 2015
%   cleaned code and increased speed

% Version: 1.0.6
% Date: 16 Aug 2014

function [G, ...
          K, ...
          S, ...
          dich_error, ...
          smooth_error, ...
          is_valid] = patch_dichromatic_decompose_(I, ...
                                                   L, ...
                                                   alpha, ...
                                                   debug)
    
    if nargin == 3
        debug = 0;
    elseif nargin > 4 || nargin < 3
        error('Incorrect input argument number');
    end 
                                                           
	[height, width, bands] = size(I);

    %Go on
    if debug >= 2
        disp('Estimation shading, specularity and diffuse reflectance for input patches');
        disp('Program starts');
        s = sprintf('Input patch size: Height = %d, Width = %d, Bands = %d', height, width, bands);
        disp(s);
    end
    frame_size = height * width;
    [U, ~, ~] = svd(reshape(I, [frame_size, bands])');
    
    % Precompute the values of w1, w2, T1, T2, Z1, Z2
    A  = U(:, 1:2);
    Z1 = A(:, 1);
    Z2 = A(:, 2);
    
    TEMPW = A \ L;
    
    w1 = TEMPW(1);
    w2 = TEMPW(2);

    T   = A \ reshape(I, [frame_size, bands])';
    T1  = reshape(T(1, :),  [height, width]);
    T2  = reshape(T(2, :),  [height, width]);
    M   = w2 * T1 - w1 * T2;    
    DMX = M(2:height, :) - M(1:height-1, :);
    DMY = M(:, 2:width) - M(:, 1:width-1);
    N   = sum(sum(DMX .* DMX)) + sum(sum(DMY .* DMY));
    
    % Now compute the optimal value of the univariate variable.
    P = zeros(height, width, bands);
    Q = zeros(height, width, bands);

    for b = 1:bands
        P(:, :, b) = I(:, :, b) * w2 - (w2 * T1 - w1 * T2) * Z1(b) - L(b) * T2;
        Q(:, :, b) = I(:, :, b) * w1 - (w2 * T1 - w1 * T2) * Z2(b) + L(b) * T1;
    end        
    
    TMPPQ = P * w1/w2 - Q;     
    y     = -sum(sum(sum((P .* P) * w1/(w2^2) - (P .* Q)/w2)))/(sum(sum(sum(TMPPQ .* TMPPQ))) + alpha * N);
    a     = (1/y + w1)/w2;
    
      
    % Having the diffuse and LProj vectors (the basis of the dichromatic
    % plane). Find all the g and k.
    G = (w2 * T1 - w1 * T2) ./ (w2 * a - w1);
    K = (a * T2 - T1) ./ (w2 * a - w1);
    
    % Compute the diffuse_ reflectance (diffuse_ = <L, S> = 1 * a + Z2)
    %diffuse = A * [a; 1];
    S = (I - repmat(K, [1, 1, bands]).*...
        repmat(reshape(L, [1, 1, bands]), [height, width, 1]))./repmat(G, [1, 1, bands])./ repmat(reshape(L, [1, 1, bands]), [height, width, 1]);
    S(isnan(S)) = 0;
    
    if mean(S) < 0
        S = -S;
        G = -G;
    end

   
    % Normalize S, g and adjust their sign if they are negative (a HACK). 
    % The cost, the diffuse component and specular component are not
    % affected.
    
    % if the sign of the diffuse and G are opposite at
    % several pixels, then we reject this patch from
    % contribution to the computation of L.
    if min(G(:)) < -1e-4
        is_valid = 0;
    else
        is_valid = 1;        
    end

    % Clip the negative values of g and K.
    G = max(0, G);
    K = max(0, K);
    
    % normalize g 
    maxg = max(G(:));
    G    = G / maxg;
    S    = S * maxg;
    
    % Recompute the cost with the current values of g and k.
    I_DICH = zeros(height, width, bands); % value of I according to the dichromatic model.
    for b = 1:bands
        I_DICH(:, :, b) = (S(:, :, b) .* G + K) * L(b);
    end
    
    I_DIFF = (I - I_DICH).* (I - I_DICH);
    dich_error = sum(I_DIFF(:));

    % Compute the smoothness error from the (normalized) shading factors.
    [GX, GY] = gradient(G);
    smooth_error = alpha * sum(sum(GX.*GX + GY.*GY));
    
    if debug >= 2
        disp('Program ends');

        s = sprintf('Input patch size: dich_error = %d, smooth_error = %d\n', dich_error, smooth_error);
        disp(s);
        if is_valid == 1
            disp('Is a valid dichomatic patch');
        else
            disp('Is not a valid dichomatic patch');
        end
    end
    
end