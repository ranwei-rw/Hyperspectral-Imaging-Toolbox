%   Function D = get_2d_gaussian_derivative_(I2D, sigma, xorder, yorder)
%   
%   is to get Gaussian derivative of a given 2D matrix using different
%   settings of sigma, and order
%   
%   Input:
%       I2D:    The input matrix
%       sigma:  value of sigma, default to 1; Sigma can not be 0.
%       xorder: order on x directrion, values between [0, 2], default to 0
%       yorder: order on y direction, values between [0, 2], default to 0
%
%   Output:
%       D:      The deriviatives of I2D
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2015 All Rights Reserved.
% Author: Ran Wei

function D = get_2d_gaussian_derivative_(I2D, sigma, xorder, yorder)

    if ~exist('yorder', 'var')
        yorder = 0;
    end
    
    if ~exist('xorder', 'var')
        xorder = 0;
    end
    
    if ~exist('sigma', 'var')
        sigma = 1;
    end
    
    if sigma == 0
        disp('sigma can not be 0. Reset it to 1');
        sigma = 1;
    end

    radius = floor(3*sigma+0.5);

    %I2D = fill_border(I2D, radius);
    [height, width] = size(I2D);
    
    TMP = zeros(height+radius*2, width+radius*2);
    %   central part
	TMP(1+radius:height+radius, 1+radius:width+radius) = I2D;
    %   four corner zones
    TMP(1:radius, 1:radius)                            = ones(radius, radius).*I2D(1, 1);
	TMP(height+radius+1:height+2*radius, 1:radius)     = ones(radius, radius).*I2D(height, 1);
	TMP(1:radius, width+radius+1:width+2*radius)       = ones(radius, radius).*I2D(1, width);
	TMP(height+radius+1:height+2*radius, width+radius+1:width+2*radius) = ones(radius, radius).*I2D(height, width);
	%   four edges
	TMP(1:radius, radius+1:width+radius)                        = ones(radius, 1)*I2D(1, :);
	TMP(height+radius+1:height+2*radius, radius+1:width+radius) = ones(radius, 1)*I2D(height, :);
	TMP(radius+1:height+radius, 1:radius)                       = I2D(:, 1)*ones(1, radius);
	TMP(radius+1:height+radius, width+radius+1:width+2*radius)  = I2D(:, width)*ones(1, radius);
    I2D = TMP;
    R = -radius:1:radius;

    G0 = 1/(sqrt(2 * pi) * sigma) * exp((-0.5*R.^2)/sigma^2 );

    switch xorder
        case 0
            G_x = G0/sum(G0);
        case 1
            G_x = -(R/sigma^2).*G0;
            G_x = G_x./sum(R.*G_x);
        otherwise
            G_x = (R.^2/sigma^4 - 1/sigma^2).*G0;
            G_x = G_x - sum(G_x)/size(R, 2);
            G_x = G_x/sum(0.5*R.^2.*G_x);
    end

    D = filter2(G_x, I2D);

    switch yorder
        case 0
            G_y = G0/sum(G0);
        case 1
            G_y = -(R/sigma^2).*G0;
            G_y = G_y./sum(R.*G_y);
        otherwise
            G_y = (R.^2/sigma^4-1/sigma^2).*G0;
            G_y = G_y-sum(G_y)/size(R, 2);
            G_y = G_y/sum(0.5*R.^2.*G_y);
    end
    
    D = filter2(G_y', D);

    D = D(radius+1:size(D, 1) - radius, radius+1:size(D, 2) - radius);

end