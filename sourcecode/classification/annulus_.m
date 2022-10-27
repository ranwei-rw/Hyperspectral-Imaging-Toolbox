%   Usage:
%       data = annulus(cen_num, peri_num, annu_num, step, spread);
%
%   Generates data for svm training and testing which comprises a set of
%   annulea in two dimensions.
%
%	Inputs: 
%       cen_num:  Number of data points in the centre of each annulus
%       peri_num: Number of points in the aureola
%       annu_num: Number of annulea, default to 1
%       step:     Displacement in the subsequently generates annulea,
%                 default to 2;
%       spread:   Spread of the annulus, default to 5;
%
%	Output:
%       data: a data struct with the following members
%                .X: [cen_num x peri_num x annu_num by 2] Training vectors.
%                .Y: [cen_num x peri_num x annu_num by 1] Labels (1 or 2).
%
%   Example:
%       trn = annulus(200, 150, 1, 2, 5);
%       where trn.X contrains 350 training data units while trn.Y contains
%       350 labels.
%
%   See also: milkyway, gaussianclouds 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = annulus_(cen_num, peri_num, annu_num, step, separation)

    if ~exist('separation', 'var')
        separation = 5;
    end
    
    if ~exist('step', 'var')
        step = 2;
    end
    
    if ~exist('annu_num', 'var')
        annu_num = 1;
    end

    for i = 1:annu_num
        center    = [step*(i-1)+separation*randn(cen_num, 1), step*(i-1)+separation*randn(cen_num, 1)];
        rho       = separation*step*i+annu_num*abs(randn(peri_num, 1));
        theta     = pi*randn(peri_num, 1);
        periphery = [step*(i-1)+separation*rho.*sin(theta), step*(i-1)+separation*rho.*cos(theta)];
        if i == 1
            data.X = cat(1, center, periphery);
            data.Y = cat(1, i*ones(cen_num, 1), 2*ones(peri_num, 1));
        else
            data.X = cat(1, data.X, center, periphery);
            data.Y = cat(1, data.Y, (2*(i-1)+1)*ones(cen_num, 1), (2*i)*ones(peri_num, 1));
        end
    end

end