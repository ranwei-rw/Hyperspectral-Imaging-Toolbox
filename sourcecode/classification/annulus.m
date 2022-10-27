%   Usage:
%       data = annulus(cen_num, peri_num, annu_num, step, spread);
%
%   Generates data for svm training and testing which comprises a set of
%   annulea in two dimensions.
%
%	Inputs: 
%       cen_num:  Number of data points in the centre of each annulus
%       peri_num: Number of points in the aureola
%       annu_num: Number of annulea
%       step:     Displacement in the subsequently generates annulea
%       spread:   Spread of the annulus
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

function data = annulus(cen_num, peri_num, annu_num, step, separation)

    switch nargin
        case 5
            data = annulus_(cen_num, peri_num, annu_num, step, separation);
        case 4
            data = annulus_(cen_num, peri_num, annu_num, step);
        case 3
            data = annulus_(cen_num, peri_num, annu_num);
        case 2
            data = annulus_(cen_num, peri_num);
        otherwise
            error('Incorrect number of inputs');
    end
end