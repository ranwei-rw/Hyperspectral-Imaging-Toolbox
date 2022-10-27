%   Usage
%       disp_spectral_model(model);
%
%   Description:
%       Plots the probability of classes given the model delivered by
%       spectral_svmtrain
%
%   Inputs:
%       model: One or two class classifier model structure
% 
%   Example:
% 
%       trn = annulus(200, 150, 1, 2, 5);
%       options.class = 1;
%       model = spectral_svmtrain(trn, options);
%       disp_spectral_model(model);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function demo_2d_spectral_model(model)

    demo_2d_spectral_model_(model);

end