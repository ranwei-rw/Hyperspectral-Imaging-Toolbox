% Usage
%   model = spectral_svmtrain(traindata, options);
%
% Description:
%
%   to train a svm model using given training data and options. 
%
%   Inputs:
%
%       traindata: Training data Set. This is a structure with the
%           following fields 
%           X:  Data points to be trained upon. Each data point should be
%               a row vector
%           Y:  Ground truth labels (Column vector) with the same length
%               as X does
%           W:  Weight matrix used when the graph is not to be computed from
%               the data points.
%
%       options: Determines whether the graph has to be computed or the
%           weight matrix is provided. This is a structure with the
%           following fields: 
%           compute_graph:   If this is 1, the graph is recomputed using a
%               Delaunay triangulation, if it is false, the input adjacency
%               matrix is used.
%           t:  Temperature used for the Gibbs field
%           h:  Bandwidth used for the computation of the weight matrix, 
%               default value is 10
%       class:  Type of train. class == 1 is for one class train while
%               class == 2 is for two class train. Default to 2: two class
%               training model. Class > 2 is for multi class training.
%             
%   Output:
%
%       model:  Structure containing the training model. This structure
%           contains the following fields
%
%           Alpha: Vector of Lagrange multipliers for the support
%               vectors. The sign is positive for support vectors
%               belonging to the first group and negative for
%               support vectors belonging to the second group.
%           Bias: Intercept of the hyperplane that separates
%               the two groups. Note: when 'autoscale' is false, this field
%               corresponds to the original data points in TRAINING. When
%               'autoscale' is true, this field corresponds to shifted and
%               scaled data points. 
%           SupportVectorIndices: A column vector indicating the indices of
%               support vectors.
%           TRI: N-D Delaunay triangulation
%           t:  Temperature for the Gibbs field
%           X:  Data points to be trained
%           Y:  Ground truth labels
%           type: Train type: 'TwoClass' or 'OneClass'
%           compute_graph: '1' when the graph was computed. '0', input was
%               used. 
%           h:  bandwidth used.
% 
% Example:
%
%   traindata = annulus_(200, 150, 1, 2, 5);
%   options = struct('compute_graph', 1, 'h', 2, 'class', 2);
%   model = spectral_svmtrain(traindata, options);
%   testset = annulus_(20, 50, 1, 2, 5);
%   [pred, prob] = spectral_svmtest(testset, model);
%  
%   model.compute_graph = false;
%   points = cat(1, testset.X, model.X);
%   TRI = delaunayn(points);
%   testset.W = get_delaunay_weight(TRI, points, model.h);
%   testset.Indx = cat(1, ones(70, 1), zeros(350, 1));
%   [pred, prob] = spectral_svmtest(testset, model);
%   disp_spectral_model(model);
%  
%   or
%   traindata = gaussianclusters(4, 100, 25, 4);
%   options = struct('compute_graph', 1, 'h', 2, 'class', 4);
%   model = spectral_svmtrain(traindata, options);
%
% See also 
%   spectral_svmtest, annulus, get_delaunay_weight, 
%   get_laplacian_matrix, delaunayn  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (c) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function model = spectral_svmtrain(traindata, options)
 
    switch nargin
        case 2
            model = spectral_svmtrain_(traindata, options);
        otherwise
            model = spectral_svmtrain_(traindata);
    end

end
