
% Usage
% [PRED, PROB, C, model] = spectral_svmtest(TESTDATA, model)
%
% Description:
%     Computes the labels and PROBabilities using a spectral classifier
%     model.
%
% Inputs:
%     TESTDATA: Testing data. This is a structure with the following fields
%           W: The adjacency matrix used if the graph is not recomputed, i.e.
%               model.compute_graph == false. The matrix has to be indexed
%               such that the alpha variables in model.alpha is
%               sequentially consistent with the row and column indeces of
%               the training nodes of the graph.
%           X:  The data POINTS to be tested if the graph is recomputed, i.e.
%               model.compute_graph == 1.
%           Indx: If the adjacency matrix is used, i.e.
%               model.compute_graph == false, the Indx vector is used to separate
%               testing from training nodes. For training nodes, Indx(i) == 0
%               testing nodes are set to 1.
%     model: Spectral two-class classifier model structure
% 
% Output:
%
%     PRED: Predicted labels for the input data
%     PROB: Probability of the input data of belonging to either of the two
%           classes
%     C:    Centrality value. This is analogous to the discriminative function
%           value in SVMs.
%     model: the model used in the test procedure. It will contain both W
%           and Indx if model.compute_graph is false.
%
% Example:
%
%   training = annulus(200, 150, 1, 2, 5);
%   TESTDATA = annulus(20, 50, 1, 2, 5);
%   model = spectral_svmtrain(training);
%   model.compute_graph = true;
%   [PRED, PROB, C] = spectral_svmtest(TESTDATA, model);
% or
%   model.compute_graph = false;
%   POINTS = cat(1, TESTDATA.X, model.X);
%   TRI = delaunayn(POINTS);
%   testset.W = get_delaunay_weight(TRI, POINTS, model.h);
%   testset.Indx = cat(1, ones(70, 1), zeros(350, 1));
%   [pred, prob] = spectral_svmtest(TESTDATA, model);
%   disp_spectral_model(model);
% or 
%   traindata = gaussianclusters(4, 100, 25, 4);
%   options = struct('compute_graph', 1, 'h', 2, 'class', 4);
%   disp_spectral_model(model);
%   testdata = gaussianclusters(4, 100, 25, 4);
%   [ypred, prob, c] = spectral_svmtest(testdata, model);
%
% See also 
%
%   spectral_svmtrain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This computer code is subject to copyright: (C) National ICT Australia
%   Limited (NICTA) 2015 All Rights Reserved. 
% Author: Antonio Robles-Kelly and Ran Wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PRED, PROB, C, model] = spectral_svmtest(TESTDATA, model)

    [PRED, PROB, C, model] = spectral_svmtest_(TESTDATA, model);

end