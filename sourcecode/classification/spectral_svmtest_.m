
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
%           X:  The data points to be tested if the graph is recomputed, i.e.
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
%   points = cat(1, TESTDATA.X, model.X);
%   TRI = delaunayn(points);
%   testset.W = get_delaunay_weight(TRI, points, model.h);
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

%   modified from Antonio's code SpectralTwoClassTest.m

function [PRED, PROB, C, model] = spectral_svmtest_(TESTDATA, model)

    if ~isfield(model, 'class')
        model.class = 2;
    end
    
    if strcmpi(model.type, 'oneclass')
        model.class = 1;
    elseif strcmpi(model.type, 'twoclass')
        model.class = 2;
    else
        model.class = 3;
    end
    
    if model.class > 2
        if ~isfield(model, 'Cindx')
            error('Field Cindx is missing from model');
        end
        [cluster_num, ~] = size(model.Cindx);
        [testsize, ~] = size(TESTDATA.X);
        C = ones(testsize, cluster_num);
        P = C;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute the classification results for all models
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:cluster_num
            for j = 1:cluster_num
                if i ~= j
                   [ypred_temp, prob_temp, c_temp] = ...
                       spectral_svmtest_(TESTDATA, model.model{abs(model.Cindx(i, j))});
                       Y = zeros(length(ypred_temp), 1);
                       if sign(model.Cindx(i, j)) == -1
                           Y(ypred_temp == 1) = 1;
                           P(1:testsize, i) = Y(1:testsize).*prob_temp(1:testsize).*P(1:testsize, i);
                           C(1:testsize, i) = Y(1:testsize).*c_temp(1:testsize).*C(1:testsize, i);
                       else
                           Y(ypred_temp == 2) = 1;
                           P(1:testsize, i) = Y(1:testsize).*prob_temp(1:testsize).*P(1:testsize, i);
                           C(1:testsize, i) = Y(1:testsize).*c_temp(1:testsize).*C(1:testsize, i);
                       end
                end
            end
        end
        [PROB, PRED] = max(P, [], 2);
        C = max(C, [], 2);
        PRED = PRED.*sign(max(P, [], 2));
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute the Delaunay graph and its adjancency matrix if
        %model.compute_graph == true.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if model.compute_graph

            [n, ~] = size(TESTDATA.X);
            [m, ~] = size(model.X);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Compute the probabilities and labels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PRED = zeros(n, 1);
            PROB = zeros(n, 1);
            C    = zeros(n, 1);
            Y    = (model.Y-1.5)*2;
            ALPHA = zeros(size(Y));
            ALPHA(model.SupportVectorIndices) = model.Alpha;
            X1 = model.X(:, 1);
            X2 = model.X(:, 2);
            for i = 1:n
                %   [~, indx] = ismember(model.X, TESTDATA.X(i, :));
                indx1 = X1 == TESTDATA.X(i, 1);
                indx2 = X2 == TESTDATA.X(i, 2);
                indx = indx1 & indx2;
                if sum(indx) == 0
                    %   when current test data is not presented in training
                    %   data
                    points = cat(1, TESTDATA.X(i, :), model.X);
                    TRI = delaunayn(points);
                    L = get_delaunay_weight_(TRI, points, model.h);
                    L = get_laplacian_matrix_(L);
                    C(i) = 1/sum(heaviside(abs(model.Alpha)))*(sum(Y.*abs(ALPHA).*reshape(L(2:m+1, 1), m, 1)));
                else
                    C(i) = -Y(indx == 1);
                end 

                PRED(i) = sign(C(i))*0.5+1.5;
                PROB(i) = 1-exp(-1/model.t*abs(C(i)));
            end
            %   normalise the possibilities
            PROB = PROB/max(PROB); 
        else

            % when the weight is not computed, input TESTDATA should contain
            % the element W
            if ~isfield(TESTDATA, 'W')
                error('adjacency matrix W is missing from input TESTDATA');
            end
            model.W = TESTDATA.W;

            if ~isfield(TESTDATA, 'Indx')
                error('Indx vector is missing from input TESTDATA');
            end

            model.Indx = TESTDATA.Indx;

            [n, ~] = size(TESTDATA.X);
            [m, ~] = size(model.X);
            L = laplacian_matrix(TESTDATA.W);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Compute the Probabilities and labels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            PRED  = zeros(n, 1);
            PROB  = zeros(n, 1);
            C     = zeros(n, 1);
            Y     = (model.Y-1.5)*2;
            L_IND = TESTDATA.Indx == 0;
            ALPHA = zeros(size(Y));
            ALPHA(model.SupportVectorIndices) = model.Alpha;
            for i = 1:n
                %   Include the labels into the testing for the one-class type
                t = 1/sum(heaviside(abs(model.Alpha)))*(sum(Y.*abs(ALPHA).*reshape(L(i, L_IND), [m, 1])));
                PRED(i) = sign(t)*0.5+1.5;
                C(i) = t;
                PROB(i) = 1-exp(-1/model.t*abs(t));
            end
            PROB = PROB/max(PROB); 
        end
    end

end