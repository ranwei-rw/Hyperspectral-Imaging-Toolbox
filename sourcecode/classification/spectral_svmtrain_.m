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


function model = spectral_svmtrain_(traindata, options)
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   varify input variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if ~exist('options', 'var') || (exist('options', 'var') && ~isfield(options, 'W'))
        options.compute_graph = 1;
    end
    
    if ~isfield(options, 't') || options.t < 0
        options.t = 1;
    end
    
    if ~isfield(options, 'h') || options.h < 0
        options.h = 10;
    end
    
    if ~isfield(options, 'class') || options.h < 0
        options.class = 2;
    end
    
    if ~isfield(traindata, 'Y')
        error('Missing ground truth data Y in training data.');
    end

    if ~isfield(traindata, 'X')
        error('Missing training data X in input.');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the Delaunay graph and its adjancency matrix if the graph is
    %supposed to be computed. Use the provided weight matrix otherwise.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.compute_graph
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute the Delaunay triangulation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TRI = delaunayn(traindata.X);
        model.TRI = TRI;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute the adjacency matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W = get_delaunay_weight_(TRI, traindata.X, options.h);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute the lapalacian matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L = get_laplacian_matrix_(W);
        model.W = W;
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute the lapalacian matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L = laplacian_matrix(traindata.W);
        model.W = traindata.W;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.class > 2
        model.type = 'MultiClass';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compute the one-versus all models
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clusters = max(traindata.Y);
        t = 1;
        for i = 1:clusters
            for j = i+1:clusters
                Cindx(i, j) = t;   %    t indicates the index of the cell, whereas its 
                Cindx(j, i) = -t;  %    sign signals a reversal of the classifier bound.
                INDXI = traindata.Y == i;
                INDXJ = traindata.Y == j;
                temp_data.Y = cat(1, ones(sum(INDXI), 1), 2*ones(sum(INDXJ), 1));
                temp_data.X = cat(1, traindata.X(INDXI, :), traindata.X(INDXJ, :));
                op = options;
                op.class = 2;
                model.model{t} = spectral_svmtrain_(temp_data, op); 
                t = t+1;
            end
        end
        model.Cindx = Cindx;
        model.type = 'MultiClass';

    elseif options.class == 2
        %   two class training
        X = mds_(L);
        %   train the model using given data
        model = svmtrain(X, traindata.Y);
        model.Y = traindata.Y;
        model.type = 'TwoClass';        
    else
        [V, D] = eig(L);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Recover the second smallest eigenvalue of the Laplacian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~, IDX] = min(diag(D));
        D(IDX, IDX) = max(diag(D));
        [~, IDX] = min(diag(D));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Put the model together
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        model.Y = (sign(V(:, IDX))/2+1.5);     %Use a 1, 2 label set
        model.Alpha = abs(V(:, IDX));   
        model.Bias = 0;
        model.type = 'OneClass';
        model.SupportVectorIndices = (1:length(model.Y))';
    end

    model.X = traindata.X;
    model.compute_graph = options.compute_graph;
    model.t = options.t;
    model.h = options.h;

end
