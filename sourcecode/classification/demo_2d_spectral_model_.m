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

function demo_2d_spectral_model_(model)


    if strcmpi(model.type, 'oneclass') || strcmpi(model.type, 'twoclass')
        %In case the model is given as an adjacency matrix, compute the MDS
        %coordinates
        [X, Y] = meshgrid(min(model.X(:, 1)):(max(model.X(:, 1))-min(model.X(:, 1)))/33:max(model.X(:, 1)), ...
                 min(model.X(:, 2)):(max(model.X(:, 2))-min(model.X(:, 2)))/33:max(model.X(:, 2)));
        if ~model.compute_graph
            if ~isfield(model, 'W')
                error('input data structure does not contain Weigth matrix');
            end
            trn.X = mds(model.W, 2);
        end

        %Do the meshgrid

        testdata.X = cat(2, reshape(X, length(X)^2, 1), reshape(Y, length(Y)^2, 1));
        %Compute the adjacency matrix for the testing coordinates
        TRI = delaunayn(testdata.X);
        testdata.W = get_delaunay_weight_(TRI, testdata.X, model.h);
        %Do the labelling and probability computation
        [PRED, PROB] = spectral_svmtest_(testdata, model);
        %Compute the final c-map

        y = zeros(1, length(PRED));
        y(PRED == 2) = 1;
        y(PRED == 1) = -1;
        figure;
        h = pcolor(X, Y, reshape(PROB(1:length(X)^2), size(X)));
        set(h, 'EdgeColor', 'none')
        figure;
        h = pcolor(X, Y, reshape(y(1:length(X)^2), size(X)));
        set(h, 'EdgeColor', 'none')
        figure;
        h = pcolor(X, Y, reshape(y(1:length(X)^2), size(X))+reshape(PROB(1:length(X)^2), size(X))/2);
        set(h, 'EdgeColor', 'none')
    else
        if ~isfield(model, 'Cindx')
            error('Missing field - Cindx from mode');
        end
        
        [cluster_num, ~] = size(model.Cindx);
        min_x1 = min(model.model{1}.X(:, 1));
        max_x1 = max(model.model{1}.X(:, 1));
        min_x2 = min(model.model{1}.X(:, 2));
        max_x2 = max(model.model{1}.X(:, 2));
        
        for i = 2:cluster_num
            if min(model.model{i}.X(:, 1))<min_x1
                min_x1=min(model.model{i}.X(:, 1));
            end
            if max(model.model{i}.X(:, 1))>max_x1
                max_x1=max(model.model{i}.X(:, 1));
            end
            if min(model.model{i}.X(:, 2))<min_x2
                min_x2=min(model.model{i}.X(:, 2));
            end
            if max(model.model{i}.X(:, 2))>max_x2
                max_x2=max(model.model{i}.X(:, 2));
            end
        end    
        %Do the meshgrid

        [X, Y] = meshgrid(min_x1:(max_x1-min_x1)/100:max_x1, min_x2:(max_x2-min_x2)/100:max_x2);

        testdata.X = cat(2, reshape(X, [length(X)^2, 1]),reshape(Y, [length(Y)^2, 1]));
%         TRI = delaunayn(testdata.X);
%         testdata.W = get_delaunay_weight_(TRI, testdata.X, model.h);
        %Do the labelling and probability computation
        [PRED, PROB] = spectral_svmtest_(testdata, model);
        
        %Compute the final c-map
        figure;
        h = pcolor(X, Y, reshape(PROB(1:length(X)^2), size(X)));
        set(h, 'EdgeColor','none')
        figure;
        h = pcolor(X, Y, reshape(PRED(1:length(X)^2), size(X)));
        set(h, 'EdgeColor','none')
        figure;
        h = pcolor(X, Y, reshape(PRED(1:length(X)^2), size(X))+reshape(PROB(1:length(X)^2), size(X))/2);
        set(h, 'EdgeColor','none')
    end

end