% Syntax
%
%   [Labels, Prob] = randwalk(G);
%
% Description:
%   Implementation of the generalised random walk with side information for labelling as presented by Z. Fu and A Robles-Kelly in the paper
%   entitled "A Quadratic Programming Approach to Image Labeling" (IET Computer Vision, 2008).
% 
% Input: 
% 
%   G: Structure containing the graph. Note that idxCls may be empty (this allows for the use of I2Graph with a fixed
%      number of clusters (G.nCls). This structure includes the following members:
%
%       width, height: width and height of the input image.
%       ii, jj, vv:    Triplet for the graph edge term. These is such that the weighted adjacency matrix of the graph is
%                      given by W = sparse(G.ii, G.jj, G.vv, N, N), where N=G.width*G.height;
%       e:             Data term with the probability of the vertex belonging to the cluster. This is an N x c array
%                      where N is the number of vertices as before and c the number of classes.
%       idxCls:        Cell of order c (number of classes) where idxCls{c} contains the indeces for the labelled pixels
%                      in the class c. 
%       mCls:          Mean vectors for the labelled pixels. This is a nCls x bands array where the row index
%                      corresponds to the class. 
%       nCls:          Number of labelled classes in the mask file.
% 
% Output: 
%
%   LABELS: Index of the most likely label for each entry entry.
%   PROB:   Vector whose entry indexed i, c denotes the probability of the node i belonging to the class c. 
%   
%
% Example:
%
%       imdata = imread('ngc6543a.jpg');
%       mask = imread('ngc6543a_mask.tif');
%       G = Mask2Graph(imdata, mask);
%       [Labels, Prob] = randwalk(G);
%       imtool(reshape(Labels, G.height, G.widthh)/G.nCls);
%
%   See also:
%
%       I2Graph, Data2Graph
%
%   Notes: 
%
%   This routine requires the SuiteSparse package, which can be downloaded from 
%   http://www.cise.ufl.edu/research/sparse/SuiteSparse/
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014-2015 All Rights Reserved.
%  Author: Zhouyu Fu, Antonio Robles-Kelly and Ran Wei
%   Version: 1.0.3
%   Last update date: 19 June 2015

%  Version: 1.0.2
%  Last Update Date: 4 Aug 2014

function [LABELS, PROB] = randwalk_(G)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Setup the variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N   = G.height * G.width;
    G.e = 1-G.e;

    if size(G.e, 2) ~= G.nCls
        error('data on G must be consistent with the number of classes');
    end
    
    if  ~isvector(G.ii) || ~isvector(G.jj) || ~isvector(G.vv) ...
        || length(G.ii)~=length(G.jj) || length(G.ii)~=length(G.vv)
        error('ii, jj, vv must be vectors with the same dimension');
    end
    
    idx_cls_present = 1;
    if ~isfield(G, 'idxCls')
        idx_cls_present = 0;
%         for i = 1:G.nCls
%             G.idxCls{i} = [];
%         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Perform the computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W     = sparse(G.ii, G.jj, G.vv, N, N);
    W     = W + W';
    D     = sparse((1:N)', (1:N)', full(sum(W, 2)));
    L     = D - W;
    IDX_U = (1:N)';
    IDX_L = [];
    YL    = [];
    if idx_cls_present
        for i = 1:G.nCls
            IDX_L = [IDX_L;G.idxCls{i}];
            YL    = [YL; i * ones(length(G.idxCls{i}), 1)];
        end

        IDX_U(IDX_L) = [];
        YL = lblvec2mat_(YL, G.nCls);
    end

    nunlabel = length(IDX_U);
    sume     = sum(G.e, 2);
    C        = cs_sparse((1:nunlabel)', (1:nunlabel)', sume(IDX_U));
    E        = sume*ones(1, G.nCls) - G.e;
    if idx_cls_present
        E = E(IDX_U, :) - L(IDX_U, IDX_L)*YL;
    else
        E = E(IDX_U, :);
    end
    
    for i = 1:G.nCls
        LABELS(:, i) = cs_cholsol(L(IDX_U, IDX_U)+C, E(:, i));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get the probabilities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PROB           = zeros(N, G.nCls);
    LABELS         = (LABELS-min(LABELS(:)))/(max(LABELS(:))-min(LABELS(:)));
    PROB(IDX_U, :) = LABELS;
    
    if idx_cls_present
        PROB(IDX_L, :) = YL;
    end

    PROB = PROB./repmat(sum(PROB, 2), 1, G.nCls);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get the labels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, LABELS] = max(PROB, [], 2);    

%   end of function randwalk_
end