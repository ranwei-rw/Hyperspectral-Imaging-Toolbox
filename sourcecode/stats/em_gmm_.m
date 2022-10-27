%%  Syntax
%       model = em_gmm(I);
%       model = em_gmm(I, init_model);
%       model = em_gmm(I, init_model, INDEX);
%       
%%  Description:
%       This function is used to compute the General Mixture Model based on the paper by 
%       Ravindran Kannan et al. entitled "The Spectral Method for General
%       Mixture Models" (COLT 2005) using the EM algorithm. It employs at input the model delivered by
%       the routine init_gmm_model and a Gaussian prior.
% 
%%  Input: 
%       I:          the radiance image (stored as a 3D array with size height x width x bands).
%       init_model: Initial model used for the General Mixture. If this is not provided, the method will employ
%                   init_gmm_model to recover one. 
%       INDEX:      Indices for the valid pixels in the image. If Indx is not provided, all pixels
%                   are considered to be valid.
% 
%%  Output: 
%       model: Structure containing the final general mixture model. It has the following members:
%           q:          Number of components for the mixture
%           T_idx:      Cell array containing the indeces for the pixels classified based 
%                       upon the distance to the mixture means.
%           centres_w:  Means for the weighted mixture weights
%           covars_w:   Covariance matrices for the weighted mixtures.
%           centres:    Means for the unweighted mixture weights
%           covars:     Covariance matrices for the unweighted mixtures.
%           Uq:         Principal components for the sampled pixels. This is computed using SVD.
%           EM_options: Set of options used for the em_gmm routine
%           prioL:      Prior for each of the components over the data (image pixels).
%           tip_record: Value of the pixel (instance) in the data which corresponds to the simplex 
%                       that is better aligned to the mixture.
%           pLx_map:    Map of labels for the data accroding to the GMM and the
%                       Markov random field used for the train_gmm_model.
%           pLx_prob:   Map of probabilities for the data as delivered by the
%                       Markov random field.
%
%   See also:
%
%       init_gmm_model, trian_gmm_model, em_kde, randwalk
%
%   Notes:
%
%   This routine requires the use of the randwalk function and, hence, 
%   the SuiteSparse package, which can be downloaded from 
%   http://www.cise.ufl.edu/research/sparse/SuiteSparse/
%
%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei, Lin Gu and Antonio Robles-Kelly
% Version: 1.0.2
% Last Update Date: 29 July 2014

function model = em_gmm_(I, init_model, INDEX)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [height_, width_, bands_] = size(I);

    if ~exist('INDEX', 'var')
        INDEX = 1:height_*width_;
    end

    if ~exist('init_model', 'var')
        init_model = init_gmm_model_(I, INDEX);
    end

        if bands_ == 3
            CENTRES_   = init_model.centres_w;
            COVARS_    = init_model.covars_w;
            ALTERED_I_ = reshape(I, height_*width_, bands_);    
        else
            if init_model.EM_options.LowDimension_on == 1
                CENTRES_   = init_model.centres_w;
                COVARS_    = init_model.covars_w;
                ALTERED_I_ = reshape(I, height_*width_, bands_);
                I          = ALTERED_I_*init_model.Uq;
                bands_     = size(I, 2);
                I          = reshape(I, height_, width_, bands_);
            else
                CENTRES_   = init_model.centres;
                COVARS_    = init_model.covars;
                ALTERED_I_ = reshape(I, height_*width_, bands_)*init_model.Uq;
            end    
        end

        ALTERED_I_(~INDEX, :) = [];

        q           = init_model.q;
        PRIO_L_     = repmat(init_model.prioL, [height_*width_, 1]);
        I_2D_       = reshape(I, height_*width_, bands_);
        I_2D_CAT_   = I_2D_(INDEX, :);
        number_     = size(I_2D_, 1);
        cat_number_ = size(I_2D_CAT_, 1);
        PLX_Q_      = zeros(q, 1);   %   define P(L|x) = P(L)P(x|L)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 1:init_model.EM_options.niters
        % E step, compute the expectation of each points get the P(L|x) = P(L)P(x|L)
        PX_L_ = zeros(number_, q);

        for i = 1:q
            DIFF_       = I_2D_ - repmat(CENTRES_(i, :), [number_, 1]);
            PX_L_(:, i) = exp(-0.5*sum((DIFF_ / COVARS_(:, :, i)) .* DIFF_, 2)) / ((2*pi) ^ (q/2)*det(COVARS_(:, :, i)));
        end

        PL_X_ = normalise_probability_(PX_L_, PRIO_L_);

        if init_model.EM_options.MRF_on == 1
            % if init_model.EM_options(1) = 1, then MRF is implemented for every E step
            PLX_MAP_      = reshape(PL_X_, height_, width_, q);
            [PLX_MAP_, ~] = mrf_prob_(PLX_MAP_);
            PL_X_         = reshape(PLX_MAP_, height_*width_, q);
        end

        PLX_MAP_ = reshape(PL_X_, height_, width_, q); %  PLX_MAP_ is a 3D matrix        
        % M step
        PLX_CAT_ = PL_X_(INDEX, :);        
        PR_      = sum(PLX_CAT_);

        if init_model.EM_options.PrioL_on == 1
            %   implement the prio probablity
            [X, Y] = meshgrid(1:width_, 1:height_);            
            XY_2D  = [X(:) Y(:)];            
            XY_2D(~INDEX, :) = [];            
            CENTRES_2D_ = PLX_CAT_'*XY_2D ./ (PR_'*ones(1, 2));

            for i = 1 : q
                DIFF_XY_ = (XY_2D - ones(cat_number_, 1)*CENTRES_2D_(i, :)) .* (PLX_CAT_(:, i)*ones(1, 2));
                COV_XY_(:, :, i) = DIFF_XY_'*(DIFF_XY_ .* (PLX_CAT_(:, i)*ones(1, 2))) / PR_(i);
            end

            for i = 1 : q
                [X, Y]        = meshgrid(1 : width_, 1 : height_);                
                X             = X - CENTRES_2D_(i, 1);
                Y             = Y - CENTRES_2D_(i, 2);
                XY_2D         = [X(:) Y(:)];
                PRIO_L_(:, i) = sqrt((det(COV_XY_(:, :, i))))*exp(- 0.5*sum(XY_2D / COV_XY_(:, :, i) .* XY_2D, 2));
            end
            clear DIFF_XY_;
            clear COV_XY_;
            clear CENTRES_2D_;

            PRIO_L_ = PRIO_L_ ./ repmat(max(sum(PRIO_L_, 2), 10 ^ -24), [1 q]);
        else
            PRIO_L_ = ones(height_*width_, 1)*PR_;
        end

        FORMER_CENTERS_ = CENTRES_;
        CENTRES_        = PLX_CAT_'*I_2D_CAT_ ./ (PR_'*ones(1, bands_));
        diff_centres    = sum(abs(CENTRES_ - FORMER_CENTERS_)) / sum(abs(CENTRES_));
        if diff_centres < 0.001            
            break;
        end

        for i = 1 : q
            if mean2(PLX_MAP_(:, :, i)) < 0.05
                PLX_Q_(i) = 1;                
            end
            DIFF_ = (I_2D_CAT_ - (ones(cat_number_, 1)*CENTRES_(i, :))) .* (sqrt(PLX_CAT_(:, i))*ones(1, bands_));
            COVARS_(:, :, i) = (DIFF_'*DIFF_) / PR_(i);
        end
        clear DIFF_;
    end
    
    sprintf('Programme iterated for %d times', n)

    % Eliminate the smallest cluster

    PLX_Q_(find(sum(isnan(CENTRES_), 2))) = 1;
    PLS_NULL_ID_                = find(PLX_Q_ == 1);
    CENTRES_(PLS_NULL_ID_, :)   = [];
    COVARS_(:, :, PLS_NULL_ID_) = [];    
    q = q - length(PLS_NULL_ID_);
    
    if q > 0
        PRIO_L_(:, PLS_NULL_ID_)  = [];
        PLX_CAT_(:, PLS_NULL_ID_) = [];
        PR_(PLS_NULL_ID_)         = [];

        % In this step, assign each pixel with the possibility of certain light
        clear PX_L_;
        for i = 1 : q
            DIFF_       = I_2D_ - repmat(CENTRES_(i, :), [number_, 1]);
            PX_L_(:, i) = exp(-0.5*sum((DIFF_ / COVARS_(:, :, i)) .* DIFF_, 2)) / ((2*pi) ^ (q / 2)*det(COVARS_(:, :, i)));
        end

        PL_X_ = normalise_probability_(PX_L_, PRIO_L_);

        CENTERS_ALTERED_ = PLX_CAT_'*ALTERED_I_ ./ (PR_'*ones(1, size(ALTERED_I_, 2)));

        for i = 1 : q
            DIFF_                 = ALTERED_I_ - (ones(cat_number_, 1)*CENTERS_ALTERED_(i, :));
            DIFF_                 = DIFF_ .* (sqrt(PLX_CAT_(:, i))*ones(1, size(ALTERED_I_, 2)));
            COV_ALTERED_(:, :, i) =  (DIFF_'*DIFF_) / PR_(i); 
        end

        if init_model.EM_options.LowDimension_on == 1
            init_model.centres_w = CENTRES_;
            init_model.covars_w  = COVARS_;
            init_model.centres   = CENTERS_ALTERED_;
            init_model.covars    = COV_ALTERED_;
        else
            init_model.centres   = CENTRES_;
            init_model.covars    = COVARS_;
            init_model.centres_w = CENTERS_ALTERED_;
            init_model.covars_w  = COV_ALTERED_;
        end
        close all hidden;

       [PL_X_, q, PRIO_L_] = shrink_cluster_(PL_X_, q, PRIO_L_);
       PLX_MAP_ = reshape(PL_X_, height_, width_, q);
        if q < 2
            PLX_MAP_ = ones(height_, width_);
            PRIO_L_  = ones(height_, width_);
        else
            [~, LABEL_INDEX_] = max(PLX_MAP_, [], 3);
            PRIO_L_ = reshape(PRIO_L_, height_, width_, q);
            for i = 1:q
                PLX_MAP_(:, :, i) = LABEL_INDEX_ == i;
            end
        end
        %   omitted the part of small_patch_remove_v2

        init_model.pLx_map        = PLX_MAP_;
        init_model.pLx_prob       = PRIO_L_;
        init_model.iterated_times = n - 1;
        init_model.q              = q;
        init_model.n              = n;  
        init_model.pLx_prob       = init_model.pLx_prob./repmat(sum(init_model.pLx_prob(:,:,:),3),[1 1 init_model.q]);
    else
        init_model.pLx_map        = ones(height_, width_);
        init_model.pLx_prob       = ones(height_, width_);
        init_model.iterated_times = 1;
        init_model.q              = 1;
        init_model.n              = 1;
    end

    model = init_model;
    
%   end of function    
end