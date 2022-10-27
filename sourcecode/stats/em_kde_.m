%%  Syntax
%       model = em_kde(I);
%       model = em_kde(I, init_model);
%       model = em_kde(I, init_model, INDEX);
%       
%%   Description:
%       This function is used to compute the General Mixture Model based on the paper by 
%       Ravindran Kannan et al. entitled "The Spectral Method for General
%       Mixture Models" (COLT 2005) using the EM algorithm. It employs at input the model delivered by
%       the routine init_gmm_model and a prior compted via the application of kernel
%       density estimation (KDE).
% 
%%   Input: 
%       I:          The radiance image (stored as a 3D array with size height x width x bands).
%       INDEX:      Indices for the valid pixels in the image. If Indx is not provided, all pixels are considered to be valid.
%       init_model: Initial model used for the General Mixture. If this is not provided, the method will employ
%                   init_gmm_model to recover one. 
% 
%%   Output: 
%       model: Structure containing the final general mixture model. It has the following members.
%              q:          Number of components for the mixture
%              T_idx:      Cell array containing the indeces for the pixels classified based upon the distance to the
%                          mixture means. 
%              centres_w:  Means for the weighted mixture weights
%              covars_w:   Covariance matrices for the weighted mixtures.
%              centres:    Means for the unweighted mixture weights
%              covars:     Covariance matrices for the unweighted mixtures.
%              Uq:         Principal components for the sampled pixels. This is computed using SVD.
%              EM_options: Set of options used for the em_gmm routine
%              prioL:      Prior for each of the components over the data (image pixels).
%              tip_record: Value of the pixel (instance) in the data which corresponds to the simplex that is better
%                          aligned to the mixture. 
%              pLx_map:    Map of labels for the data accroding to the GMM and the Markov random field used for the
%                          train_gmm_model. 
%              pLx_prob:   Map of probabilities for the data as delivered by the Markov random field.
%
%%   See also:
%       init_gmm_model, trian_gmm_model, em_gmm, randwalk
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

function model = em_kde_(I, init_model, INDEX)

    [height_, width_, bands_] = size(I);

    if ~exist('INDEX', 'var')
        INDEX = 1:height_*width_;
    end

    if ~exist('init_model', 'var')
        init_model = init_gmm_model_(I, INDEX);
    end
        
    if bands_ == 3
        CENTERS_ = init_model.centres_w;
        COVARS_  = init_model.covars_w;
    else
        if init_model.EM_options.LowDimension_on == 1
            CENTERS_ = init_model.centres_w;
            COVARS_  = init_model.covars_w;
            I        = reshape(I, height_*width_, bands_)*init_model.Uq;
            bands_   = size(I, 2);
            I        = reshape(I, height_, width_, bands_);
        else
            CENTERS_ = init_model.centres;
            COVARS_  = init_model.covars;
        end
    end

    framesize_ = height_*width_;
    PRIO_L_    = repmat(init_model.prioL, [framesize_, 1]);
    I_2D_      = reshape(I, framesize_, bands_);
    PXL_       = zeros(framesize_, init_model.q);

    for i = 1:init_model.q
        DIFF_     = I_2D_ - repmat(CENTERS_(i, :),[framesize_, 1]);
        PXL_(:, i) = exp(-0.5*sum((DIFF_ / COVARS_(:, :, i)) .* DIFF_, 2)) / ((2*pi) ^ (init_model.q / 2)*det(COVARS_(:, :, i)));
    end

    PL_X_   = normalise_probability_(PXL_, PRIO_L_);
    PL_X_   = reshape(PL_X_, height_, width_, init_model.q);
    PRIO_L_ = spatial_prio_(PL_X_, INDEX, init_model.EM_options.PrioL_on);
    [PL_X_, init_model.q, PRIO_L_] = shrink_cluster_(PL_X_, init_model.q, PRIO_L_);

    [~, q] = size(PL_X_);

    I_2D_ = reshape(I, height_*width_,[]);

    if q > 1
        for n = 1 : init_model.EM_options.niters
            % E step,  compute the expectation of each points
            % Here use the P(x|L) comes from the KDE
            q = size(PL_X_, 3);

            if q ~= 1
                PL_X_ = reshape(PL_X_, [], q);
            else
                q = size(PL_X_, 2);
            end

            PL_X_F_ = PL_X_(logical(INDEX), :);
            I_2D_F_ = I_2D_(logical(INDEX), :);
            PX_L_   = zeros(size(PL_X_));

            for i = 1 : q
                P_KDE_B_ = zeros(size(PL_X_));
                MIN_     = min(I_2D_, [], 1);
                MAX_     = max(I_2D_, [], 1);

                for j = 1:size(I_2D_, 2)
                    PROB_TABLE_    = ksdensity(I_2D_F_(:, j), MIN_(j) : (MAX_(j) - MIN_(j)) / 1000 : MAX_(j), ...
                                               'weight', PL_X_F_(:, i));
                    QUAN_          = max(ceil((I_2D_(:, j) - MIN_(j)) / (MAX_(j) - MIN_(j)) * 1000), 1);
                    P_KDE_B_(:, j) = PROB_TABLE_(QUAN_);
                end
                
                PX_L_(:, i) = prod(P_KDE_B_, 2);
            end

            PX_L_(isnan(PX_L_)) = 1/q;
            
            PL_X_ = normalise_probability_(PX_L_, PRIO_L_);

            if init_model.EM_options.MRF_on == 1
                %   MRF is implemented for every E step
                PL_X_MAP_           = reshape(PL_X_, height_, width_, q);
                [PL_X_MAP_,PRIO_L_] = mrf_prob_(PL_X_MAP_);
                PRIO_L_             = reshape(PRIO_L_,height_ * width_, q);
            else
                % Otherwise use a simple spatial prior estimation
                PRIO_L_ = spatial_prio_(PL_X_MAP_, INDEX, init_model.EM_options.PrioL_on);
            end        
                
            % M step       
 
            PL_X_     = reshape(PL_X_MAP_, height_*width_, q);             
            PL_X_CAT_ = PL_X_(INDEX, :);
            NEW_PR_   = sum(PL_X_CAT_);

            if ~prod(NEW_PR_)
                break;
            end

            NEW_CENTERS_   = PL_X_CAT_' * I_2D_(INDEX, :) ./ (NEW_PR_' * ones(1, size(I_2D_, 2)));

            if sum(size(NEW_CENTERS_(:))) == sum(size(CENTERS_(:))) && ...
                    sum(abs(NEW_CENTERS_ - CENTERS_)) / sum(abs(NEW_CENTERS_)) < 0.001
                break;
            else
                CENTERS_ = NEW_CENTERS_;
            end

        end

    end

    % Eliminate the smallest cluster
    [PL_X_, q, PRIO_L_] = shrink_cluster_(PL_X_, q, PRIO_L_);

    if q < 2
        PL_X_MAP_ = ones(height_, width_);
        PRIO_L_   = ones(height_, width_);
    else
        PL_X_MAP_ = reshape(PL_X_, height_, width_, q);
        PRIO_L_   = reshape(PRIO_L_, height_, width_, q);
    end

    
    init_model.pLx_map  = PL_X_MAP_;
    init_model.pLx_prob = PRIO_L_;
    init_model.q        = q;
    init_model.n        = n;
    init_model.pLx_prob = init_model.pLx_prob./repmat(sum(init_model.pLx_prob(:, :, :), 3), [1 1 init_model.q]);
    model               = init_model;
    close all hidden;
    
end