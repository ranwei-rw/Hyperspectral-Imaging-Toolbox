% Syntax
%       model = init_gmm_model(I);
%       model = init_gmm_model(I, INDEX);
%       model = init_gmm_model(I, INDEX, ndims);
%       model = init_gmm_model(I, [], ndims);
%       model = init_gmm_model(I, INDEX, ndims, ngmm);
%       model = init_gmm_model(I, [], ndims, ngmm);
%       
% Description:
%   This function is used to initialise the General Mixture Model based on the paper by 
%   Ravindran Kannan et al. entitled "The Spectral Method for General
%   Mixture Models" (COLT 2005).
% 
% Input: 
% 
%   I:     the radiance image (stored as a 3D array with size height x width x bands).
%   INDEX: Indices for the pixels in the image to be sampled for the estimation process. 
%          If INDEX is empty, all pixels are considered to be valid.
%   ndims: Number of dmensions in which the data is to be subspace-projected. It is set, by default, to six.
%   ngmm:  Number of mixtures used for the initialisation. This is only used if no initmodel is provided and, by
%          default, is set to 3; 
% 
% Output: 
%
%   model: Structure containing the initial general mixture model. This contains the following fields:
%          ngmm:       Number of components for the mixture
%          T_idx:      Cell array containing the indeces for the pixels classified based upon the distance to the
%                      mixture means. 
%          centres_w:  Means for the weighted mixture weights
%          covars_w:   Covariance matrices for the weighted mixtures.
%          centres:    Means for the unweighted mixture weights
%          covars:     Covariance matrices for the unweighted mixtures.
%          Uq:         Principal components for the sampled pixels. This is computed using SVD.
%          EM_options: Set of options used for the em_gmm routine
%          prioL:      Prior for each of the components over the data (image pixels).
%          tip_record: Value of the pixel (instance) in the data which corresponds to the simplex that is better aligned
%                      to the mixture. 
%
% See also:
%       train_gmm_model
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Lin Gu, Ran Wei, Antonio Robles-Kelly
% Version: 1.0.2
% Last Update Date: 29 July 2014
    
function model = init_gmm_model_(I, INDEX, ndims, ngmm)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup the variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [height_, width_, bands_] = size(I);
    
    if ~exist('ngmm', 'var')    
        ngmm = bands_;
    end
    
    if ~exist('ndims', 'var') || numel(ndims) > ngmm 
       if bands_ > 6
           ndims = 6;
       else
           ndims = bands_;
       end
    end
    
    if ~exist('INDEX', 'var')
        INDEX = 1 : height_*width_;
    end
    
    ep_ = 1/ngmm;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape input data into 2D mode
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_2D_ = reshape(I, height_*width_, bands_);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove zero elements according to Indx
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_2D_(~INDEX, :)   = [];
    I_INDEXES_         = 1:height_*width_;
    I_INDEXES_(~INDEX) = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first sample a subset S to find the SVD then discard it, however, it is
    % useless to do this step for the RGB image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if bands_ == 3
        Uq_ = diag(ones(3, 1));
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % use about 10% valid pixels of data to do SVD.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        valid_pixel_number_       = max(INDEX);
        SAMPLE_INDEX_             = unique(ceil(1 + rand(ceil(0.1*valid_pixel_number_), 1)*(valid_pixel_number_-1)));
        I_SAMPLED_                = I_2D_(SAMPLE_INDEX_, :);
        [U_, ~, ~]                = svd(I_SAMPLED_');
        Uq_                       = U_(:, 1:ndims);
        I_2D_(SAMPLE_INDEX_, :)   = [];% get rid of points that already being sampled
        I_INDEXES_(SAMPLE_INDEX_) = [];
    end
    
    PROJECTED_SUBSPACE_ = I_2D_*Uq_;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check whether the separation is enough
    % randomly choose some points in set to form the clustering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k_ = 1:ngmm
        
        pj_height_ = size(PROJECTED_SUBSPACE_, 1); 
        
        CLUSTER_SAMPLE_INDEX_ = unique(ceil(1 + rand(ceil(0.02*pj_height_), 1)*(pj_height_ - 1)));
        
        for m_ = 1 : length(CLUSTER_SAMPLE_INDEX_)
            POINT_              = PROJECTED_SUBSPACE_(CLUSTER_SAMPLE_INDEX_(m_), :);
            distance_           = sum((PROJECTED_SUBSPACE_ - repmat(POINT_, [pj_height_, 1])) .^ 2, 2);
            sorted_distance_    = sort(distance_);
            distance_threshold_ = sorted_distance_(ceil(ep_ * k_ * pj_height_ / 2));
            CHOPPED_ITW_        = PROJECTED_SUBSPACE_(distance_ < distance_threshold_, :);
            CHOPPED_ITW_        = CHOPPED_ITW_ - repmat(mean(CHOPPED_ITW_), [size(CHOPPED_ITW_, 1), 1]);
            pccoeff             = princomp(CHOPPED_ITW_);
            POINT_SIGMA_(m_)    = std(CHOPPED_ITW_ * pccoeff(:, 1));
        end
       
        [~, POINT_SIGMA_INDEX_] = sort(POINT_SIGMA_);
        point_index_            = POINT_SIGMA_INDEX_(ceil(length(POINT_SIGMA_INDEX_) * 0.98));
        clear POINT_SIGMA_;

        if k_ == ngmm
            [ITW_IDX_, C_, SUMD_] = kmeans(PROJECTED_SUBSPACE_, ngmm, 'EmptyAction', 'singleton', 'OnlinePhase', 'off');
            
            for n_ = 1 : max(ITW_IDX_) 
                KMEAN_SUM_(n_) = sum(ITW_IDX_ == n_);
            end
            
            SUMD_          = SUMD_ ./ KMEAN_SUM_';
            [~, SUMD_IDX_] = sort(SUMD_);
            TIP_POINT_     = C_(SUMD_IDX_(1), :);
        else
            TIP_POINT_ = PROJECTED_SUBSPACE_(CLUSTER_SAMPLE_INDEX_(point_index_), :);
        end
        
        TIP_RECORD_(k_, :)          = TIP_POINT_;
        PROJECT_T_                  = PROJECTED_SUBSPACE_ * null(TIP_POINT_);        
        NULL_POINT_DISTANCE_        = sum(PROJECT_T_ .^ 2,2);
        SORTED_NULL_POINT_DISTANCE_ = sort(NULL_POINT_DISTANCE_);
        
        if ep_ < 0.4
            DISTANCE_WITHIN_ = SORTED_NULL_POINT_DISTANCE_(ceil(ep_ * 2 * pj_height_));
        else
            DISTANCE_WITHIN_ = SORTED_NULL_POINT_DISTANCE_(min(ceil(ep_ * 1.5 * pj_height_), pj_height_));
        end
        
        t_distance_         = SORTED_NULL_POINT_DISTANCE_(ceil(ep_ * k_ * 0.5 * pj_height_));
        ITW_CHOPPED_        = PROJECTED_SUBSPACE_(NULL_POINT_DISTANCE_ < t_distance_, :);
        T_IDX_{k_}          = I_INDEXES_(NULL_POINT_DISTANCE_ < t_distance_);
        T_MEAN_(k_, :)      = mean(ITW_CHOPPED_);        
        T_COVARI_(:, :, k_) = cov(ITW_CHOPPED_);
        PROJECTED_SUBSPACE_(NULL_POINT_DISTANCE_ < DISTANCE_WITHIN_, :) = [];
        I_INDEXES_(NULL_POINT_DISTANCE_ < DISTANCE_WITHIN_) = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish off with the model parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SUM_I_          = reshape(sum(I, 3), height_*width_, 1);
    q_threshold_low = min(SUM_I_);
    I_2D_           = reshape(I, height_*width_, bands_);   %   Re-generate 2D input data cube

    for p_ = 1 : ngmm
        T_IDX_{p_}(SUM_I_(T_IDX_{p_}) < q_threshold_low) = [];
        T_MEAN_FULL_(p_, :)      = mean(I_2D_(T_IDX_{p_}, :));
        T_COVARI_FULL_(:, :, p_) = cov(I_2D_(T_IDX_{p_}, :));
        L_PRIO_(p_)              = length(T_IDX_{p_});
    end

    L_PRIO_                          = L_PRIO_ / sum(L_PRIO_);
    model.q                          = ngmm;
    model.T_idx                      = T_IDX_;
    model.centres_w                  = T_MEAN_;
    model.covars_w                   = T_COVARI_;
    model.centres                    = T_MEAN_FULL_;
    model.covars                     = T_COVARI_FULL_;
    model.Uq                         = Uq_;
    model.EM_options.illustrate_on   = 0;
    model.EM_options.PrioL_on        = 1;
    model.EM_options.MRF_on          = 1;
    model.EM_options.niters          = 3;
    model.EM_options.LowDimension_on = 1;
    model.EM_options.Bounddetect_on  = 1;
    model.prioL                      = L_PRIO_;
    model.tip_record                 = TIP_RECORD_;

end