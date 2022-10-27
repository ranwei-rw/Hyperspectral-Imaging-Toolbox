% Syntax
%   model = train_gmm_model(I);
%   model = train_gmm_model(I, initmodel);
%   model = train_gmm_model(I, initmodel, options);
%
% Description:
%   This function is used to recover the General Mixture Model (GMM) based on the paper by 
%   Ravindran Kannan et al. entitled "The Spectral Method for General
%   Mixture Models" (COLT 2005).
% 
% Input: 
%   I: the radiance image (stored as a 3D array with size height x width x bands).
%   initmodel: Initial value of the general Mixture Model.
%   options: Structure containing the following fields
%       use_gaussian: If this variable is set ot unity, the GMM model is
%           recovered using a Gaussian prior. If this field is set to zero, a
%           kernel density estimator (KDE) is used instead. The default is
%           options.use_gaussian = 0.
%       ndims: Number of dmensions in which the data is to be subspace-projected.
%           If ndims is not provided, it is automatically set to 6.
%       ngmm: Number of mixtures used for the initialisation. This is only
%           used if no initmodel is provided and, by default, is set to 3;
% 
% Output: 
%   model: Structure containing the initial general mixture model. This contains the following fields
%           q: Number of components for the mixture
%           T_idx: Cell array containing the indeces for the pixels classified based 
%             upon the distance to the mixture means.
%           centres_w: Means for the weighted mixture weights
%           covars_w: Covariance matrices for the weighted mixtures.
%           centres: Means for the unweighted mixture weights
%           covars: Covariance matrices for the unweighted mixtures.
%           Uq: Principal components for the sampled pixels. This is computed using SVD.
%           EM_options: Set of options used for the em_gmm routine
%           prioL: Prior for each of the components over the data (image pixels).
%           tip_record: Value of the pixel (instance) in the data which corresponds to the simplex 
%             that is better aligned to the mixture.
%           pLx_map: Map of labels for the data accroding to the GMM and the
%             Markov random field used for the train_gmm_model.
%           pLx_prob: Map of probabilities for the data as delivered by the
%             Markov random field.
%
%   See also:
%
%       init_gmm_model, em_kde, em_gmm, randwalk
%
%   Notes:
%
%   This routine requires the use of the randwalk function and, hence, 
%   the SuiteSparse package, which can be downloaded from 
%   http://www.cise.ufl.edu/research/sparse/SuiteSparse/
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Lin Gu
% Version: 1.0.1
% Last Update Date: 15 Jan 2014

function model = train_gmm_model_(I, initmodel, options)

I = double(I);
[height_, width_, bands_] = size(I);

if ~exist('options', 'var') || ~isfield(options,'use_gaussian')
    options.use_gaussian = 0;
end
if ~isfield(options,'ndims')
    if bands_>6
        options.ndims = 6;
    else
        options.ndims = bands_;
    end
end
if ~isfield(options,'ngmm')
    options.ngmm = 3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup the variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do pre segment function pre_segment_multi(I)
SUM_I_        = sum(I, 3);
frame_        = height_ * width_;
SUM_I_1D_     = reshape(SUM_I_, frame_, 1);
SORTED_I_1D_  = sort(SUM_I_1D_);
thresh_low_   = SORTED_I_1D_(ceil(0.05 * frame_));
thresh_high_  = SORTED_I_1D_(ceil(0.95 * frame_));
SORTED_INDEX_ = (SUM_I_1D_ < thresh_high_) & (SUM_I_1D_ > thresh_low_);

clear SORTED_I_1D_;
clear SUM_I_1D_;
clear SUM_I_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise the model if no initmodel is provided
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('initmodel', 'var') || isempty(initmodel)       
        initmodel = init_gmm_model_(I, SORTED_INDEX_, options.ndims, options.ngmm);                     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the refinement of the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.use_gaussian
    model = em_gmm_(I, initmodel, SORTED_INDEX_);
else
    model = em_kde_(I, initmodel, SORTED_INDEX_);
end