% function [SHADE, DATA, WEIGHT] = shape_from_shading_(I, MASK, debug)
%
%% Syntax:
%    [SHADE, DATA, WEIGHT] = shape_from_shading(I, MASK, debug) or
%    [SHADE, DATA, WEIGHT] = shape_from_shading(I)
% 
%% Description:
%    Recover shapes of objects using only shading information. The algorithm used is half blind signal separation. 
% 
%% Input:
%    I:     Hyperspectral data cube of size [height, width, bands].
%    MASK:  the foreground matrix (optional) where 0 for background and 1 for foreground. 
%           Alternatively, if a single numeric value v is given, it will be used as the threshold 
%           of calculating foreground. In that case, MASK = mean(I, 3) > v;
%    debug: Level of debugging information. Ranges from 1 to 3 while 1 is minimum
%    
%% Output:
%    SHADE:  result of shape estimated.
%    DATA:   data structure, it also contains MASK used in this function
%    WEIGHT: The MASK used (or generated) in this function. If MASK is given as an input, its output will be the same;
% 
%%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013-2014 All Rights Reserved.
% Author: Ran Wei and Lin Gu
% Version: 1.0.1
% Last Update Date: 1 Sep 2014

function [SHADE, DATA, weight] = shape_from_shading_(I, MASK, debug)

    %    check input arguments
    if ~exist('debug', 'var') || numel(debug) ~= 1 || debug(1) < 0
        debug = 1; %    debug level 1 to 3. 1 is minimum
    end

    % get the dimensions of input data cube
    [height_ width_ bands_] = size(I);
    frame_size_ = height_*width_;
    
    disp('Start to estimate shape from shading');
    
    if debug > 2
        s = sprintf('Input data cube dimensions: height = %d, width = %d, bands = %d', height_, width_, bands_);
        disp(s);
    end
    
    if ~exist('MASK', 'var') 
    %   mask is missing, then we need to estimate the foreground.  use 1 here. later it should be set to
    %   other formula
        threshold_ = (max(I(:)) - min(I(:)))*0.05 + min(I(:));
        MASK = mean(I, 3) > threshold_;
        if debug > 2
            s = sprintf('Foreground threshold is %d', threshold_);
            disp(s);
        end
    else 
        if isnumeric(MASK)    %    when MASK is given,
            [mh_, mw_, mb_] = size(MASK);
            if mh_ * mw_ * mb_ == 1
                % a threshold of mask is given
                if debug > 2
                    s = sprintf('Foreground threshold is %d', MASK);
                    disp(s)
                end
                MASK = mean(I, 3) > MASK;
            else
                if mh_ ~= height_ || mw_ ~= width_ || mb_ ~= bands
                    error('Size of argument MASK is not correct. Please check.');
                end
            end
        else
            error('Input MASK is not correct. Please check.');
        end
    end

    if ~isfloat(I)
        % if input is not of float-point type, convert it.
        I = double(I);
        if debug > 3
            disp('Converted input data to double');
        end
    end
    % now we have all input argument ready
    disp('Start to estimate relative shading\n');
    
    %     minus average values from I
    chro_       = I - repmat(mean(I, 3), [1 1 bands_]);
    std_        = std(chro_, 0, 3);
    std_(~MASK) = 1;
    std_        = std_(:);
    RY          = reshape((std_ ./ mean(std_)) .* MASK(:), [], 1)/10;

    if debug >= 2 
        disp('Finished to estimate relative shading');
    end
    %   tested, code generates identical results to original code

    % generate A and D matrix
    if debug >= 2
        disp('Start to compute peripheral matrix');
    end    
    sample_number_ = 40;
    sample_range_  = 5;

    if debug >= 3
        disp('Start randomly selecting coordinates');
    end
    % random sample the pixels either from neighbour
    [xt_ yt_]  = ind2sub([height_ width_], 1 : frame_size_);
    per_xt_    = reshape(repmat(xt_, [sample_number_ 1]), [], 1) + randi(sample_range_, sample_number_ * frame_size_, 1) - ceil(sample_range_ / 2);
    per_xt_    = min(max(per_xt_, 1), height_);
    x_sampled_ = double(reshape(repmat(uint32(1 : frame_size_), sample_number_, 1),[],1));
    per_yt_    = reshape(repmat(yt_, [sample_number_ 1]), [], 1) + randi(sample_range_, sample_number_ * frame_size_, 1) - ceil(sample_range_ / 2);
    per_yt_    = min(max(per_yt_, 1), width_);
    y_sampled_ = double(sub2ind([height_, width_], per_xt_, per_yt_));
    f_dis_     = [MASK(:)' * 100 ; [xt_; yt_] / sqrt(height_^2 + width_^2)+ rand(2, frame_size_) * 1e-6];
    normal_I_  = sqrt(sum(I .^ 2, 3));
    normal_I_(normal_I_ == 0) = eps;
    normalised_I_ = reshape((I ./ repmat(normal_I_, [1 1 bands_])), [], bands_);
    weight_ref_   = 1 - (acos(min(sum(normalised_I_(x_sampled_, :) .* normalised_I_(y_sampled_, :), 2), 1)) / pi * 180)/100.0;
    weight_ref_(weight_ref_ < 0.9) = 0;
    weight_d_ = (max(1 - sum(abs(f_dis_(2 : 3, x_sampled_) - f_dis_(2 : 3, y_sampled_))) / 0.1, 0))';

    if debug >= 3
        disp('Finished selecting coordinates');
    end

    sampled_I_ = RY;
    sampled_I_(sampled_I_ == 0) = 1;
    A_coeff_ = weight_d_ ./ (sampled_I_(x_sampled_) .* (1 - weight_ref_) + weight_ref_ .* sampled_I_(y_sampled_));
    A_coeff_(isnan(A_coeff_)) = 1;

    A_sparsed_ = sparse(double(x_sampled_), double(y_sampled_), weight_d_, frame_size_, frame_size_);
    Aw_coeff_  = full(sum(A_sparsed_, 2) ./ sampled_I_);
    Aw_coeff_(isnan(Aw_coeff_)) = 1;

    Dwp_ = spdiags(Aw_coeff_, 0, frame_size_, frame_size_);
    Ap_  = sparse(double(x_sampled_), double(y_sampled_), A_coeff_, frame_size_, frame_size_);

    if debug >= 2
        disp('Finished to compute peripheral matrix');
    end
    %   Ap_ and Dwp_ have been tested. They are of same value as original code does, from peripheral_matrix_v4.m
    
    rand_i_list_ = find(MASK(:));
    y_sampled_ = reshape(rand_i_list_(randi(length(rand_i_list_), sample_number_, frame_size_)), [], 1);
    y_sampled_ = double(y_sampled_);
    x_sampled_ = reshape(repmat(uint32(1 : frame_size_), sample_number_, 1), [], 1);
    x_sampled_ = double(x_sampled_);
    normal_I_  = sqrt(sum(I .^ 2, 3));
    normal_I_(normal_I_ == 0) = eps;
    normalised_I_ = reshape((I ./ repmat(normal_I_, [1 1 bands_])), [], bands_);
    weight_ref_   = 1 - (acos(min(sum(normalised_I_(x_sampled_, :) .* normalised_I_(y_sampled_, :), 2), 1)) / pi * 180)/100.0;
    weight_ref_(weight_ref_ < 0.9) = 0;
    
    A_coeff_ = weight_ref_ ./ (sampled_I_(x_sampled_) .* (1 - weight_ref_) + weight_ref_ .* sampled_I_(y_sampled_));
    A_coeff_(isnan(A_coeff_)) = 1;
    Ag_ = sparse(double(x_sampled_), double(y_sampled_), A_coeff_, frame_size_, frame_size_);
    
    Aw_coeff_ = full(sum(sparse(double(x_sampled_), double(y_sampled_), weight_ref_, frame_size_, frame_size_), 2) ./ sampled_I_);
    Aw_coeff_(isnan(Aw_coeff_)) = 1;
    Dwg_  = spdiags(Aw_coeff_, 0, frame_size_, frame_size_);
    D_b_  = spdiags(double(MASK(:)), 0, frame_size_, frame_size_);
    D_bc_ = spdiags(double(~MASK(:)), 0, frame_size_, frame_size_);

    % mask out all the pixels outside the foreground
    A  = Ap_ + Ag_;
    Dw = Dwp_ + Dwg_;

    DATA.A    = D_b_ * A + D_bc_;
    DATA.Dw   = D_b_ * Dw + D_bc_;
    DATA.sI   = RY;
    DATA.MASK = MASK;

    weight.BL    = 10;
    weight.Shade = 0;
    weight.Prec  = 0;
    weight.sm    = 10;
    weight.maxit = 200;

    SHADE = get_shade_(I, DATA, 1, weight);
    
%    end of function
end