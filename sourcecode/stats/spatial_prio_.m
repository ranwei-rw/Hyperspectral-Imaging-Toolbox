function PRIO = spatial_prio_(PLX, FILTERED_INDEX_, USE_PRIO)

% check number of inputs

if nargin < 3
    error('Error in input arguments. Please check');
end

if length(size(PLX)) < 3
    error('Input should be a 3D matrix');
end

[height_, width_, q_] = size(PLX);

%   reshape input into 2D matrix

PLX_CAT_ = reshape(PLX, [], q_);
PLX_CAT_(~FILTERED_INDEX_, :) = [];
PR_ = sum(PLX_CAT_);

if(prod(PR_))
    [X, Y] = meshgrid(1:width_, 1:height_);
    XY_2D_ = [X(:) Y(:)];
    XY_2D_(~FILTERED_INDEX_, :) = [];
    XY_CENTERS_ =  PLX_CAT_' * XY_2D_./ (PR_' * ones(1, 2));

    for i = 1 : q_
        DIFF_            = (XY_2D_ - (ones(size(XY_2D_, 1), 1) * XY_CENTERS_(i, :))) .* (PLX_CAT_(:, i) * ones(1, 2));
        XY_COV_          = DIFF_' * (DIFF_ .* (PLX_CAT_(:, i) * ones(1, 2))) / PR_(i);
        [X, Y]           = meshgrid(1:width_, 1:height_);
        X                = X - XY_CENTERS_(i, 1);
        Y                = Y - XY_CENTERS_(i, 2);
        XY_2D_T_         = [X(:) Y(:)];
        PRIO(:, i)       = sqrt((det(XY_COV_))) * exp(- 0.5 * sum(XY_2D_T_ / XY_COV_ .* XY_2D_T_, 2));
    end
    
    if(USE_PRIO)
        PRIO = PRIO ./ repmat(max(sum(PRIO, 2), 10 ^ -24), [1 q_]);
    else
        PRIO = ones(height_ * width_, 1) * PR_;
    end
    
else
    PR_  = ones(1, q_) / q_;
    PRIO = ones(height_ * width_, 1) * PR_;
end

PRIO(isnan(PRIO)) = 1 / q_;