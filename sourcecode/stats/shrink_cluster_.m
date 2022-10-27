% Syntax
%   [PLX, q, L_PRIO, CENTERS] = shrink_cluster(PLX, q, L_PRIO)
%       
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.2
% Last Update Date: 4 Aug 2014

function  [PLX, q, L_PRIO] = shrink_cluster_(PLX, q, L_PRIO)

    if size(PLX, 3) > 1
        PLX = reshape(PLX, [], q);
    end

    PLX_Q_ = zeros(1, q);

    for i = 1:q
        if mean(PLX(:, i)) < 0.05
            PLX_Q_(i) = 1;
        end
    end

    q = q - length(find(PLX_Q_));
    LGCQ_ = logical(PLX_Q_);
    L_PRIO(:, LGCQ_) = [];
    PLX(:, LGCQ_) = [];

    %   end of function: shrink_cluster_
end