% Syntax
%   PL_X = normalise_probability(PXL_, PRIO_L);
%       
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Antonio Robles-Kelly
% Version: 1.0.1
% Last Update Date: 9 Jan 2014

function PL_X = normalise_probability_(PXL_, PRIO_L)

    PL_X  = PXL_ .* PRIO_L;
    EVID_ = sum(PL_X, 2);
    
    if any(EVID_ == 0)
        ZR_          = find(EVID_ == 0);
        EVID_        = EVID_ + (EVID_ == 0);
        PL_X(ZR_, :) = PRIO_L(logical(ZR_), :);
    end

    PL_X = PL_X ./ repmat(EVID_, [1 size(PL_X, 2)]);
    
end
