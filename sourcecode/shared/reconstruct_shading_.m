% Syntax:
%   G = reconstruct_shading(HSZ);
% 
% Description:
%   Returns the shading as stored on the HSZ structure
% 
% Input:
%   HSZ: NICTA pipeline structure
% 
% Output:
%   G: Shading factor
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.5
% Last Update Date: 29 Oct 2013


function G = reconstruct_shading_(HSZ)

    G = HSZ.S.Factor/max(max(HSZ.S.Factor));
    
end