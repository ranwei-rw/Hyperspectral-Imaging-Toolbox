% Syntax:
%     K = reconstruct_specular_coefficient(HSZ)
% 
% Description:
%     Returns the specular coefficient stored on the HSZ struct
% 
% Input:
%     HSZ: NICTA pipeline structure
% 
% Output:
%     K: Specular coefficient
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.5
% Last Update Date: 29 Oct 2013


function K = reconstruct_specular_coefficient_(HSZ)

    K = HSZ.K.Factor/max(max(HSZ.K.Factor));
    
end