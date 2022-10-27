% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Angle in degrees
function deg = angle_(L1, L2)
    deg = acos(abs(dot(L1, L2))/(norm(L1)*norm(L2)))/pi * 180;
end