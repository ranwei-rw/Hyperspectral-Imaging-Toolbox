%   function ae = angular_diff_(V1, V2)
%       gives the angular diffenence between given vectors V1 and V2. 
%   Input:
%       V1 and V2: input vectors
%
%   Outtput:
%       ad:  the error value, not in degree value
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei

function ad = angular_diff_(V1, V2)

	ad = acos( abs(dot(V1, V2))/(max(1e-4, norm(V1)*norm(V2))) );

end
