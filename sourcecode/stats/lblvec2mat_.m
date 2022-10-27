%   Syntax
%       Y = lblvec2mat(y, nc)
%
%   subroutine to convert 
%		1. label vector R_n in {0,1,..,K} to label matrix R_{nxK} (lblvec2mat)
%		2. label image hxw to label matrix image hxwxk
%       this is the inverse of routine lblmat2vec
%       make sure y is an integer array/matrix
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Antonio Robles-Kelly
% Version: 1.0.1
% Last Update Date: 15 Jan 2014

function Y = lblvec2mat_(y, nc)

    if isempty(y)
        Y = zeros(0, nc);
        return;
    end
    
    if ~exist('nc', 'var')
        nc = max(y(:));
    end
  
    id = eye(nc);
    idxnon0 = find(y);
    Y = zeros(prod(size(y)),nc);
    Y(idxnon0,:) = id(y(idxnon0),:);
    
    if ~isvector(y),
        Y = reshape(Y,size(y,1),size(y,2),nc);
    end

end