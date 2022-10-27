%   this is a function used to do mean shift clustering on a given set of
%   data. For example, a set of neurons for neural network applications.
%   The task is to cluster these data into a set of clusters.
%
%   Description:
%       [MS_MODES, MAPPING] = meanshift(MODES, WEIGHTS)
%  
%   Input: 
%   
%       MODES: mode set to be clusterred. Each element is a vector while the
%              number of elements is determined by the column number (size:
%              depth by mode_number
%       WEIGHTS: the number of pixels that each mode represents
%
%   Output: 
%       MS_MODES:  the final mode matrix after meanshift clustering (size:
%                  depth by final_mode_num
%       MAPPING: the mapping between input mode matrix MODES and output MS_MODES
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2015 All Rights Reserved.
% Author: Ran Wei

function [MS_MODES, MAPPING] = meanshift(MODES, WEIGHTS)
    
    switch nargin 
        case 2
            [MS_MODES, MAPPING] = meanshift_(MODES, WEIGHTS);
        case 1
            [MS_MODES, MAPPING] = meanshift_(MODES);
        otherwise
            error('error in input arguments');
    end

end