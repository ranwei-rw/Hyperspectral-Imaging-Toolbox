%   this function is designed to test whether a text string begins with a
%   specified sub-string.
%
%   Description:
%       value = beginswith(str, suffix, casesensitive)
%  
%   Input: 
%   
%       str: the string to be tested
%       suffix: substr this function is looking for
%       casesensitive: case sensitive or not. default to 0 (no)
%
%   Output: 
%       value: true when str is ended by suffix or false, otherwise
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.0
% Last Update Date: 10 Feb 2015
%

function value = beginswith_(str, suffix, casesensitive)

    if ~exist('casesensitive', 'var')
        casesensitive = 0;
    end
    
    n = length(suffix);
    
    if length(str) < n || n == 0
        value =  false;
    else
        if casesensitive
            value = strcmp(str(1:n), suffix);
        else
            value = strcmpi(str(1:n), suffix);
        end
    end

end