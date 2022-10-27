%   this function is designed to test whether a text string is ended by a
%   specified sub-string.
%
%   Description:
%       value = endswith(str, suffix, casesensitive)
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
% Version: 1.0.1
% Last Update Date: 20 Nov 2014
%
function value = endswith_(str, suffix, casesensitive)

    if ~exist('casesensitive', 'var')
        casesensitive = 0;
    end
    
    n = length(suffix);
    
    if length(str) < n || n == 0
        value =  false;
    else
        if casesensitive
            value = strcmp(str(end-n+1:end), suffix);
        else
            value = strcmpi(str(end-n+1:end), suffix);
        end
    end

end