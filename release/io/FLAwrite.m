% Write a hyperspectral FLA (flat) file to disk
%
% Syntax:
%       FLAwrite(filename, DATA, overwrite)
%
% Description:
%   This function is designed to output hyperspectral image data into a standard ENVI 
%   file (fla/hdr pair). Postfixes such as like RAW, IMG or DAT are also acceptable.
%
%   Inputs:
%       filename:  Name (including the path and extension) of the FLA file 
%                  being saved to disk
%       DATA:      Structure containing the image cube and the header. If
%                  DATA is a matrix only, header information will be deduced
%                  from it.
%       overwrite: whether to overwrite any existins files. value is
%       true/false, or 1/0
%
%  Note: the image data cube DATA.I is always assumed in BSQ format for
%        input.
%
% See also
%
%   FLAread, HSZwrite, HSZread, SLZwrite, SLZread
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2015 All Rights Reserved.
% Author: Ran Wei

function FLAwrite(filename, DATA, overwrite)

    if nargin > 2
        export_fla_(filename, DATA, overwrite);
    else
        export_fla_(filename, DATA);
    end

end
