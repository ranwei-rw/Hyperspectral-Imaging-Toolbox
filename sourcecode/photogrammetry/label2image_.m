% Description 
%   function I = label2image(S, CNRPOINTS, MAPPING, LABELS, ph, pw)
%   
% This function is to assign label to pitches of an image using recovered
% illuminants
%   
% Input: 
%   S:          the size information of original image, a vector contains
%               height, width and bands dimensions
%   CNRPOINTS:  coordinates of top left corners of pitches for labelling
%   MAPPING:    mapping matrix indicates which label a pitch belongs to
%   LABELS:     illuminants labels.
%   ph, pw:     height and width of pitches, by default they are 20.
%
% Output:
%   I:          image filled with label illuminants, which is of size S
%
% This computer code is subject to copyright: (c) National ICT Australia
% Limited (NICTA) 2015 All Rights Reserved. 
% Author: Ran Wei

function I = label2image_(S, CNRPOINTS, MAPPING, LABELS, ph, pw)

    if ~exist('pw', 'var')
        pw = 20;
    end
    
    if ~exist('ph', 'var')
        ph = 20;
    end
    
	% Project labels for patch nodes into illuminant estimate for entire image
    if length(S) < 2
        error('Error in image size S');
    elseif length(S) < 3
        S(3) = 1;
    end
    
	I = zeros(S);
    ind = 1;
    [blockheight, blockwidth] = size(MAPPING);
    
    if size(LABELS, 2) < S(3)
        bands = size(LABELS, 2);
    else
        bands = S(3);
    end
    
	for i = 1:blockheight
		for j = 1:blockwidth
            %   go through pitches defined by CNRPOINTS
            %   get the coordinates for each pitch
            cnr = CNRPOINTS(ind, :);
            ht = cnr(1);
            wl = cnr(2);
            hb = min(ht + ph - 1, S(1));
            wr = min(wl + pw - 1, S(2));
            pitch_label = LABELS(MAPPING(i, j), :);
            for k = 1:bands
                I(ht:hb, wl:wr, k) = pitch_label(k);
            end
            ind = ind + 1;
		end
	end

end



