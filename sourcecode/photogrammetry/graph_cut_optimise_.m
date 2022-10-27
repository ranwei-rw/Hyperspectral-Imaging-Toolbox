% Description
%
%   function [L_PAIR,  MAPPING, opt_lambda, opt_theta, opt_q] = graph_cut_optimise(I, LABELS, KPMASK, SPMATRIX, CNRPOINTS, pitch_height, pitch_width, THETA, LAMBDA, Q)
%   This function is to find the pair of labels in given Labels using
%   specified values of theta, lambda and Q. 
%
%   Input:
%       I:          Input image data cube
%       LABELS:     given label set
%       KPMASK:     the mask of pitch. 
%       SPMATRIX:   the specularity matrix of given image
%       CNRPOINTS:  the coordinates of left bottom corner points of pitches 
%       pitch_height and pitch_width: the dimension of pitches, both
%                   default to 20
%       THETA:      values of theta candidates, default to 0.001
%       LAMBDA:     values of lambda candidates, default to 0.7
%       Q:          values of q candidates, default to 1
%   Output:
%       L_PAIR:     selected pair of labels
%       MAPPING:    the mapping matrix for each pitch to labels
%       opt_lambda: optimal value of lambda
%       opt_q:      optimal value of q
%       opt_theta:  optimal value of theta
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2015 All Rights Reserved.
% Author: Ran Wei

function [L_PAIR, MAPPING, opt_lambda, opt_theta, opt_q] = graph_cut_optimise_(I, LABELS, KPMASK, SPMATRIX, CNRPOINTS, pitch_height, pitch_width, THETA, LAMBDA, Q)

    if ~exist('Q', 'var')
        Q = 1;
    end
    
    if ~exist('LAMBDA', 'var')
        LAMBDA = 0.7;
    end
    
    if ~exist('THETA', 'var')
        THETA = 1;
    end
    
    if ~exist('pitch_width', 'var')
        pitch_width = 20;
    end
    
    if ~exist('pitch_height', 'var')
        pitch_height = 20;
    end

    % Perform graph-cuts minimisation for pairs of current_label_pair
	optimal_energy = 1e10;
	% Generate combinations of pairs to run minimisation over
	%   pairs = combinations(LABELS); 
    % The following code is tested identical to the results from
    % function combinations
    %   start to replace function pairs = combinations(LABELS)
    lbsindx = 1:size(LABELS, 1);
	pairs = combvec(lbsindx, lbsindx)';
	% Remove pairs with the same index
    pairs = pairs(pairs(:, 1) ~= pairs(:, 2), :);

	% Now remove permutations, i.e. to have unique combinations
	pairs = unique(sort(pairs, 2), 'rows');
    %   end of replacing combinations
	
    if isempty(pairs)
        error('No valid label was found');
    end
    
    %   now pairs is row-format
    %   Labels are in n by bands matrix
    
    lambdasize  = length(LAMBDA);
    qsize       = length(Q);
    
    %   compose the 3D matrix for lambda and Q values, both lambda and Q
    %   are row vector at this moment
    [height, width, bands] = size(I);
    block_no_h = ceil(height/pitch_height);
    block_no_w = ceil(width/pitch_width);
    weights_s  = zeros(block_no_h, block_no_w);
    weights_p  = zeros(block_no_h, block_no_w);
    phi_s      = zeros(block_no_h, block_no_w, 2);
    phi_p      = zeros(block_no_h, block_no_w, 2);
    avgclr     = zeros(block_no_h, block_no_w, bands);
    
    r = 0;
    for i = 1:pitch_height:height
        r = r + 1;
        c = 0;
        for j = 1:pitch_width:width
            
            ht = i;                             %   height top
            hb = min(i+pitch_height-1, height); %   height bottom
            wl = j;                             %   width left
            wr = min(j+pitch_width-1, width);   %   width right
            c  = c+1;                           %   the x_th PATCH on x axle

            ph = hb - ht + 1;
            pw = wr - wl + 1;

            if ph < pitch_height
                %   reached the bottom of original image
                ht = height - pitch_height+1;
                hb = height;
                %ph = hb - ht + 1;
            end

            if pw < pitch_width
                %   reached right side of original image
                wl = width - pitch_width + 1;
                wr = width;
                %pw = wr - wl + 1;
            end
            %   re-crop image            
            
            PATCH = I(ht:hb, wl:wr, :);

            % 1. Parameters for phi_s
            % Get the local average colour of patches
            ii = sum(sum(PATCH));
            % Normalise the vector and store
            avgclr(r, c, :) = ii./norm(ii(:));
            weights_s(r, c) = sum(ii);
            if SPMATRIX(r, c) == true
                weights_p(r, c) = 1;
            end
        end
    end
    opt_q = 0;
    opt_lambda = 0;

    for p = 1:size(pairs, 1)
        %   Do optimisation for a pair of illuminant current_label_pair
        %   Using the pair information to compose a pair of current_label_pair which
        %   are picked out from all current_label_pair set. Then go through all
        %   possible pairs to get the one with minimum energy as the
        %   optimal pair of current_label_pair
        current_label_pair = [LABELS(pairs(p, 1), :); LABELS(pairs(p, 2), :)];		

        % 1. Compute the datacost term
        %   [Dc, avgclr] = unary(I, pitch_height, pitch_width, current_label_pair, SPMATRIX, KPMASK, li, qi);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   replacing content of unary 

        % Row and column indices of patches
        r = 0; 
        for i = 1:pitch_height:height
            r = r + 1;
            c = 0;
            for j = 1:pitch_width:width
                c = c + 1;
                %   continuing here
                % Compute phi_s	and phi_p for each label
                av = reshape(avgclr(r, c, :), [1, bands]);
                kpv = reshape(KPMASK(r, c, :), [1, bands]);
                for l = 1:size(current_label_pair, 1)
                    phi_s(r, c, l) = (1 - exp(-angular_diff_(av, current_label_pair(l, :))^2/(2*2.5^2) ))*weights_s(r, c);
                    phi_p(r, c, l) = (1 - exp(-angular_diff_(kpv, current_label_pair(l, :))^2/(2*2.5^2) ) )*weights_p(r, c);
                end
            end
        end

        % Compute datacost
        SCM = ones(2) - eye(2);
        [VCM, HCM] = pairwise_(I, avgclr, CNRPOINTS, pitch_height, pitch_width);
        for li = 1:lambdasize
            for qi = 1:qsize
                DATACOST = (1 - LAMBDA(li)).*phi_s.^(Q(qi)-1) + LAMBDA(li).*phi_p;

                %   end of unary.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 2. Compute the spatially varying smoothness terms

                % 3. Set up the spatially invariant smoothness term

                % 4. Compute graph-cut optimisation (Eqn 4)
                for k = 1:length(THETA)
                    theta = THETA(k);
                    gch = GraphCut('open', DATACOST, theta*SCM, VCM, HCM);
                    [gch, L_CAND] = GraphCut('expand', gch);
                    % Get the energy
                    [gch, energy] = GraphCut('energy', gch);
                    if energy < optimal_energy
                        MAPPING = L_CAND;
                        optimal_energy = energy;
                        opt_lambda = LAMBDA(li);
                        opt_q = Q(qi);
                        opt_theta = theta;
                        L_PAIR = current_label_pair;
                    end
                end
                GraphCut('close', gch);
            end
        end
    end

    MAPPING = MAPPING + 1;
end


% 
% function [vC, hC] = pairwise_(im, avgclr, cnrpoints, pH, pW)
% 	% Compute pairwise potential 
% 	[h, w, ~] = size(im);
% 	sz = size(avgclr);
% 	
% 	% Make vC and hC respectively
% 	vC = zeros(sz(1:2));
% 	hC = vC;
% 	for r = 1:sz(1)
% 		for c = 1:sz(2)
% 			% Note that we switch row and column indices since Matlab uses column major
% 			p1ix = sub2ind(sz(1:2), r, c);
% 			% Get the correponding row and column indices of the entire patch
% 			rix1 = cnrpoints(p1ix,1):min(cnrpoints(p1ix,1)+pH, h);
% 			cix1 = cnrpoints(p1ix,2):min(cnrpoints(p1ix,2)+pW, w);
% 
% 			if r > 1
% 				p2ix = sub2ind(sz(1:2), r-1, c);
% 				
% 				rix2 = cnrpoints(p2ix,1):min(cnrpoints(p2ix,1)+pH, h);
% 				cix2 = cnrpoints(p2ix,2):min(cnrpoints(p2ix,2)+pW, w);			
%                 n1 = sum(ismember(rix1+1, rix2));
%                 n2 = sum(ismember(cix1, cix2));
%                 n3 = sum(ismember(rix1, rix2));
%                 n4 = sum(ismember(cix1+1, cix2));
%                 h = n1*n2+n3*n4;
%                 if h ~= 0
%                     if avgclr(r, c) == avgclr(r-1, c)
%                         vC(r,c) = 0;
%                     else
%                         vC(r,c) = h;
%                     end
%                 end
% 
% 			end
% 			if c > 1
% 
% 				p2ix = sub2ind(sz(1:2), r, c-1);
% 				
% 				rix2 = cnrpoints(p2ix,1):min(cnrpoints(p2ix,1)+pH,h);
% 				cix2 = cnrpoints(p2ix,2):min(cnrpoints(p2ix,2)+pW,w);		
%                 n1 = sum(ismember(rix1+1, rix2));
%                 n2 = sum(ismember(cix1, cix2));
%                 n3 = sum(ismember(rix1, rix2));
%                 n4 = sum(ismember(cix1+1, cix2));
%                 h = n1*n2+n3*n4;
%                 if h ~= 0
%                     if avgclr(r, c) == avgclr(r, c-1)
%                         hC(r,c) = 0;
%                     else
%                         hC(r,c) = h;
%                     end
%                 end
% 
% 			end
% 		end
% 	end
% 
% end
function [vC, hC] = pairwise_(im, avgclr, cnrpoints, pH, pW)
	% Compute pairwise potential 
	[h, w, b] = size(im);
	sz = size(avgclr);
	
	% Make vC and hC respectively
	vC = zeros(sz(1:2));
	hC = vC;
	for r = 1:sz(1)
		for c = 1:sz(2)
			% Note that we switch row and column indices since Matlab uses column major
%			p1ix = sub2ind(sz(1:2), c, r);
			p1ix = sub2ind(sz(1:2), r, c);
			% Get the correponding row and column indices of the entire patch
			rix1 = cnrpoints(p1ix,1):min(cnrpoints(p1ix,1)+pH,h);
			cix1 = cnrpoints(p1ix,2):min(cnrpoints(p1ix,2)+pW,w);

			if r > 1
				p2ix = sub2ind(sz(1:2), r-1, c);
				
				rix2 = cnrpoints(p2ix,1):min(cnrpoints(p2ix,1)+pH,h);
				cix2 = cnrpoints(p2ix,2):min(cnrpoints(p2ix,2)+pW,w);			

				h = boundary_length(rix1,cix1, rix2, cix2);
				if h ~= 0
                    if avgclr(r, c) == avgclr(r, c-1)
                        hC(r,c) = 0;
                    else
                        hC(r,c) = h;
                    end
                end
			end
			if c > 1
				p2ix = sub2ind(sz(1:2), r, c-1);
				
				rix2 = cnrpoints(p2ix,1):min(cnrpoints(p2ix,1)+pH,h);
				cix2 = cnrpoints(p2ix,2):min(cnrpoints(p2ix,2)+pW,w);			

				h = boundary_length(rix1,cix1, rix2,cix2);
                if h ~= 0
                    if avgclr(r, c) == avgclr(r, c-1)
                        hC(r,c) = 0;
                    else
                        hC(r,c) = h;
                    end
                end
			end
		end
	end

end


function h = boundary_length(rix1, cix1, rix2, cix2)

	% Generate coordinates for the two patches
	p1ix = zeros(length(rix1)*length(cix1), 2);
	p2ix = zeros(length(rix2)*length(cix2), 2);

	p = 0;
	for r = 1:length(rix1)
		for c = 1:length(cix1)
			p = p+1;
			p1ix(p,:) = [rix1(r), cix1(c)];
		end
	end

	p = 0;
	for r = 1:length(rix2)
		for c = 1:length(cix2)
			p = p+1;
			p2ix(p,:) = [rix2(r), cix2(c)];
		end
	end

	% Run through first patch pixels to find if they have a shared boundary
	% Note that we need to test for all possible permutations
	h = 0;
	for p = 1:size(p1ix)
		temp_p = p1ix;
		% 1. Row neighbour
		if any(true == ismember( [temp_p(1)+1, temp_p(2)], p2ix, 'rows'))
			h = h+1;
		end
		% 2. Column neighbour
		if any(true == ismember( [temp_p(1), temp_p(2)+1], p2ix, 'rows'))
			h = h+1;
		end
	end

end