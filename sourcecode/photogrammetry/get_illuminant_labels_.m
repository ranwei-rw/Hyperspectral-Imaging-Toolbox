
function [LABELS, KPMASK, SPECULARITY, CORNORS] = get_illuminant_labels_(I, pitch_height, pitch_width, method, order)
	
    if ~exist('order', 'var')
        order = 1;
    end
    
    if ~exist('method', 'var')
        method = 'grey_edge';
    else
        method = lower(method);
        if ~strcmp(method, 'grey_edge') && ~strcmp(method, 'grey_world') && ...
           ~strcmp(method, 'maxrgb') && ~strcmp(method, 'shades_of_grey')
            error('Unknown colour constancy algorithm.');
        end
    end

	% Compute initial statistics and physics based estimates
	% [Ks, Kp, rootN, SPECULARITY, CORNORS] = initialestimates(I, pitch_height, pitch_width, method, order);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   initialestimates starts
    %   function [LIGHTMASK, KPMASK, rootN, SPECULARITY, CORNORS] = initialestimates(I, pitch_height, pitch_width, method, order)	

	[height, width, bands] = size(I);
    
	LIGHTMASK   = zeros(height, width, bands);  %   the statistics-based estimate is the LIGHTMASK
	KPMASK      = zeros(height, width, bands);
	CCMASK      = zeros(height, width, bands);  %   light corrected data for this image
	SPECULARITY = false(ceil(height/pitch_height), ceil(width/pitch_width));
    CORNORS     = zeros(ceil(height/pitch_height)*ceil(width/pitch_width), 2);  %   this a  vector storing the coordinates of the upperleft corner of PATCHes
    PATCH_no    = 0;	
    patch_count = 0;
    specu_count = 0;
    r           = 0; %   here r and c record the PATCH position number in blocks
    n           = 10;
    eror = 0;
    %   iterates through given image from left to right, top to bottom
	for i = 1:pitch_height:height
		r = r+1;                            %   the r_th block on vertical direction
        c = 0;
		for j = 1:pitch_width:width
            %   get the coordinates for current PATCH            
            ht = i;                                 %   height top
            hb = min(i+pitch_height-1, height);     %   height bottom
            wl = j;                                 %   width left
            wr = min(j+pitch_width-1, width);       %   width right
            c  = c+1;                               %   the x_th PATCH on x axle
            ph = hb - ht + 1;
            pw = wr - wl + 1;

            if ph < pitch_height
                %   reached the bottom of original image
                ht = height - pitch_height+1;
                hb = height;
                ph = hb - ht + 1;
            end

            if pw < pitch_width
                %   reached right side of original image
                wl = width - pitch_width + 1;
                wr = width;
                pw = wr - wl + 1;
            end
            %   re-crop image
                
            PATCH = I(ht:hb, wl:wr, :);
            framesize   = pitch_height*pitch_width;
            PATCH_no    = PATCH_no + 1;
            % Get the corner points of the PATCH
            
            CORNORS(PATCH_no, 1:2) = [i, j];
      
			%   get the illuminant estimation for this PATCH.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   replace [pLr, pLg, pLb, pOut] = rungeneral_cc(PATCH, method, order);

            switch method
                case 'grey_world'
                    [L, C] = estimate_lightsource_(PATCH, 0, 1, 0);
                case 'maxrgb'
                    [L, C] = estimate_lightsource_(PATCH, 0, -1, 0);
                case 'shades_of_grey'
                    [L, C] = estimate_lightsource_(PATCH, 0, 5, 0);
                otherwise% 'grey_edge'
                    [L, C] = estimate_lightsource_(PATCH, order, 5, 2);
            end					

            %   ends replacement of function rungeneral_cc(PATCH, method, order);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
			% Assign light source estimate to image mask, along with the
			% light corrected data for this PATCH 
            if isempty(L)
                error('Could not find estimations of illuminants');
            end
            
            CCMASK(ht:hb, wl:wr, :)	= C;
			PATCH_MASK = ones(ph, pw);
            eror = eror + 1;

            for k = 1:bands
                LIGHTMASK(ht:hb, wl:wr, k) = PATCH_MASK .* L(k);
            end
			
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   [SPMAP, sufspecular] = highlightpixels(PATCH);
            %   starts replacement of function highlightpixels
            %   function [specularMask, sufspecular] = highlightpixels(PATCH, t_b, t_s, t_sp)
            %   This function determines if a PATCH is sufficiently specular
            % We use default here as the value used in the paper, Beigpour 2014
            
            % Collect specular pixels
            P2D   = reshape(PATCH, [framesize, bands]);
            BR    = sum(P2D, 2)./bands;
            maxbr = max(BR(:));
            BR(BR < (0.2*maxbr)) = -1;
            MINP  = min(P2D, [], 2)./BR;
            SPMAP = zeros(framesize, bands);
            IDX   = MINP > (0.2*maxbr);
            SPMAP(IDX, :) = P2D(IDX, :);
            SPMAP = reshape(SPMAP, [ph, pw, bands]);
            
            % Determine if the PATCH is sufficiently specular
            sufspecular = false;    %   indicate how specular this patch is
            if sum(SPMAP(:))/bands >= n
                sufspecular = true;
            end
            %   ends replacement of function highlightpixels
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
			% Increase PATCH counts
			patch_count = patch_count + 1;
			if sufspecular
				specu_count = specu_count + 1;
                %   update light estimation using new patch data with
                %   specularity points removed
				                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   was lnew = IICestimate(SPMAP);
                %   replace function IICestimate(SPMAP); 
                %function [lr, lg, lb] = IICestimate(specpix, stepsize)
                stepsize = n;
                iic = zeros(ph, pw, bands);

                % 1. Project specular pixels to inverse-intensity chromaticity space
                %   this operation could be quite sparse, given the step
                %   size and the patch size
                for r = 1:stepsize:ph
                    for c = 1:stepsize:pw
                        chroma = sum(SPMAP(r, c, :)); 
                        if chroma > 0
                            iic(r, c, :) = SPMAP(r, c, :)./chroma;
                        else
                            continue;
                        end
                    end
                end

                % 2. Project clusters to Hough space and compute intersection
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   starts replacement of function houghspace
                %   the call was like the below
                %	lr = houghspace(iic(:,:,1), specpix);
                %   lg = houghspace(iic(:,:,2), specpix);
                %   lb = houghspace(iic(:,:,3), specpix);
                %   LNEW = houghspace(iic, SPMAP);
                %   function lc = houghspace(iic, SPMAP)
                LNEW = ones(bands, 1);
                % y = m*x + b
                b_start = 100;
                m_max   = 3000;
                b_end   = 990;
                params  = zeros(m_max, b_end-b_start+1);
                pa_idx  = zeros(bands*(b_end-b_start+1), 2);
                count   = 0;

                SPMAP   = sum(reshape(SPMAP, [framesize,  bands]), 3);
                SPMAP(SPMAP > 0) = SPMAP(SPMAP > 0)/10000;
                B       = b_start:b_end;
                
                for iicd = 1:bands
                    IIC_1D = reshape(iic(:, :, iicd), [framesize, 1]);

                    for p = 1:framesize

                        x = SPMAP(p);
                        y = IIC_1D(p);

                        if x < 0
                            continue;
                        end                    

                        M = abs(floor(2*(B/(1000*x) - y/x) + 0.5));
                        
%                         indx = M < m_max;
%                         K = M(indx)+1;
%                         params(K, indx) = params(K, indx) + 1;
%                         pos = params == n;
%                         cc   = sum(pos);
%                         if cc ~= 0
%                             pos = reshape(pos, [m_max*b_end, 1]);
%                             pos = pos_id(pos);
%                             [ki, kj] = ind2sub(size(params), pos);
%                             pa_idx(count+1, count+c) = [ki', kj'];
%                             count = count + cc;
%                         end

                        MB = M(M < m_max);
                        
                        for b = 1:length(MB)
                            k = MB(b) + 1;
                            params(k, b) = params(k, b)+1;
                            if params(k, b) == n
                                count = count+1;
                                pa_idx(count, :) = [k, b];
                            end
                        end
                        
%                         for b = 1:length(B)
%                             if M(b) < m_max
%                                 k = M(b) + 1;
%                                 params(k, b) = params(k, b)+1;
%                                 if params(k, b) == n
%                                     count = count+1;
%                                     pa_idx(count, :) = [k, b];
%                                 end
%                             end
%                         end
                        
                    end
                    
                    chroma = 0;
                    if count > 0
                        BHIST = zeros(1000, 1);
                        for p = 1:count
                            n1 = pa_idx(p, 1);
                            n2 = pa_idx(p, 2);
                            BHIST(n2) = BHIST(n2) + params(n1, n2);
                        end
                        [~, IB] = sort(BHIST, 'descend');
                        %   then get indexes for n biggest numbers in BHIST
                        if length(IB) > n
                            IB = IB(1:n);
                        else
                            n = length(IB);
                        end
                        
                        chroma = sum(IB);
                    end
                    
                    LNEW(iicd) = chroma/n/1000;
                    
                end
                %   ends replacement of function houghspace
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %   ends replacement of function IICestimate
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
				% Normalise result to PATCH and store result in KPMASK
                for kl = 1:bands
                    KPMASK(ht:hb, wl:wr, kl) = PATCH_MASK*LNEW(kl);
                end
				
			end
			% Store specular logical matrix
			SPECULARITY(r, c) = sufspecular;
		end
	end
	% Return number of cluster centres
	centers = [floor(sqrt(specu_count)), floor(sqrt(patch_count))];
    %   initialestimates ends
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Construct labels
    %   from here, Ks should be LIGHTMASK and Kp should be KPMASK

	% The number of cluster centres is rootN
	cdata = reshape(LIGHTMASK, [height*width, bands]);
	[~, ligthmap_cnrs] = kmeans(cdata, centers(1), 'Distance', 'cosine', 'EmptyAction', 'drop', 'MaxIter', 50);

    % Remove zero points
    cdata = reshape(KPMASK, [height*width, bands]);
    cdata = cdata(sum(cdata, 2) >= 0.001, :);

    % In the case of zero illuminant(s), then cdata is empty
    kp_cnrs = [];
    if centers(2) > 0.01 && ~isempty(cdata)
        [~, kp_cnrs] = kmeans(cdata, centers(2), 'Distance', 'cosine', 'EmptyAction', 'drop', 'MaxIter', 50);
    end
    
    % Get the ambient illuminant estimate, L_0
	LE = estimate_lightsource_(I, method, order);
    
 	% Mean values of the K clusters become the labels?
	% Combine illuminant clusters
    %   here K will be a n by bands matrix
    
    if ~isempty(kp_cnrs)
        K = [ligthmap_cnrs; kp_cnrs];
        K = K(~isnan(sum(K, 2)), :);
        LABELS = [K; LE];
    else
        LABELS = LE;
    end
    
	% Merge labels which differ by less than 0.5 degrees
	%   LABELS = mergelabels(LABELS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   start to replace LABELS = mergelabels(LABELS);
    % Generate combinations to test for angular error
    %   pairs = combinations(LABELS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   start to replace function pairs = combinations(LABELS)
    lbsindx = 1:size(LABELS, 1);
	pairs = combvec(lbsindx, lbsindx)';
	% Remove pairs with the same index
    pairs = pairs(pairs(:, 1) ~= pairs(:, 2), :);

	% Now remove permutations, i.e. to have unique combinations
	pairs = unique(sort(pairs, 2), 'rows');
    %   end of replacing combinations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lbsindx = true(1, size(LABELS, 1));
    threshold = 0.5*pi/180;
    for p = 1:size(pairs, 1)
        V1 = LABELS(pairs(p, 1),:);
        V2 = LABELS(pairs(p, 2),:);

        if acos( abs(dot(V1, V2))/(max(1e-4, norm(V1)*norm(V2))) ) < threshold
            % Remove the second index (which should be larger than the first)
            lbsindx(pairs(p, 2)) = false;
        end
    end
    % Update label set
    LABELS = LABELS(lbsindx, :);

    %   ends replacing LABELS = mergelabels(LABELS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Remove saturated illuminants (Section IV-D)
	
    [ln, lb] = size(LABELS);
    V2 = 1/sqrt(lb)*ones(1, lb);
	threshold = 15*pi/180;
    if size(LABELS, 2) > 3
        %   double the threshold for hyperspectreal images
        threshold = threshold*2;
    end
	% Exclude any illuminants that exceed this limit
	lbsindx = true(1, ln);
	for p = 1:ln
        V1 = LABELS(p, :);
		if angular_diff_(V1, V2) >= threshold
			lbsindx(p) = false;
		end
	end
	% Update labels
	LABELS = LABELS(lbsindx, :);

end