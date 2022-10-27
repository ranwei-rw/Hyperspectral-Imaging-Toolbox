%   this is a function used to do meanshift clustering on a given set of
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
% Version: 1.0.6
% Last Update Date: 22 Oct 2015
%

function [MS_MODES, MAPPING] = meanshift_(MODES, WEIGHTS)

    if ~exist('MODES', 'var')
        error('error in input arguments');
    end
    
    %   get the size of data set
    [band, num] = size(MODES);
    
    if ~exist('WEIGHTS', 'var')
        WEIGHTS = ones(num, 1);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %   Step 0. Rescale input modes to range [0, 65535]
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %min_m = min(MODES(:));
    %max_m = max(MODES(:));
    %MODES = (MODES - min_m)/(max_m - min_m)*65535;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %   Step 1. Find initial window and weightdp2 for each of the neurons
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   
    %   this part is from the computePilotPoint
    %   the underlying algorithm is: for any particular mode, to get its
    %   current adaptive window, follow steps below
    %       1. Create a buffer of distance BD and initialised to 0
    %       2. Compute the distance from this mode to all other mode.
    %          The distance will be treated as integers.
    %       3. If the distance d > d_threshold, discard it because it's too
    %          far away (this d_threshold is also the size of BD). If the
    %          distance is less than d_threshold, then BD[d] will be
    %          increassed by 1.
    %       4. Go through all neurons and compute all distances and save the
    %          count in BD
    %       5. Starting from BD[1], count the total number n of neurons that
    %          are within distance i to current mode. 
    %               n = n+BD[i];
    %       6. When the value of n reaches a threshold, say 20, then, teh
    %          search stops and the value of i becomes the new window for
    %          this mode. If n could not reache the threshold after go
    %          through all distance buffer BD, then the new window size
    %          will reach maximum value (d_threshold).
    %       7. Repeat step 3 to 6 for all neurons.
    %   this step was based on void
    %   FAMS::ComputePilotPoint::operator()(const tbb::blocked_range<int>
    %   &r) (line 99, mfams.cpp)
    
    WINDOW    = zeros([1, num]);
    WEIGHTDP2 = zeros([1, num]);
    
    fams_float_shift   = 100000;
    fams_alpha         = 1;
    distance_threshold = 400;
    
    wjd    = 300;
    thresh = 20;
    weight_threshold = 50;
    
    for i = 1:num
        diff = zeros(size(MODES));
        %   go through every point
        temp = zeros(distance_threshold, 1);
        P    = MODES(:, i);
        
        for j = 1:band
            diff(j, :) = MODES(j, :) - P(j);
        end
        
        %   D T and n are defined to to get L1 norm
        D = ceil(sum(abs(diff))/wjd);
        
        [M, N] = hist(D, 0:distance_threshold+1);
        M = M(2:distance_threshold+1);
        N = N(2:distance_threshold+1);
        temp(N) = temp(N) + M';

        indicator = 0;
        acc = 0;
        for k = 1:distance_threshold
            acc = acc + temp(k);
            if acc > thresh
                indicator = k;
                break;
            end
        end
        
        if indicator == 0
            indicator = distance_threshold;
        end
        
        WINDOW(i) = indicator*wjd;
        WEIGHTDP2(i) = power(fams_float_shift/WINDOW(i), (band+2)*fams_alpha);
    end 
    %   What the above section does is to compute the cfams.PrepareFAMS(NULL);
    %   To get all neurons' window and weightdp2 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %   Step 2. Compute converging modes for each neuron
    %
    %   In this step, we'll go through every mode(neuron) from its
    %   beginning location to its final convergency location. In each
    %   iteration of convergence, we'll compute its new position and new
    %   updating radius. We'll stop when the new updating radius equals to
    %   zero, or the maximum iteration number is reached.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  Now the updating window and weights are available. 
    %   Suppose that we have 1000 modes while each of them is associated
    %   with a weight (WEIGHTS)
    max_iteration_num = 100;

    %   define the variable holding final converging position and window
    tmode.MODES  = zeros(size(MODES));
    tmode.WINDOW = zeros([num, 1]);
    for i = 1:num
        %   get the original location of i-th mode
        CUR_P = MODES(:, i);
        OLD_P = zeros(size(CUR_P));
        cur_radius = WINDOW(i); %   get the initial updating radius

        for j = 1:max_iteration_num
            
            if sum(abs(OLD_P - CUR_P)) < 1e-09
                %   stops when mode location stops changing
                break;
            end
            %   backup current position
            OLD_P = CUR_P;
            
            %%%%%%%%
            %   compute new position and new radius
            %   this part was the content of function
            %   fams::DoMSAdaptiveIteration (line 273 mfams.cpp)
            %%%%%%%%
            
            D = zeros(size(MODES));
            
            for k = 1:length(CUR_P)
                D(k, :) = abs(MODES(k, :) - CUR_P(k));        
            end
            %   get the L1 distance
            L1D = sum(D);
            L1D_IND = find(L1D < cur_radius);
            X = ones(size(L1D_IND));
            X = X - L1D(L1D_IND)./WINDOW(L1D_IND);
            W = WEIGHTDP2(L1D_IND) .* (X.*X);
            total_weight = sum(W);
            SEL_M = MODES(:, L1D_IND);

            for n = 1:length(CUR_P)
                CUR_P(n) = sum(SEL_M(n, :) .* W, 2);
            end

            if total_weight == 0
                %   stop when no more radius;
                cur_radius = 0;
                CUR_P = OLD_P;
            else
                [~, d_i]   = min(L1D(L1D_IND));
                cur_radius = WINDOW(L1D_IND(d_i));
                CUR_P      = CUR_P./total_weight;
            end
            %   end of fams::DoMSAdaptiveIteration
        end
        
        %   save the location
        tmode.MODES(:, i) = CUR_P;
        tmode.WINDOW(i)   = cur_radius;
    end
    
    
    %   now tmode is holding all the modes that initial points converge to.
    %   It's tiem to classify or prune these modes because moany of them
    %   are identical or very close. This was caused by the fact that
    %   multiple initial points converge to the same final position.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %   Step 3. prune modes
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   Define several local variables for mode pruning
    mode_num_threshold = 50; %   this is the max number of modes should reach during pruning,
    invalid_mode       = zeros([num, 1]);           %   whether each mode is invalid or valid
    mode_size          = zeros([num, 1]);           %   how many hits this each mode has
    mode_weight        = zeros([num, 1]);           %   the weight for each modes. It gets initial values from WEIGHTS but not necessarily equal to WEIGHTS later
    total_mode_num     = 1;                         %   count how many modes are generated
    mode_weight(1)     = WEIGHTS(1);                %   initialise values for the first mode
    mode_size(1)       = 1;
    updated_modes      = tmode.MODES;               %   create a buffer for updated modes
    cmodes             = zeros(size(tmode.MODES));
    cmodes(:, 1)       = tmode.MODES(:, 1);
    LNK_INDEX          = zeros([num, 1]);           %   This is the list recording initial mode[i] is linked to which mode in cmode (final modes);
    LNK_INDEX(1)       = 1;
    %   in each circle of the following loop, it looks for a mode which is
    %   of the closest weighted distance to mode(i) and record its index
    for i = 2:num
        %   define a variable for distance
        min_distance  = realmax;    %   this equal to realmax('double')
        min_dis_index = -1;         %   the index of the mode which has current minimum distance to the candidate mode
        %   compute the weight distance from this 
        weighted_dis = 0;
        
        for j = 1:total_mode_num
            if invalid_mode(j)
                %   when this mode is tagged invalid, jump to next one
                continue;
            end
            weighted_dis = sum(abs(cmodes(:, j)/mode_size(j) - tmode.MODES(:, i)));
            if weighted_dis < min_distance
                min_distance = weighted_dis;
                min_dis_index = j;
                if weighted_dis <= 1e-09
                    %   means current mode is already with the shortest
                    %   distance with MODES(:, i); There is no need to do the
                    %   rest comparing
                    break;
                end
            end
        end
        
        %   here we get the minimum distance from mode min_dis_index to ith
        %   mode, then we see whether it's in the radius of its windows
        if min_distance < tmode.WINDOW(i)/2
            %   in this case this mode is with the update radius. 
            cmodes(:, min_dis_index) = cmodes(:, min_dis_index) + tmode.MODES(:, i);
            updated_modes(:, i) = cmodes(:, min_dis_index);
            %   record this mode
            mode_size(min_dis_index) = mode_size(min_dis_index) + 1;
            mode_weight(min_dis_index) = mode_weight(min_dis_index) + WEIGHTS(i);
            
            %   normalise (to its weight) this modes in updated_modes
            updated_modes(:, i) = updated_modes(:, i)/mode_size(min_dis_index);
            LNK_INDEX(i) = min_dis_index;
        else
            %   even closest mode is out of radius of current modes. need
            %   to create a new mode
            total_mode_num = total_mode_num + 1;
            cmodes(:, total_mode_num) = tmode.MODES(:, i);
            mode_size(total_mode_num) = 1;
            mode_weight(total_mode_num) = WEIGHTS(i);
            LNK_INDEX(i) = total_mode_num;
        end
        
        %   check the number of modes to prevent too many modes.
        
        if total_mode_num > mode_num_threshold
            %   check for modes which have fewer members. by setting the
            %   invalidation array to 1
            for k = 1:total_mode_num
                if mode_weight(k) < 3
                    invalid_mode(k) = 1;
                end
            end
        end
    end
    %   This corresponds to Line 562 in mfams.cpp. Now all initial modes
    %   are mapped into modes in cmode.
    
    %   now to sort the mode weights
    mode_weight = mode_weight(1:total_mode_num);
    mode_size   = mode_size(1:total_mode_num);
    [mode_weight, mode_weight_sorted_index] = sort(mode_weight,  'descend');
    
    %   resize this sorted modes into a shorter version, allows only modes
    %   which has weigth larger than a threshold

    n = 2;
    for i = n:total_mode_num
        %   here we don't need to worry about invalid modes. Because the
        %   threshold for invalid modes is weight < 3 which is much smaller
        %   than the weight_threshold        
        if mode_weight(i) < weight_threshold
            n = i-1;
            break;
        else
            n = i;
        end
            
    end
    %   This line corresponds code above line 587 in mfams.cpp
    
    
    %   now n contains the number of modes which has weight larger than
    %   given size
    if n < total_mode_num
        mode_weight = mode_weight(1:n);
        mode_weight_sorted_index = mode_weight_sorted_index(1:n);
        mode_size = mode_size(mode_weight_sorted_index);
        cmodes = cmodes(:, mode_weight_sorted_index);
        %   update current max mode number
        total_mode_num = n;
        %   line 607 in mfams.cpp
    end
    
    %   From now on, we have the short-listed final modes which are saved
    %   in cmodes. Then, it's time to test them againest initial modes
    %   start from final mode 1, find the closest mode for it
    min_distance  = realmax;
    min_dis_index = -1;
    init_mode_index = zeros([num, 1]);
    for k = 1:total_mode_num
        %   get the distance
        weighted_dis = sum(abs(cmodes(:, k)/mode_size(k) - tmode.MODES(:, 1)));
        if weighted_dis < min_distance
            min_distance = weighted_dis;
            min_dis_index = k;
        end
    end
    
    if min_dis_index ~= -1
        init_mode_index(1) = min_dis_index;
    else
        init_mode_index(1) = 0;
    end
    
    %   process the rest of initial modes
    for i = 2:num
        min_distance  = realmax;
        min_dis_index = -1;
        for j = 1:total_mode_num
            weighted_dis = 0;
            %   compute the L1 distance
            weighted_dis = sum(abs(cmodes(:, j)/mode_size(j) - tmode.MODES(:, i)));
            if weighted_dis < min_distance
                min_distance = weighted_dis;
                min_dis_index = j;
                if weighted_dis <= 1e-09
                    %   means current mode is already with the shortest
                    %   distance with MODES(:, i); There is no need to do the
                    %   rest comparing
                    break;
                end
            end
        end
        
        %   get the distance, check the radius
        if min_dis_index >= 1
            cmodes(:, min_dis_index) = cmodes(:, min_dis_index) + tmode.MODES(:, i);
            updated_modes(:, i) = cmodes(:, min_dis_index);
            init_mode_index(i) = min_dis_index;
            mode_size(min_dis_index) = mode_size(min_dis_index) + 1;
            if mode_size(min_dis_index) > 1
                updated_modes(:,i) = updated_modes(:,i) / mode_size(min_dis_index);
            end
        else
            init_mode_index(i) = 10001;
        end
        %   end of content of PruneModes at Lin 677 in mfams.cpp
    end
    %   now, updated_modes contains all the final mode positions
    %   corresponding to initial modes. We need to get these averaged
    %   total_mode_num position using this information.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %   Step 4. Link input with final modes
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   actually, upon last line, we've got all the information we need to
    %   index all input modes, which is saved in variable init_mode_index.
    %   It has num (1000 for now) elements and each of them
    %   init_mode_index(i) stors the final mode that the i-th input modes
    %   is pointing to.
    MAPPING = init_mode_index;
     
    
    MS_MODES = zeros([band, total_mode_num]);
    ave_count = zeros([total_mode_num, 1]);
    for i = 1:num
        pos = init_mode_index(i);
        if pos ~= 0 && pos ~= 10001
            MS_MODES(:, pos) = MS_MODES(:, pos) + updated_modes(:, i);
            ave_count(pos) = ave_count(pos) + 1;
        end
    end
    
    for i = 1:total_mode_num
        MS_MODES(:, i) = MS_MODES(:, i)/ave_count(i);
    end

%   end of function
end