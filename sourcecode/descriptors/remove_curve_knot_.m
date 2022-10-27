function [t, NEW_KNOTS, NEW_CP_REFLEC, NEW_CP_WAVE] = ...
    remove_curve_knot_(degree, KNOTS, CP_REFLEC, CP_WAVE, index, mul, num, tolerance)

%function [t, NEW_KNOTS, NEW_CP_REFLEC, NEW_CP_WAVE] = 
%    remove_curve_knot(degree, KNOTS, CP_REFLEC, CP_WAVE, index, mul, num, tolerance)
%   This function is designed according to Algorithm A5.8 of the NURBS book (Page 185). It decides
%   whether a knot is removable and how many times. 
% 
%   Input:
%
%       degree:       degree of the curve.
%
%       KNOTS:     knot vector (row vector).
%
%       CP_REFLEC: control point coordinates in the reflectance dimension (stored as a height x
%                  width x (n + 1) 3D array while n + 1 is the number of control points per spectrum).
%
%       CP_WAVE:   control point coordinates in the wavelength dimension (stored as a (n + 1) x 1
%                  column vector, where n + 1 is the number of control points per spectrum). These
%                  control points are all the same for all the pixels. 
%
%       index:     index of the knot (in the knot vector) to be removed (index starts at 1).
%
%       mul:       the multiplicity of the knot to be removed.
%
%       num:       number of times intended to remove the knot.
%
%       tolerance: tolerance of the new control points displacement, in terms of distance.
%
%   Output:
%
%       t:             the number of times the indicated knot can be removed.
%
%       NEW_KNOTS:     the new knot vector after removing the knot t times.
%
%       NEW_CP_REFLEC: the new control point coordinates in the reflectance dimension after removing
%                      the knot t times, stored as a 3D array. 
%
%       NEW_CP_WAVE:   the new control point coordinates in the wavelength dimension after removing
%                      the knot t times, stored as a column vector. 
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
%
%   Version 1.0.1
%
%   Date: 2012.11.30

%   this function was adapted from RemoveCurveKnot3D.m

    %%%%%%%%%%%%%   code starts %%%%%%%%%%%%%%%%%%%%%%%%
    %   height, width: spatial dimensions of the image.
    [height, width, n] = size(CP_REFLEC);
    
    n = n - 1;          %   n: the number of control points - 1.
    m = n + degree + 1;    %   m: the number of knots - 1.
    
    ord   = degree + 1;    %   ord: number of coefficients
    fout  = floor((2 * index - mul- degree)/2);    %   First control point to be removed.
    last  = index - mul;
    first = index - degree;

    removed_knot = KNOTS(1, index); %   the knot to be removed.

    t = 0;
    %   don't allow remove more number of times than the multiplicity
    while t <= num-1 && t <= mul-1
        %   index difference (offset)
        offset = first-1; 
        TEMP_CP_REF(:, :, [1, last + 2 - offset]) = CP_REFLEC(:, :, [offset, last + 1]);
        TEMP_CP_WAVE([1, last + 2 - offset], 1) = CP_WAVE([offset, last + 1], 1);
        
        i  = first;
        j  = last;
        ii = 2; %   will have to access ii-1 later
        jj = last - offset + 1; 
        remflag = 0;

        while j - i > t
            %   Compute new control points for one removal step.
            alfi = (removed_knot - KNOTS(1, i))/(KNOTS(1, i + ord + t) - KNOTS(1, i));
            alfj = (removed_knot - KNOTS(1, j - t))/(KNOTS(1, j + ord) - KNOTS(1, j-t));
            
            TEMP_CP_REF(:, :, ii) = (CP_REFLEC(:, :, i) - (1 - alfi) * TEMP_CP_REF(:, :, ii-1))/alfi;
            TEMP_CP_REF(:, :, jj) = (CP_REFLEC(:, :, j) - alfj * TEMP_CP_REF(:, :, jj+1))/(1 - alfj);
            
            TEMP_CP_WAVE(ii, 1) = (CP_WAVE(i, 1) - (1 - alfi) * TEMP_CP_WAVE(ii-1, 1))/alfi;
            TEMP_CP_WAVE(jj, 1) = (CP_WAVE(j, 1) - alfj * TEMP_CP_WAVE(jj+1, 1))/(1 - alfj);
            
            i  = i + 1;
            ii = ii + 1;
            j  = j - 1;
            jj = jj- 1;
        end

        %   Check if the knot is removable.
        if j - i < t
            %   if the maximum difference in the new control point coordinates resulting from
            %   different linear equations 
            if (max(max(abs(TEMP_CP_REF(:, :, ii-1) - TEMP_CP_REF(:, :, jj+1)))) + ...
                abs(TEMP_CP_WAVE(ii-1, 1) - TEMP_CP_WAVE(jj+1, 1)) <= tolerance)
                remflag = 1;
            end
        else
            alfi = (removed_knot - KNOTS(1, i))/(KNOTS(1, i + ord + t) - KNOTS(1, i));
            
            if (max(max(abs(CP_REFLEC(:, :, i) - ...
                (alfi * TEMP_CP_REF(:, :, ii+t+1) + (1-alfi) * TEMP_CP_REF(:, :, ii-1))))) + ...
                abs(CP_WAVE(i, 1) - ...
                (alfi * TEMP_CP_WAVE(ii+t+1, 1) + (1-alfi) * TEMP_CP_WAVE(ii-1, 1))) <= tolerance)            
                remflag = 1;
            end
        end

        if remflag == 0
            %   Cannot remove any more knots. get out of for loop.
            break; 
        else 
            %   successful removal. Save new control points.
            i = first;
            j = last;
            while (j-i > t)
                CP_REFLEC(:, :, i) = TEMP_CP_REF(:, :, i+1 - offset);
                CP_REFLEC(:, :, j) = TEMP_CP_REF(:, :, j+1 - offset);   
                
                CP_WAVE(i, 1) = TEMP_CP_WAVE(i+1 - offset, 1);
                CP_WAVE(j, 1) = TEMP_CP_WAVE(j+1 - offset, 1);   
                
                i = i+1;
                j = j-1;
            end
        end

        first = first - 1;
        last  = last + 1;
        %   update t (the number of times the knot can be removed).
        t = t + 1; 
    end 

    if t == 0
        %   There is no change to the existing knots and control points.
        NEW_KNOTS     = KNOTS;
        NEW_CP_REFLEC = CP_REFLEC;
        NEW_CP_WAVE   = CP_WAVE;
                
        return;
    end    
    
    %   Store new knots and control points.
    NEW_KNOTS     = zeros(1, m + 1 - t); %   (t) knots will be removed
    NEW_CP_REFLEC = zeros(height, width, n+1-t);
    NEW_CP_WAVE   = zeros(n+1-t, 1);
    
    NEW_KNOTS(1, 1:index-t) = KNOTS(1, 1:index-t);
    NEW_KNOTS(1, [index+1:m+1]-t) = KNOTS(1, [index+1:m+1]);  %   Shift knots    

    %   Control points with index from i through to j will be overwritten.
    j = fout;
    i = j;
    for k = 1:t-1
        if mod(k, 2) == 1
            i = i + 1;
        else
            j = j -1;
        end
    end

    NEW_CP_REFLEC(:, :, 1:j-1)   = CP_REFLEC(:, :, 1:j-1);    
    NEW_CP_REFLEC(:, :, j:j+n-i) = CP_REFLEC(:, :, i+1:n+1);
    
    NEW_CP_WAVE(1:j-1, 1)   = CP_WAVE(1:j-1, 1);    
    NEW_CP_WAVE(j:j+n-i, 1) = CP_WAVE(i+1:n+1, 1);
end 

