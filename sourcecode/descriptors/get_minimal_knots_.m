%     function [MIN_KNOTS, MIN_CP_REFLEC, MIN_CP_WAVE] = ...
%            get_minimal_knots(I, WAVELENGTHS, degree, KNOTS, CP_REFLEC, CP_WAVE, alpha, target_knotnum)
%
%   This function is the implementation of the knot removal algorithm for a hyperspectral image.
%   It's coded according to the `Algorithm 1' in the paper `A B-Spline Representation of Spectral
%   Images' by Cong Phuoc Huynh and Antonio Robles-Kelly.
%
%   Input:
%
%       I:              the input set reflectance image - a 3D array of size (height x width x bands).
%
%       WAVELENGTHS:    the wavelengths of the multi-spectral image.
%
%       degree:         the degree of the basis functions.
%
%       KNOTS:          the original knot vector (in the form of a row vector).
%
%       CP_REFLEC:      the initial coordinates of control points in the reflectance dimension,
%                       where P_ref(h, w, :) are the control point coordinates across all the bands
%                       at pixel coordinate (h, w). It is of size height x width x bands. Using
%                       function global_interpolation.m can generate such control points.
%
%       CP_WAVE:        a 1 x bands row vector - the initial control point coordinates in the
%                       wavelength dimension. These coordinates in this dimension is the same for
%                       all the pixels in the given image. Using global_interpolation.m can generate
%                       such control points.
%
%       alpha:          balance factor between 0 and 1.
%
%       target_knotnum: the target number of KNOTS remained.
%
%       tolerance:      tolerance on the control point displacement. If not set, it will be set to INF to allow as many knots to be removed as specified.
%
%    Output:
%
%       MIN_KNOTS:      the minimal knot vector, in the form of a row vector.
%
%       MIN_CP_REFLEC:  the minimal control point coordinates in the reflectance dimension - in the
%                       form of a 3D array.
%
%       MIN_CP_WAVE:    the minimal control point coordinates in the wavelength dimension - in the
%                       form of a row vector.
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 - 2014 All Rights Reserved.
% Author: Ran Wei
% Version: 1.0.7
% Last Update Date: 9 June 2015
%   improved routine. 

% Version: 1.0.6
% Last Update Date: 29 Aug 2014


% Version: 1.0.5
% Last Update Date: 29 Oct 2013
%   version 1.0.2 is based on MinimalKnots3D_v2.m instead of v3 in 1.0.1

%   Version 1.0.1
%
%   Date: 2012.11.15
%   this function is modified from Cong's MinimalKnots3D_v3.m
%

%   From now on, the code is modified from MinimalKnots3D_v2. Mainly only the argument names are changed.

function [MIN_KNOTS, MIN_CP_REFLEC, MIN_CP_WAVE] = ...
    get_minimal_knots_(I, WAVELENGTHS, degree, KNOTS, CP_REFLEC, CP_WAVE, alpha, target_knotnum, tolerance, debug)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Code starts:
    if ~exist('debug', 'var')
        debug = 1;
    end

    if ~exist('tolerance', 'var')
        tolerance = Inf;
    end

    %   tolerance of error in the wavelength while searching for a parametric point t corresponding to a wavelength.
    tolerance_wave = 1e-4;

    %   get the image dimensions and pixel number for each band
    [height, width, ~] = size(I);

    %   reduction in cost -- initially the curve fits every data point.
    previous_sse        = 0;
    candidate_sse       = 0;
    REMAINING_KNOTS     = KNOTS;
    REMAINING_CP_REFLEC = CP_REFLEC;
    REMAINING_CP_WAVE   = CP_WAVE;
    TMP_REM_KNOTS       = KNOTS;
    TMP_REM_CP_REFLEC   = CP_REFLEC;
    TMP_REM_CP_WAVE     = CP_WAVE;

    %   initialise the iteration counter to 0
    iter = 0;
    while (true)
        iter = iter + 1;

        %   get knots number that are currently existing
        remaining_knots_num = length(REMAINING_KNOTS);

        %   when the remaining knots is less or equal to required knot number, break.
        if remaining_knots_num <= target_knotnum
            break;
        end

        % initialse the variable which will record the max reduction in current iteration
        max_reduction = -1;
        % initialse the variable which will record the index of the knot to be removed
        knot_index_remove = -1;

        for i = degree + 3 : remaining_knots_num - degree - 2
            %   only remove internal KNOTS (excluding the first knot other than zero and the last
            %   knot other than one). Try removing the knot index i resulting from the previous
            %   iteration.
            [rem_time, NEW_KNOTS, NEW_CP_REFLEC, NEW_CP_WAVE] ...
                = remove_curve_knot_(degree, REMAINING_KNOTS, REMAINING_CP_REFLEC, REMAINING_CP_WAVE, i, 1, 1, tolerance);

            if rem_time > 0 %  if the knot is removable.
                tnew        = find_parapoints_(WAVELENGTHS, tolerance_wave, degree, NEW_KNOTS, NEW_CP_WAVE);
                [TEMP_R, ~] = reconstruct_curve_(NEW_KNOTS, tnew, NEW_CP_REFLEC, NEW_CP_WAVE, degree);
                DIFF_R      = abs(TEMP_R - I);
                new_sse     = sum(DIFF_R(:))/(height * width);
                reduction   = (1 - alpha) + alpha * (previous_sse - new_sse);
                
                if max_reduction < reduction
                    max_reduction     = reduction;
                    candidate_sse     = new_sse;
                    knot_index_remove = i;
                    TMP_REM_KNOTS     = NEW_KNOTS;
                    TMP_REM_CP_REFLEC = NEW_CP_REFLEC;
                    TMP_REM_CP_WAVE   = NEW_CP_WAVE;
                end

            end
        end

        if max_reduction <= 0
            break;
        end

        %   if a candidate knot is selected.
        if knot_index_remove >= degree + 3 &&  knot_index_remove <= remaining_knots_num - degree - 2
            % Store the SSE for the candidate knot, which will also be the old
            %   SSE for the next iteration.
            % Note: One can compute the cost after knot removal directly.
            % However, in cases where there are a lot of bands, we can
            % instead compute the cost reduction by making use of the local
            % support property around the neighbouring bands of the removed
            % one.
            previous_sse = candidate_sse;

            REMAINING_KNOTS = TMP_REM_KNOTS;
            REMAINING_CP_REFLEC = TMP_REM_CP_REFLEC;
            REMAINING_CP_WAVE = TMP_REM_CP_WAVE;
            if debug >= 2
                s = sprintf('Knot %d is removed', knot_index_remove);
                disp(s);
            end
        else
            if debug >= 2
                disp('Scyllarus:descriptors:get_minimal_knots - No removeable knot was detected');
            end
            break;
        end

    end

    MIN_KNOTS     = REMAINING_KNOTS;
    MIN_CP_REFLEC = REMAINING_CP_REFLEC;
    MIN_CP_WAVE   = REMAINING_CP_WAVE;

end