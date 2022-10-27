% Compute the NURBS representation for spectral data
%
% Syntax
%   [KNOTS, CP_REF, CP_WAVE] = get_nurbs(R, WAVES, degree, knot_thresh, alpha, iter)
%
% Description
%
%       This function is designed to a generate NURBS representation for a 
%       given set of hyperspectral data.
%
% Input:
%
%       R:           Input spectral data, stored as a 3D array of size (height x width x bands).
%       WAVES:       Wavelength vector for the spectra.
%       degree:      Degree of the basis (polynominal) functions. By default degree = 3;
%       alpha:       A balancing factor between 0 and 1, which is the weight of the data 
%                    closeness term (reconstruction error) 
%                    with respect to the original data. By default alpha = 0.1
%       knot_thresh: The threshold for the number of minimal knots. By default 
%                    knot_thresh = band - 2;
%       iter:        the maximum number of iterations. By default iter = 10
%       debug:       Whether to show debugging information. 1 or 2 for yes and 0 for no (default)
%   
% Output:
%
%       KNOTS:   Final minimal knot vector - row column.
%       CP_REF:  Minimal control point coordinates in the reflectance
%                dimension - in the form of a 3D array of size (height x width x ctrl_pts).
%       CP_WAVE: Minimal control point coordinates in the wavelength
%                dimension - in the form of a row vector (1 x ctrl_pts).
%
% Example
%
%   [KNOTS, CP_REF, CP_WAVE] = get_nurbs(HSZ.S.Elements, HSZ.HDR.wavelength, 2);
%            
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version 1.0.6
% Last Update Date: 28 Aug 2014

% Version: 1.0.5
% Last Update Date: 29 Oct 2013

function [KNOTS, CP_REF, CP_WAVE] = get_nurbs_(R, WAVES, degree, knot_thresh, alpha, iter, debug)

    if ~exist('iter', 'var')
        iter = 10;
    end

    if ~exist('alpha', 'var')
        alpha = 0.5;
    end

    if ~exist('knot_thresh', 'var')
        n = size(R, 3);
        knot_thresh = n - 2;
    end

    if ~exist('degree', 'var')
        degree = 2;
    end

    if ~exist('debug', 'var')
        debug = 0;
    end

    iter_num    = 0;
    r_sample    = R; %   the resampled reflectance in each iteration.
    wave_sample = WAVES; %   the resampled wavelengths in each iteration.

    while (true)
        iter_num = iter_num + 1;
            
        if debug >= 1
            s = sprintf('Iteration %d ...\n', iter_num);
            disp(s);
        end
        
        %   Initial global interpolation
        [knots0, ~, cp_ref_0, cp_wave_0] = global_interpolation_(r_sample, wave_sample, degree);

        %   Remove knots
        [KNOTS, CP_REF, CP_WAVE] = ...
            get_minimal_knots_(r_sample, wave_sample, degree, knots0, cp_ref_0, cp_wave_0, alpha, knot_thresh);

        %allRemovedKnots = [allRemovedKnots; removed_knots];

        init_knot_num = size(knots0, 2);  %   initial number of knots.
        curr_knot_num = size(KNOTS, 2);   %   the minimal number of knots after removal

        %   if no further knots are removed, we can stop here.
        if curr_knot_num == init_knot_num
            break;
        end

        %   Check to see if we can break here without further computation.
        if curr_knot_num <= knot_thresh || iter_num >= iter
            break;
        end

        params(1, 1) = KNOTS(1, 1 + degree);
        %           params(1, curr_knot_num - 2*degree + 1) = KNOTS(1, curr_knot_num - degree);
        for i = 2 : curr_knot_num - 2*degree
            params(1, i) = (0 * KNOTS(1, i + degree - 1) + 2 * KNOTS(1, i + degree))/2;
        end

        [r_sample, wave_sample] = reconstruct_curve_(KNOTS, params, CP_REF, CP_WAVE, degree); %   resample.

    end

    if debug >= 2
        s = sprintf('Final number of knots: %d\n', size(KNOTS, 2));
        disp(s);
    end

%  end of function
end


