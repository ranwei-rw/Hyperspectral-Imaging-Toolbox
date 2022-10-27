%% Re-sampling routine for spectral data
%
%% Syntax:
%       [Q NEW_WAVE] = translate_spectra(S, SOURCE_WAVE, TARGET_WAVE);
% 
%% Description:
%       Routine for evaluating the spectra in a SOURCE_WAVE domain different to that
%       in which it was originally sampled. This is done using NURBS
% 
%% Input:
%       S:           Matrix of spectra, where the last dimension corresponds to the SOURCE_WAVE domain.
%       SOURCE_WAVE: Vector containing wavelengths in which the spectra was originally sampled.
%       TARGET_WAVE: Vector containing wavelengths in which the spectra is to be re-sampled.
%     
%% Output:
%       Q: Matrix of re-sampled spectra
%       final_wavelength: Vector of wavelengths used to re-sample Q
% 
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Cong Phuoc Huynh
% Version 1.0.2
% Date: 28 July 2014

function [Q, NEW_WAVE] = translate_spectra_(S, SOURCE_WAVE, TARGET_WAVE)

    [NEW_WAVE, ~] = wavelength_subset_(TARGET_WAVE, SOURCE_WAVE);
    dims_s = ndims(S);
    
    if dims_s == 2
        [rows bands]=size(S);
        R = reshape(S, rows, 1, bands);
    else
        R = S;
    end
    
    [KNOTS, CP_REF, CP_WAVE] = get_nurbs_(R, SOURCE_WAVE, 2);
    
    R = eval_nurbs_(KNOTS, NEW_WAVE, CP_REF, CP_WAVE, 2);
    
    if dims_s == 2
        Q = reshape(R, rows, length(NEW_WAVE));
    else
        Q = R;
    end

end