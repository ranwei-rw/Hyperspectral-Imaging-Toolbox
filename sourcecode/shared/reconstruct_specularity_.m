%% Evaluate the specular highlights from an HSZ image file
%
%% Syntax:
% K = reconstruct_specularity(HSZ);
% 
%% Description:
% Recover the specularity from a HSZ structure
% 
%% Input:
% HSZ: HSZ structure
%     
%% Output:
% K: 3D array corresponding to the specularities
%
%% See also:
%
% Scyllarus, reconstruct_illuminant, reconstruct_reflectance
%
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2013 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly.
% Version: 1.0.6
% Last Update Date: 25 July 2014

function K = reconstruct_specularity_(HSZ)

    bands = HSZ.HDR.bands;
    
    switch upper(HSZ.HDR.EncodingK)
        case 'NURBS'
            disp('K encoding method: NURBS');
            if HSZ.HDR.IndexedK == 0
                [rows, cols, cps] = size(HSZ.K.Elements);
                clear HSZ.K.Elements;
                HSZ.K.Elements = eval_nurbs_(HSZ.HDR.KnotsElementsK, ...
                                             HSZ.HDR.wavelength, ...
                                             HSZ.K.Elements, ...
                                             reshape(HSZ.K.Elements(rows, cols, :), [1, cps])', ...
                                             2); 
            else
                %   HSZ.HDR.IndexedK == 1
                [mats, cols] = size(HSZ.K.Elements);
                clear HSZ.K.Elements;
                Q = eval_nurbs_(HSZ.HDR.KnotsElementsK, ...
                                HSZ.HDR.wavelength, ...
                                reshape(HSZ.K.Elements(1:mats-1, :), [mats-1 1 cols]), ...
                                HSZ.K.Elements(mats,:)', ...
                                2);   
                HSZ.K.Elements = reshape(Q, [mats-1 bands]);
            end
        case 'GMM'
            disp('K encoding method: GMM');
            if HSZ.HDR.IndexedK == 0
                clear HSZ.K.Elements;
                HSZ.K.Elements = eval_gaussian_mixture_(HSZ.K.Elements(:,:,1:HSZ.HDR.numGMMsK),...
                                                        HSZ.K.Elements(:, :, HSZ.HDR.numGMMsK+1:2*HSZ.HDR.numGMMsK), ...
                                                        HSZ.K.Elements(:, :, 2*HSZ.HDR.numGMMsK+1:3*HSZ.HDR.numGMMsK),...
                                                        HSZ.HDR.wavelength);
            else
                %   IndexedK == 1
                [mats, ~] = size(HSZ.K.Elements);
                clear HSZ.K.Elements;
                Q = eval_gaussian_mixture_(reshape(HSZ.K.Elements(:, 1:HSZ.HDR.numGMMsK),[mats 1 HSZ.HDR.numGMMsK]),...
                                           reshape(HSZ.K.Elements(:, HSZ.HDR.numGMMsK+1:2*HSZ.HDR.numGMMsK), [mats 1 HSZ.HDR.numGMMsK]), ...
                                           reshape(HSZ.K.Elements(:, 2*HSZ.HDR.numGMMsK+1:3*HSZ.HDR.numGMMsK), [mats 1 HSZ.HDR.numGMMsK]),...
                                           HSZ.HDR.wavelength);
                HSZ.K.Elements = reshape(Q, [mats bands]);                
            end
        case 'RAW'
            disp('K encoding method: RAW');
        otherwise
            error('Unknown K encoding method.');
    end
    
    HSZ.HDR.EncodingK = 'RAW'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reconstruct the specularity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if HSZ.HDR.IndexedK == 1
        if ndims(HSZ.K.ElementAbundances) == 2
            [rows cols] = size(HSZ.K.ElementAbundances);
            K = mix_spectra_(reshape(HSZ.K.ElementAbundances, [rows*cols 1]),...
                             reshape(HSZ.K.ElementAbundanceIndexes, [rows*cols 1]),...
                             HSZ.K.Elements, ...
                             HSZ.HDR.numElementsK);
            K = reshape(K, [rows cols bands]);
        else
            K= mix_spectra_(HSZ.K.ElementAbundances,...
                            HSZ.K.ElementAbundanceIndexes,...
                            HSZ.K.Elements, ...
                            HSZ.HDR.numElementsK);
        end
    else
        K = HSZ.K.Elements;
    end
    
    K = K.*HSZ.K.Factor(:, :, ones(bands, 1));

end
