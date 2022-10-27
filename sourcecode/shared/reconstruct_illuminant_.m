function L = reconstruct_illuminant_(HSZ)
% Syntax:
%    L = reconstruct_illuminant_(HSZ);
% 
% Description:
%    Recover the Lin Guilluminant from an HSZ structure 
% 
% Input:
%    HSZ: HSZ structure (see NICTApipeline help for more details)
%    
% Output:
%    L: 3D array corresponding to the illuminant
% 
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly.
% Version: 1.0.6
% Last Update Date: 6 Feb 2014

    bands=HSZ.HDR.bands;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Evaluate the HSZ structure if necesary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do the NURBS if necesary 
    %Encode the Specularity
    if strcmp(HSZ.HDR.EncodingL, 'NURBS') && HSZ.HDR.IndexedL == 0
        HSZ.HDR.EncodingL = 'RAW';
        % [rows, cols, cps] = size(HSZ.L.Elements);
        % clear HSZ.L.Elements;
        % HSZ.L.Elements = eval_nurbs(HSZ.L.ElementsKnots, HSZ.HDR.wavelength, ...
                                    % HSZ.L.Elements, reshape(HSZ.L.Elements(rows, cols, :), [1, cps])', 2); 
        HSZ.L.Elements = eval_nurbs(HSZ.L.ElementsKnots, HSZ.HDR.wavelength, HSZ.L.Elements, HSZ.L.ElementsCP, 2);                             
    end
    if strcmp(HSZ.HDR.EncodingL, 'NURBS') && HSZ.HDR.IndexedL == 1 && HSZ.HDR.EndmemberIndexedL == 0   
        HSZ.HDR.EncodingL = 'RAW';
        [mats, cols] = size(HSZ.L.Elements);
        clear HSZ.L.Elements;
        % Q = eval_nurbs(HSZ.L.ElementsKnots, HSZ.HDR.wavelength, ...
                       % reshape(HSZ.L.Elements(1:mats-1, :), [mats-1 1 cols]), HSZ.L.Elements(mats, :)', 2);   
        % HSZ.L.Elements = reshape(Q, [mats-1 bands]);
        Q = eval_nurbs(HSZ.L.ElementKnots, HSZ.HDR.wavelength, HSZ.L.Elements, HSZ.L.ElementCP, 2);   
        HSZ.L.Elements = reshape(Q, [mats bands]);
    end
    
    if strcmp(HSZ.HDR.EncodingL, 'NURBS') && HSZ.HDR.EndmemberIndexedL == 1
        HSZ.HDR.EncodingL = 'RAW';
        [mats, cols] = size(HSZ.L.Endmembers);
        clear HSZ.L.Endmembers;
        % Q = eval_nurbs(HSZ.HDR.KnotsEndmembersL, HSZ.HDR.wavelength, ...
           % reshape(HSZ.L.Endmembers(1:mats-1, :), [mats-1 1 cols]), HSZ.L.Endmembers(mats, :)', 2);  
        % HSZ.L.Endmembers = reshape(Q, [mats-1 bands]);
        Q = eval_nurbs(HSZ.L.ElementKnots, HSZ.HDR.wavelength, HSZ.L.Elements, HSZ.L.ElementCP, 2);   
        HSZ.L.Elements = reshape(Q, [mats bands]);        
    end
    
    %Do the GMMs if required
    if strcmp(HSZ.HDR.EncodingL, 'GMM') && HSZ.HDR.IndexedL == 0
        HSZ.HDR.EncodingL = 'RAW';
        clear HSZ.L.Elements;
        % HSZ.L.Elements = eval_gaussian_mixture_(HSZ.L.Elements(:, :, 1:HSZ.HDR.numGMMsL), ...
           % HSZ.L.Elements(:, :, HSZ.HDR.numGMMsL+1:2*HSZ.HDR.numGMMsL), ...
           % HSZ.L.Elements(:, :, 2*HSZ.HDR.numGMMsL+1:3*HSZ.HDR.numGMMsL), ...
           % HSZ.HDR.wavelength);
        HSZ.L.Elements = eval_gaussian_mixture_(HSZ.L.Elements, HSZ.L.ElementsMean, HSZ.L.ElementsStd, HSZ.HDR.wavelength);
     end
     if strcmp(HSZ.HDR.EncodingL, 'GMM')  && HSZ.HDR.IndexedL == 1 && HSZ.HDR.EndmemberIndexedL == 0
        HSZ.HDR.EncodingL = 'RAW';
        [mats, ~] = size(HSZ.L.Elements);
        clear HSZ.L.Elements;
        % Q = eval_gaussian_mixture_(reshape(HSZ.L.Elements(:, 1:HSZ.HDR.numGMMsL), [mats 1 HSZ.HDR.numGMMsL]), ...
           % reshape(HSZ.L.Elements(:, HSZ.HDR.numGMMsL+1:2*HSZ.HDR.numGMMsL), [mats 1 HSZ.HDR.numGMMsL]), ...
           % reshape(HSZ.L.Elements(:, 2*HSZ.HDR.numGMMsL+1:3*HSZ.HDR.numGMMsL), [mats 1 HSZ.HDR.numGMMsL]), ...
           % HSZ.HDR.wavelength);
        Q = eval_gaussian_mixture_(HSZ.L.Elements, HSZ.L.ElementsMean, HSZ.L.ElementsStd, HSZ.HDR.wavelength);   
        HSZ.L.Elements = reshape(Q, [mats bands]);
    end
    if strcmp(HSZ.HDR.EncodingL, 'GMM') && HSZ.HDR.EndmemberIndexedL == 1
        HSZ.HDR.EncodingL = 'RAW';
        [mats, ~] = size(HSZ.L.Endmembers);
        clear HSZ.L.Endmembers;
        % Q = eval_gaussian_mixture_(reshape(HSZ.L.Endmembers(:, 1:HSZ.HDR.numGMMsL), [mats 1 HSZ.HDR.numGMMsL]), ...
           % reshape(HSZ.L.Endmembers(:, HSZ.HDR.numGMMsL+1:2*HSZ.HDR.numGMMsL), [mats 1 HSZ.HDR.numGMMsL]), ...
           % reshape(HSZ.L.Endmembers(:, 2*HSZ.HDR.numGMMsL+1:3*HSZ.HDR.numGMMsL), [mats 1 HSZ.HDR.numGMMsL]), ...
           % HSZ.HDR.wavelength);
        Q = eval_gaussian_mixture_(HSZ.L.Elements, HSZ.L.ElementsMean, HSZ.L.ElementsStd, HSZ.HDR.wavelength); 
        HSZ.L.Endmembers = reshape(Q, [mats bands]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reconstruct the specularity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if HSZ.HDR.EndmemberIndexedL == 1
        HSZ.L.Elements = mix_spectra(HSZ.L.EndmemberAbundances, HSZ.L.EndmemberAbundanceIndexes, HSZ.L.Endmembers);
    end
    if HSZ.HDR.IndexedL == 1
       if (ndims(HSZ.L.ElementAbundances)==2)
          [rows cols]=size(HSZ.L.ElementAbundances);
          L= mix_spectra(reshape(HSZ.L.ElementAbundances, [rows*cols 1]), ...
             reshape(HSZ.L.ElementAbundanceIndexes, [rows*cols 1]), ...
             HSZ.L.Elements, HSZ.HDR.numElementsL);
          L = reshape(L, [rows cols bands]);
       else
          L = mix_spectra(HSZ.L.ElementAbundances, ...
             HSZ.L.ElementAbundanceIndexes, HSZ.L.Elements, HSZ.HDR.numElementsL);
       end
    else
       L=HSZ.L.Elements;
    end

end