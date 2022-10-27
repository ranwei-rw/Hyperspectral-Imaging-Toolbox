%% Recover the reflectance from an HSZ data structure
%
%% Syntax:
%     R = reconstruct_reflectance(HSZ);
% 
%% Description:
%     Recover the reflectance from an HSZ file 
% 
%% Input:
%     HSZ: HSZ structure
%     
%% Output:
%     R: 3D array corresponding to the reflectance
%
%% See also:
%    reconstruct_illuminant
% 
% This computer code is subject to copyright: (c) National ICT Australia Limited (NICTA) 2014 All Rights Reserved.
% Author: Ran Wei and Antonio Robles-Kelly
% Version: 1.0.7
% Last Update Date: 24 July 2014

function R = reconstruct_reflectance_(HSZ)

    %   get band number
    bands = HSZ.HDR.bands;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Evaluate the HSZ structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch upper(HSZ.HDR.EncodingS)
        case 'NURBS'
            disp('S encoding method: NURBS');
            if HSZ.HDR.IndexedS == 0
                HSZ.S.Elements = eval_nurbs_(HSZ.S.ElementKnots, HSZ.HDR.wavelength, HSZ.S.Elements, HSZ.S.ElementCP, 2); 
            elseif HSZ.HDR.EndmemberIndexedS == 1
                [mats, ~] = size(HSZ.S.Endmembers);
                Q = eval_nurbs_(HSZ.S.EndmembersKnots, HSZ.HDR.wavelength, HSZ.S.Endmembers, HSZ.S.EndmembersCP, 2); 
                HSZ.S.Endmembers = reshape(Q, [mats bands]);
            else
                if HSZ.HDR.IndexedS == 1 && HSZ.HDR.EndmemberIndexedS == 0
                    [mats, ~] = size(HSZ.S.Elements);
                    Q = eval_nurbs_(HSZ.S.ElementKnots, HSZ.HDR.wavelength, HSZ.S.Elements, HSZ.S.ElementCP, 2); 
                    HSZ.S.Elements = reshape(Q, [mats bands]);
                else
                    error('Unknown S index option.');
                end                
            end
        case 'GMM'
            disp('S encoding method: GMM');
            if HSZ.HDR.IndexedS == 0
                HSZ.S.Elements = eval_gaussian_mixture_(HSZ.S.Elements, HSZ.S.ElementsMean, HSZ.S.ElementsStd, HSZ.HDR.wavelength);
            elseif HSZ.HDR.EndmemberIndexedS == 1
                [mats, ~] = size(HSZ.S.Endmembers);
                Q = eval_gaussian_mixture_(HSZ.S.Endmembers, HSZ.S.EndmembersMean, HSZ.S.EndmembersStd, HSZ.HDR.wavelength);
                HSZ.S.Endmembers = reshape(Q, [mats bands]);
            else
                if HSZ.HDR.IndexedS == 1 && HSZ.HDR.EndmemberIndexedS == 0
                    [mats, ~] = size(HSZ.S.Elements);
                    Q = eval_gaussian_mixture_(HSZ.S.Elements, HSZ.S.ElementsMean, HSZ.S.ElementsStd, HSZ.HDR.wavelength);
                    HSZ.S.Elements = reshape(Q, [mats bands]);
                else
                    error('Unknown S index option.');
                end                
            end
        case 'RAW'
            %   do nothing
            disp('S encoding method: RAW');
        otherwise
            error('Unknown S Encoding option.');
    end
    HSZ.HDR.EncodingS = 'RAW';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reconstruct the reflectance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if HSZ.HDR.EndmemberIndexedS == 1
       HSZ.S.Elements = mix_spectra_(HSZ.S.EndmemberAbundances, HSZ.S.EndmemberAbundanceIndexes, HSZ.S.Endmembers);
    end
    
    if HSZ.HDR.IndexedS == 1
       if ndims(HSZ.S.ElementAbundances) == 2
          [rows cols] = size(HSZ.S.ElementAbundances);
          R = mix_spectra_(reshape(HSZ.S.ElementAbundances, [rows*cols 1]), ...
                           reshape(HSZ.S.ElementAbundanceIndexes, [rows*cols 1]), ...
                           HSZ.S.Elements, ...
                          HSZ.HDR.numElementsS);
          R = reshape(R, [rows cols bands]);
       else
          R = mix_spectra_(HSZ.S.ElementAbundances, HSZ.S.ElementAbundanceIndexes, HSZ.S.Elements, HSZ.HDR.numElementsS);
       end
    else
       R = HSZ.S.Elements;
    end
    R = R.*HSZ.S.Factor(:, :, ones(bands, 1));
    
end