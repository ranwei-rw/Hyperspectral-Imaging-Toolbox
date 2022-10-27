%% Function List
%
%% Scyllarus
% 
% <http://www.scyllarus.com/URL/Scyllarus.html Scyllarus> is the main function 
% in the toolbox. It provides a means to compute a hyperspectral data
% structure (HSZ) which contains photometric parameters and material
% information. It also provides a varied set of processing options and
% spectra representations. Whenever this function is used for research purposes, 
% please cite
%%
% 
% * A. Robles-Kelly and Cong Huynh, Imaging Spectroscopy for Scene Analysis, 
% Springer, 2013.
% 
%% Spectra Representation Functions
%
% <http://www.scyllarus.com/URL/eval_gaussian_mixture.html eval_gaussian_mixture>: Evaluates 
%       spectra represented as a Gaussian mixture.
%
% <http://www.scyllarus.com/URL/eval_nurbs.html eval_nurbs>: Evaluates 
%       spectra represented using NURBS. Whenever this function is used for
%       research purposes, please cite 
%%
% 
% * C. P. Huynh and A. Robles-Kelly, "A NURBS-Based Spectral Reflectance Descriptor 
% with Applications in Computer Vision and Pattern Recognition”, In Proceedings of 
% the IEEE Conference on Computer Vision and Pattern Recognition, 2008.
%
% <http://www.scyllarus.com/URL/get_gaussian_mixture.html get_gaussian_mixture>: Represents
%       spectra using a Gaussian mixture.
%
% <http://www.scyllarus.com/URL/get_nurbs.html get_nurbs>: Represents
%       spectra B-splines (NURBS).
%%
% 
% * C. P. Huynh and A. Robles-Kelly, "A NURBS-Based Spectral Reflectance Descriptor 
% with Applications in Computer Vision and Pattern Recognition”, In Proceedings of 
% the IEEE Conference on Computer Vision and Pattern Recognition, 2008.
% 
%% Input-Output Functions
% 
% <http://www.scyllarus.com/URL/FLAread.html FLAread>: Read a flat (FLA)
%       file from disk.
% 
% <http://www.scyllarus.com/URL/FLAwrite.html FLAwrite>: Write a flat (FLA)
%       file to disk.
%
% <http://www.scyllarus.com/URL/HSZread.html HSZread>: Read a hyperspectral
%       zipped file (HSZ) from disk.
%
% <http://www.scyllarus.com/URL/HSZwrite.html HSZwrite>: Write a hyperspectral
%       zipped file (HSZ) to disk.
%
% <http://www.scyllarus.com/URL/SLZread.html SLZread>: Load a spectral
%       library file from disk.
% 
% <http://www.scyllarus.com/URL/SLZwrite.html SLZwrite>: Write a spectral
%       library file to disk.
%
%% Photogrammetry Functions
% 
% <http://www.scyllarus.com/URL/recover_dichromatic_parameters.html recover_dichromatic_parameters>: Compute 
%       the dichromatic photometric parameters for an input hyperspectral
%       image. Whenever this function is used for research purposes, please cite 
%%
% 
% * L. Gu, A. Robles-Kelly and J. Zhou, “Efficient Estimation of Reflectance Parameters from 
% Imaging Spectroscopy”, In IEEE Transactions on Image Processing. Vol 22(9), 3648-3663, 2013.
%
% <http://www.scyllarus.com/URL/recover_global_illuminant.html recover_global_illuminant>: Recover 
%       the global illuminant in the scene. Whenever this function is used for research purposes, please cite 
%%
% 
% * C. P. Huynh and A. Robles-Kelly, “A Solution of the Dichromatic Model for Multispectral 
% Photometric Invariance”, In International Journal of Computer Vision, Vol. 40(1), 1-27, 2010.
% 
% <http://www.scyllarus.com/URL/recover_materials.html recover_materials>: Recover the 
%       materials and their abundances per pixel in the scene. Whenever this function is used for 
%       research purposes, please cite 
%%
% 
% * C. P. Huynh and A. Robles-Kelly, "A Probabilistic Approach to Spectral Unmixing", In Proceedings 
% of the 13th Int Workshop on Structural and Syntactical Pattern Recognition (S+SSPR 2010), 344-353,2010.
%
% <http://www.scyllarus.com/URL/recover_multi_illuminant.html recover_multi_illuminant>: Recover the 
%       multiple illuminants in the scene.
% 
%% General Purpose Functions
% 
% <http://www.scyllarus.com/URL/encode_HSZ.html encode_HSZ>: Encode a Scyllarus data 
%       structure using spectral descriptors such as NURBS or Gaussian mixtures.
% 
% <http://www.scyllarus.com/URL/eval_HSZ.html eval_HSZ>: Evaluate a Scyllarus data 
%       structure which has been encoded using spectral descriptors such as NURBS or 
%       Gaussian mixtures.
%
% <http://www.scyllarus.com/URL/L2_unmixing.html L2_unmixing>: Perform spectral
%       unmixing using the L2 norm. Whenever this function is used for research purposes, please cite 
%%
% 
% * Y. Qian, S. Jia, J. Zhou and A. Robles-Kelly, “Hyperspectral Unmixing Via L1/2 Sparsity-constrained 
% Nonnegative Matrix Factorization”. In IEEE Transactions on Geoscience and Remote Sensing, 
% Vol. 49(11): 4282-4297, 2011. 
% 
% <http://www.scyllarus.com/URL/mix_spectra.html mix_spectra>: Compute mixed 
%       spectra from a set of end members and abundances.
%
% <http://www.scyllarus.com/URL/reconstruct_illuminant.html reconstruct_illuminant>: Illuminant evaluation 
%       routine for HSZ data structures
%
% <http://www.scyllarus.com/URL/reconstruct_image.html reconstruct_image>: Compute a 
%       flat image data structure from HSZ data.
% 
% <http://www.scyllarus.com/URL/reconstruct_reflectance.html reconstruct_reflectance>: Reflectance 
%       evaluation routine for  HSZ data.
%
% <http://www.scyllarus.com/URL/reconstruct_specularity.html reconstruct_specularity>: Specularity 
%       evaluation routine for  HSZ data.
% 
% <http://www.scyllarus.com/URL/recover_RGBImage.html recover_RGBImage>: Compute a pseudo-colour 
%       image from either an HSZ or flat data structure using camera sensitivity functions.
% 
% <http://www.scyllarus.com/URL/resize_I.html resize_I>: Resize a hyperspectral image cube.
%
% <http://www.scyllarus.com/URL/crop_I.html crop_I>: Crop a hyperspectral image cube.
% 
% <http://www.scyllarus.com/URL/resize_image.html resize_image>: Resize a hyperspectral image 
%       or HSZ data structure.
%
% <http://www.scyllarus.com/URL/crop_image.html crop_image>: Resize a hyperspectral image
%        or HSZ data structure.
%
% <http://www.scyllarus.com/URL/rgb2fla.html rgb2fla>: Convert a trichromatic or monochrome
%       image into a flat data structure (FLA). 
%
% <http://www.scyllarus.com/URL/translate_HSZ.html translate_HSZ>: Re-evaluation routine for HSZ data structures.
% 
% <http://www.scyllarus.com/URL/translate_spectra.html translate_spectra>: Re-evaluation routine for spectral data.
%
% <http://www.scyllarus.com/URL/unmix_illuminant.html unmix_illuminant>: Unmixes the illuminants in an 
%       HSZ data structure so as to refer to the library of end members.
% 
% <http://www.scyllarus.com/URL/unmix_reflectance.html unmix_reflectance>: Unmixes the reflectancde in an 
%       HSZ data structure with respect to a set of end members.
%
% <http://www.scyllarus.com/URL/wavelength_subset.html wavelength_subset>:  Return the indices and values 
%       for the subset of wavelengths that are common to two sets of spectral data
%
% <http://www.scyllarus.com/URL/HSZ2SLZ.html HSZ2SLZ>:  Create a library from an indexed HSZ data structure.
%
