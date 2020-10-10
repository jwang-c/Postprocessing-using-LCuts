% This is the demo code for calling function post_processing_using_lcuts.
% Test/run with matlab R2018a or higher. Some functions might not be
% available in earlier versions.
% Source code is updated at: https://github.com/jwang-c/Postprocessing-using-LCuts
% or it can be available upon request.

% Please cite the following papers if you are refering this code and
% related functions. And please let me know if you have any questions:
% jiewang@virginia.edu.

% [1] Wang J, Batabyal T, Zhang M, Zhang J, Aziz A, Gahlmann A, Acton ST.
% Lcuts: Linear Clustering of Bacteria Using Recursive Graph Cuts. In 2019
% IEEE International Conference on Image Processing (ICIP) 2019 Sep 22 (pp.
% 1575-1579). IEEE.

% [2] Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST, Gahlmann
% A. Non-Invasive Single-Cell Morphometry and Tracking in Living Bacterial
% Biofilms. [BioRxiv preprint]

% Jie Wang, University of Virginia, VIVA lab
% Last update: Oct-10-2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
addpath('support functions')
% results after unet:
filename = '1_SBR_2000_DR_0.94_197_seg.nii';
info = niftiinfo(filename);
% info = niftiinfo('C:/Users/jw4hr/Box Sync/GahlmannActonShare/Biofilm Simulation/underSegmented_Images/shewan_underSegmented_20190819_T1');
V_temp = niftiread(info);
V = imrotate(V_temp, 90);
V = flip(V, 1);
% origial image
V0_inputFile = '1_SBR_2000_DR_0.94_197_deconvolvedwholeexp_T1.nii';
V0info = niftiinfo(V0_inputFile);
V0_temp = niftiread(V0info);
V0 = imrotate(V0_temp, 90);
V0 = flip(V0, 1);
%% run main function
close all;
sizeLimit = [4/3*pi*10*4*4, 4/3*pi*45*4*4,2,6,45];% adjustable parameter or based on cell information, the current choices are estimates
% [volume_lowerbound, volume_upperbound, r_min, r_max, L]
sigmaD2 = 5; % adjustable parameter to control the rate of exponential decay of similarity measure on point distance
sigmaT2 = 0.2; % adjustable parameter to control the rate of exponential decay of similarity measure on point angle feature
showI = 1;
% note: 
% Two figures will be produced if showI = 1: 
% 1.undersegmented clusters (ROIs) that are selected for postprocessing, 
% 2.reconstructed ROIs by fitting geometrical models after LCuts
[post_V,V_needpost,post_radii,post_radiidx] = post_processing_using_lcuts(V,V0,sizeLimit,sigmaD2,sigmaT2,showI);
% Outputs: post_V - postproecssed result (3D image)
%          V_needpost - clusters that need post processing in the image (3D image)
%          post_radii, post_radidx are inscribed spheres radii and index for the point 

