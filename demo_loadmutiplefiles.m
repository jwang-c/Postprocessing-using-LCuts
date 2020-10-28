% Demo code to run multiple files
%% Load file
addpath('support functions');
[dataFile, datapath] = uigetfile({'*.tif';'*.*'},'Load unet result','MultiSelect', 'on');
datanum = size(dataFile,2);
[dataFilev, datapathv] = uigetfile({'*.nii';'*.*'},'Load original image','MultiSelect', 'on');

%% Parameters
sizeLimit = [4/3*pi*4*4*4, 4/3*pi*4*4*4+pi*45*4*4,2,6,45];% adjustable parameter or based on cell information
% [volume_lowerbound, volume_upperbound, r_min, r_max, L]
sigmaD2 = 5; % adjustable parameter to control the rate of exponential decay of similarity measure on point distance
sigmaT2 = 0.2; % adjustable parameter to control the rate of exponential decay of similarity measure on point angle feature
showI = 0;

%%
for d = 1:datanum
    filename = fullfile(datapath, dataFile{1,d});
    V = imread(filename, 1) ;
    info = imfinfo(filename);
    for ii = 2 : size(info, 1)
        temp = imread(filename, ii);
        V = cat(3 , V, temp);
    end

    
    V0_inputFile = fullfile(datapathv, dataFilev{1,d});
    V0info = niftiinfo(V0_inputFile);
    V0_temp = niftiread(V0info);
    V0 = imrotate(V0_temp, 90);
    V0 = flip(V0, 1);
    [post_V,V_needpost,post_radii,post_radiidx,V_noneedforcut] = post_processing_using_lcuts(V,V0,sizeLimit,sigmaD2,sigmaT2,showI);
    
    [filepath,name,ext] = fileparts(V0_inputFile);
    savename = append(name,'_lcuts.mat');
    savepath = strrep(datapath,'instance_seg_result1','LCuts result');
    savefilename = fullfile(savepath,savename);
    save(savefilename,'post_V','post_radii','post_radiidx','V_needpost','V_noneedforcut');
end

