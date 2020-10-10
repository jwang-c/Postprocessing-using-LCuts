% This is the major function for post-processing on undersegmented clusters
% using LCuts [1]. There are four major steps as described in [2].

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
% Last update: Jun-06-2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: V - 3D labeled result after initial segmentation method (e.g. u-net)
%         V0 - Initial 3D intensity image
%         sizeLimit - includes volume, radius and length limits
%         sigmaD2 - adjustable parameter to control the rate of exponential
%         decay of similarity measure on distance
%         sigmaD2 - adjustable parameter to control the rate of exponential
%         decay of similarity measure on node direction feature
% Outputs: post_V - postproessed result
%          V_needpost - clusters that need post processing in the image
%          post_radii, post_radidx are inscribed spheres radii and index for the point 
%

function [post_V,V_needpost,post_radii,post_radiidx] = post_processing_using_lcuts(V,V0,sizeLimit, sigmaD2,sigmaT2,showI)

%% global parameters based on the data info
volume_lowerbound = sizeLimit(1);%4/3*pi*5*4.5*4.5;
volume_upperbound = sizeLimit(2);%4/3*pi*60*4.5*4.5;
r_min = sizeLimit(3);%2.5; % for normal: 2; but for low intensity high density, use lower value
r_max = sizeLimit(4);%9;
L = sizeLimit(5);%60;
%% initialization
num_e = max(V(:)); % number of clusters in current result (e: experimental)
need_more_cuts = []; % indices of clusters that need more cuts using lcuts
seg2process = zeros(size(V));

%% Step 1: Data filtering for selecting the region of interest
for n = 1: num_e
    current_seg = V(:,:,:)==n;
    %figure;isosurface(current_seg,0);axis equal;
    % remove the segments with small volumes that are less than a cell size
    if sum(current_seg(:))< volume_lowerbound 
       V(current_seg==1)=0;
    % check coeff. of variance to remove noisy segments
    else 
        CV_scale = plc_findCoeffVar(double(current_seg),double(V0));
        if CV_scale <= 2
    % find undersegmented clusters based on volume limit
             if sum(current_seg(:))> volume_upperbound
                need_more_cuts = [need_more_cuts;n];
                seg2process(current_seg==1) = 1;
             end
        else
             V(current_seg==1)=0;
        end
    end
end
%% Step 2-3: Apply LCuts on each under-segmented clusters
[BW_medial,BW_idx]=bwdist(bwperim(seg2process));
if size(need_more_cuts,1)>0
    clear post_segments  post_radii post_radiidx nodes;
    count = 0;
    for n = 1:size(need_more_cuts,1)
        current_seg = V(:,:,:)== need_more_cuts(n);
        if showI==1
            figure(1),isosurface(current_seg,0);axis equal; hold on; title('Under-segmented clusters');
        end
        % Step 2: Point cloud data generation by constrained medial axis extraction.        
        nodes{n,1} = plc_extractMedialAxis(current_seg,r_min,r_max);
        % figure;plot3(nodes{n,1}(:,1),nodes{n,1}(:,2),nodes{n,1}(:,3),'.')
        % Step 3: Linear Clustering    
        current_postseg = LCuts(nodes{n,1},L,sigmaD2,sigmaT2,0);
        if ~isempty(current_postseg{1,1})
            for i = 1:size(current_postseg,2)
                post_segments{1,count+i} = current_postseg{1,i};
                post_radii{1,count+i} = plc_findMaxRadii(current_postseg{1,i},BW_medial);
                post_radiidx{1,count+i} = plc_findMaxRadii(current_postseg{1,i},BW_idx);

            end
            count = count + size(current_postseg,2); 
        end
     end

%% Step 4: Reconstruct the biofilm by fitting geometrical models
    if count>0 % if there's no output
        Bact = plc_LCuts2Surfaces(post_segments,post_radii,r_min,showI);
    else
        Bact =[];
        post_radii = [];
        post_radiidx =[];
    end
    % put them back to the image and relabel the results
    label = 1;
    post_V = zeros(size(V));
    V_needpost = zeros(size(V));
    % relabeled seg2process =  V_needpost, this is the input for
    % reconstruct the biofilm using convex hull
    label0 = 1;
    for n = 1: num_e
        current_seg = V(:,:,:)==n;
        if sum(current_seg(:))>volume_upperbound
              V_needpost(V(:,:,:)==n)= label0;
              label0 = label0+1;
        end
    end

    % relabeled post_V
    for n = 1: num_e
        current_seg = V(:,:,:)==n;
        if sum(current_seg(:))<=volume_upperbound && sum(current_seg(:))>0
            post_V(V(:,:,:)==n)= label;
            label = label+1; 
        end
    end
    %V_noneedforcut = post_V;
    
    for b = 1:size(Bact,2)
        current_Bact = round(Bact{1,b});
        current_Bact = unique(current_Bact,'row');
        current_Bact(current_Bact<=0)=1;
        current_Bact(current_Bact(:,3)>=size(V,3),3)=size(V,3);
        current_Bact(current_Bact(:,2)>=size(V,2),2)=size(V,2);
        current_Bact(current_Bact(:,1)>=size(V,1),1)=size(V,1);
        if size(current_Bact,1)> volume_lowerbound
            for x = 1:size(current_Bact,1)
                post_V(current_Bact(x,1),current_Bact(x,2),current_Bact(x,3)) = label;
            end
            label = label+1;
        end
    end
else % if not use LCuts
    % re-label for post_V
    post_V = zeros(size(V));
    label = 1;
    for n = 1: num_e
        current_seg = V(:,:,:)==n;
        if sum(current_seg(:))>0
              post_V(V(:,:,:)==n)= label;
              label = label+1;
        end
    end
    %post_V = V;
    post_segments = [];
    Bact = [];
    post_radii = [];
    post_radiidx =[];
    V_needpost = [];
    V_noneedforcut = post_V;
end
