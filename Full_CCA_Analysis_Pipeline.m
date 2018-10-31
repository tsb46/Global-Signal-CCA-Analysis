%% Set Path to Packages
addpath(genpath('/Users/tbolt/Desktop/University of Miami/BCCL Research/MATLAB Packages/matlab_cifti'));
addpath(genpath('FSLNets'))
addpath(genpath('palm-alpha105'))
addpath(genpath('NearestSymmetricPositiveDefinite'))


input_beta_folder = '/Users/tbolt/Desktop/Emory_University/GS_CCA/GS_betamap/data/HCP_S1200_GS_betamap/MSM_reg_wbsgrayordinatecortex';
beta_folders = dir(strcat(input_beta_folder,'/*'));

input_motion_folder = '/Users/tbolt/Desktop/Emory_University/GS_CCA/GS_betamap/data/HCP_S1200_FD_DV';
motion_folders = dir(strcat(input_motion_folder,'/*'));

%% Load in All Global Signal Beta Maps
V_N = 96854; % Num of Data Points
N = length(beta_folders) - 3; % Num of Subjects 
% cifti_all = zeros(V_N,N);
% for i = 4:length(beta_folders) % Start at subject 100206
%     disp(beta_folders(i).name)
%     files = dir(strcat(input_beta_folder,'/',beta_folders(i).name,'/*.nii'));
%     cifti_subject = zeros(96854,length(files));
%     for j = 1:length(files)
%         cifti = ciftiopen_mod(strcat(input_beta_folder,'/',beta_folders(i).name,'/',files(j).name));
%         cifti_subject(:,j) = cifti.cdata;
%     end
%     cifti_all(:,i) = mean(cifti_subject,2);
% end

load('cifti_all.mat');

cifti_all(:,1:3) = []; % Remove first three zero columns
%% Load in DVARS
% mean_DVARS_all = zeros(1,N);
% for i = 4:length(beta_folders) % Start at subject 100206
%     sub_folder = dir(strcat(input_motion_folder,'/',motion_folders(i).name,'/r*'));
%     motion_temp = zeros(1,length(sub_folder));
%     for j = 1:length(sub_folder) 
%         DVARS = load(strcat(input_motion_folder,'/',motion_folders(i).name,'/',sub_folder(j).name,'/DVARS.txt'));
%         motion_temp(j) = mean(DVARS(2:end)); 
%     end
%     mean_DVARS_all(i) = mean(motion_temp);
% end
% mean_DVARS_all(1:3) = [];
load('mean_DVARS_all.mat')
%% Identify NaN and non-NaN Nifti Data Points
nan_vertex_indx = find(isnan(cifti_all(:,1)));
data_vertex_indx = find(~isnan(cifti_all(:,1))); 
cifti_all = cifti_all(data_vertex_indx,:);

%% Load in Behavioral Data
[data,labels,all] = xlsread('Behavior/All_Behavior.xlsx');

%% Preprocess Behavioral Data
    % Define Family Structure Variables
    family_structure = all(:,3:7);
    % Define Subject ID Variable
    Subject_ID = all(:,1);
    Subject_ID = cell2mat(Subject_ID(2:end));
    
    %% Dummy Code Variables
        gender = labels(2:end,9); % Gender
        gender_temp = categorical(gender);
        gender_dummy = dummyvar(gender_temp);
        gender_dummy(:,end) = []; % Set reference variable

        
        Recon_Vrs = labels(2:end,11); % Reconstruction Version - HCP switched in the middle of study
        Recon_Vrs_dummy = zeros(N,1);
        for j = 1:length(Recon_Vrs)
            temp = Recon_Vrs{j};
            if isempty(temp)
                Recon_Vrs_dummy(j) = NaN;
            elseif length(temp)==9
                Recon_Vrs_dummy(j) = 0;
            elseif temp=='r227'
                Recon_Vrs_dummy(j) = 1;
            else 
                Recon_Vrs_dummy(j) = 0;
            end  
        end
        
        race = labels(2:end,12); % Race
        race_temp = categorical(race);
        race_dummy = dummyvar(race_temp);
        race_dummy(:,end) = []; % Set reference variable
        
        ethnicity = labels(2:end,13); % Ethnicity
        ethnicity_temp = categorical(ethnicity);
        ethnicity_dummy = dummyvar(ethnicity_temp);
        ethnicity_dummy(:,2) = []; % Set reference variable
        
        
        color_Vis = labels(2:end,339); % Color Vision (not much variability here)
        color_temp = categorical(color_Vis);
        color_dummy = dummyvar(color_temp);
        color_dummy(:,[1 3:4]) = []; % Set reference variable
        
        eye = labels(2:end,340);
        eye_temp = categorical(eye);
        eye_dummy = dummyvar(eye_temp);
        eye_dummy(:,2) = []; % Set reference variable
        
        %% Concatenate original and dummy code variables together
            data_indx = setdiff(1:size(data,2),[1:8 9,11,12,13 127:236 339,340]); % Remove Dummy Coded, Release Version, Acquisition, Task fMRI Performance, Family Structure Variables
            data_2 = horzcat(data(:,data_indx),gender_dummy,Recon_Vrs_dummy,...
                race_dummy,ethnicity_dummy,color_dummy,eye_dummy); % Concatenate data and dummy coded variables
            labels_2 = horzcat(labels(1,data_indx),{'gender','Recon_Vrs','race_dummy1','race_dummy2','race_dummy3',...
                'race_dummy4','race_dummy5','ethnicity_dummy1','ethnicity_dummy2','color','eye'});
        %% Identify "bad" variables e.g. because of bad outliers or not enough distinct values
        % This code is from Steve Smith's CCA Code
        % He had a missing value cut-off of 250 when he was using the 500 sample dataset
        % I increased that to 500 because we're using the 1200 dataset
        % He also uses some cut-off for too many outliers etc.
        
        badvars=[];
        for i=1:size(data_2,2)
          Y=data_2(:,i); grotKEEP=~isnan(Y);  
          grot=(Y(grotKEEP)-median(Y(grotKEEP))).^2; grot=max(grot/mean(grot));  % do we have extreme outliers?
          if (sum(grotKEEP)>500) & (std(Y(grotKEEP))>0) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95) & (grot<100)
              % the 3rd thing above is:  is the size of the largest equal-values-group too large?
            i=i; % do nothing
          else
            [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot]
            badvars=[badvars i];
          end
        end
        
        %% Define Confound Variables
        % Using Smith's definition of confound here
        confounds = [mean_DVARS_all' data_2(:,[311,9,10,17,18,138]) data_2(:,308:309).^(1/3)];
        % In order this is motion, reconstruction version, weight, height,
        % systolic, and diastolic blood pressure,HbA1c, and cube root of total brain
        % volume
        confounds(isnan(confounds))=0;  % impute missing data as zeros
        confounds=[confounds confounds(:,2:end).^2];  % add on squared terms and renormalise
        %% Get the variables we want to keep! 
        varskeep=setdiff([1:size(data_2,2)],[311,9,10,17,18,138,308:309 ... %confound variables
          1,3,4,6,8,11,12,13,14,15,16,136,137,139:154,54,55,62:73,79:83,86,88,92,134,310,115:121, 312:317 ... %remove the "undesirable" variables a la Smith et al.
          badvars]); %bad variables detected above
        % Get final data table
        data_final = data_2(:,varskeep);
        % Get final label variables
        labels_final = labels_2(varskeep);
        
%% Perform PCA of Behavioral Data with Simultaneous Missing Value Imputation
N_Components = 100; % Smith set his to 100 components
% This is the missing data imputation script that Smith uses when he
% performs his PCA.

%%% "impute" missing vars data - actually this avoids any imputation
% data_final=palm_inormal(data_2(:,varskeep)); % Gaussianise

% I (Taylor) added this normalization step here
data_final_z = bsxfun(@minus,data_final,nanmean(data_final,1));
data_final_z = bsxfun(@rdivide,data_final_z,nanstd(data_final_z,[],1));  


for i=1:size(data_final_z,2) % deconfound ignoring missing data
  grot=(isnan(data_final_z(:,i))==0); grotconf=nets_demean(confounds(grot,:)); data_final_z(grot,i)=data_final_z(grot,i)-grotconf*(pinv(grotconf)*data_final_z(grot,i));
end
varsdCOV=zeros(size(data_final_z,1));
for i=1:size(data_final_z,1) % estimate "pairwise" covariance, ignoring missing data
  for j=1:size(data_final_z,1)
    grot=data_final_z([i j],:); grot=cov(grot(:,sum(isnan(grot))==0)'); varsdCOV(i,j)=grot(1,2);
  end
end
varsdCOV2=nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix
[uu,dd]=eigs(varsdCOV2,N_Components);  % SVD (eigs actually)
Beh_PC_Scores=uu-confounds*(pinv(confounds)*uu);    % deconfound again just to be safe
    
%% Perform PCA on Global Signal Beta Maps
[coeff,GS_PC_Scores,~,~,explained] = pca(zscore(cifti_all'),'NumComponents',N_Components);

%% CCA Analysis 
    %% Set up Permutation Scheme
    % This is the script used by Smith
    %%% prepare permutation scheme using PALM - for more details see:
    %%% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/ExchangeabilityBlocks#Generating_the_set_of_permutations_only
    Nperm=10000;                                                                       % in the paper we used 100000 but 10000 should be enough
     [EB,~,missing]=hcp2blocks_mod('Behavior/RESTRICTED_family_structure.csv', [ ], false, Subject_ID); % change the filename to your version of the restricted file
%     PAPset=palm_quickperms([ ], EB, Nperm);   
    load('PAPset_10000.mat');
    %% Remove Participants w/ Missing Family Info from CCA Analysis
    missing_indx = arrayfun(@(x)find(Subject_ID==x,1),missing);
    Beh_PC_Scores(missing_indx,:) = [];
    GS_PC_Scores(missing_indx,:) = [];
    cifti_all(:,missing_indx) = [];
    data_final_z(missing_indx,:)=[];
    
    %% Perform Standard CCA Implemented in MATLAB 
    [GS_weights,Beh_weights,r,GS_CCA_scores,Beh_CCA_scores]=canoncorr(GS_PC_Scores,Beh_PC_Scores);
    
    %% Permutation Testing of CCA Components
    % This is the permutation script used by Smith
    grotRp=zeros(Nperm,N_Components); clear grotRpval;
    for j=1:Nperm
      j
      [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(GS_PC_Scores,Beh_PC_Scores(PAPset(:,j),:));
    end
    for i=1:N_Components;  % get FWE-corrected pvalues
      grotRpval(i)=(1+sum(grotRp(2:end,1)>=r(i)))/Nperm;
    end
    grotRpval
    Ncca=sum(grotRpval<0.05);  % number of FWE-significant CCA components

%% Project CCA Results Back onto Original Variables
N_Vis_Comp = 10; % # of components to look at
orig_GS_CCA_weights = zeros(N_Vis_Comp,size(cifti_all,1));
orig_Beh_CCA_weights = zeros(N_Vis_Comp,size(data_final_z,2));

for i = 1:N_Vis_Comp
    orig_GS_CCA_weights(i,:) = corr(GS_CCA_scores(:,i),cifti_all');
    orig_Beh_CCA_weights(i,:) = corr(Beh_CCA_scores(:,i),data_final_z,'rows','pairwise');
end

%% Visualize GS Weight in Cifti
cifti = ciftiopen_mod('/Users/tbolt/Desktop/Emory_University/GS_CCA/GS_betamap/data/HCP_S1200_GS_betamap/MSM_reg_wbsgrayordinatecortex/100206/100206_rfMRI_Rest1_LR_GS_betamap.dtseries.nii');
weights_cifti = zeros(N_Vis_Comp,V_N);
for i = 1:N_Vis_Comp
   temp_cifti = zeros(1,V_N);
   temp_cifti(nan_vertex_indx) = NaN;
   temp_cifti(data_vertex_indx) = orig_GS_CCA_weights(i,:);
   weights_cifti(i,:) = temp_cifti;
end
cifti.cdata = weights_cifti';
ciftisave_mod(cifti,'CCA_GS_Beh_Results/CCA_GS_Weights_SmithExclusion.dtseries.nii');





