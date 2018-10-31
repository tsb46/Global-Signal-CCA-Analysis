%% Load in Cifti Data
load('cifti_all.mat');
cifti_all(:,1:3) = []; % Remove first three zero columns

V_N = 96854; % Num of Data Points
N = size(cifti_all,2); % Num of Subjects

%% Identify NaN and non-NaN Nifti Data Points
nan_vertex_indx = find(isnan(cifti_all(:,1)));
data_vertex_indx = find(~isnan(cifti_all(:,1))); 
cifti_all = cifti_all(data_vertex_indx,:);

%% Load in Motion
load('Mean_DVARS_all.mat');


%% First Robustness Analysis: Correlate DVARS w/ GS-Betas at each voxel across subjects
     %% Correlate Mean DVARS and Global Signal Betas at each Voxel
     corr_weights = corr(mean_DVARS_all',cifti_all');
     
     %% Visualize Results
     cifti = ciftiopen_mod('/Users/tbolt/Desktop/Emory_University/GS_CCA/GS_betamap/data/HCP_S1200_GS_betamap/MSM_reg_wbsgrayordinatecortex/100206/100206_rfMRI_REST1_LR_GS_betamap.dtseries.nii');
     weights_cifti = zeros(V_N,1);
     weights_cifti(nan_vertex_indx) = NaN;
     weights_cifti(data_vertex_indx) = corr_weights;
     
     %Save
     cifti.cdata = weights_cifti;
     ciftisave_mod(cifti,'GS_DVARS_Corr.dtseries.nii');
     
%% Second Robustness Analysis: Correlate Canonical Variates w/ DVARS
    %% Remove Participants w/ missing family info to match canonical variates
    load('Behavior/family_missing_indx.mat');
    mean_DVARS_all(missing_indx) = [];
    %% Correlate DVARS and Canonical Variate Scores
    % Load Original CCA Results
    load('CCA_GS_Beh_Results/CCA_Results.mat');
    corr_weights1_GS = corr(mean_DVARS_all',GS_CCA_scores);
    corr_weights1_Beh = corr(mean_DVARS_all',Beh_CCA_scores);
    
    %Load CCA Results w/ Smith Exclusion
    load('CCA_GS_Beh_Results/CCA_Results_Smith_Exclusion.mat');
    corr_weights2_GS = corr(mean_DVARS_all',GS_CCA_scores);
    corr_weights2_Beh = corr(mean_DVARS_all',Beh_CCA_scores);
    
%% Third Robusness Analysis: Regress out GSR from Behavior and Brain
    %% Compute 'mean global signal'
    mean_GS = mean(abs(cifti_all),1);
    % Not 100% sure how you all want this calculated, but it seems to me
    % that the mean of the absolute values of the beta parameters more reflect
    % the 'globalness' of the global signal
    %% Load in Preprocessed and Confound Data
    load('Behavior/All_Behavior_Preprocessed_SmithExclusion.mat');
    load('Behavior/Confounds.mat');
    
    %% Remove Participants w/ missing family 
    mean_GS(missing_indx) = [];
    data_final(missing_indx,:) = [];
    confounds(missing_indx,:) = [];
    cifti_all(:,missing_indx) = [];

    %% Regress Mean GS out of behavior
    data_final_regress = zeros(size(data_final));
    for i = 1:size(data_final,2)
        [~,~,data_final_regress(:,i)] = regress(data_final(:,i),mean_GS');
    end
    
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
    
    %% Perform Standard CCA Implemented in MATLAB 
    [GS_weights_regress,Beh_weights_regress,r_regress,GS_CCA_scores_regress,Beh_CCA_scores_regress]=canoncorr(GS_PC_Scores,Beh_PC_Scores);
    
    %% Project CCA Results Back onto Original Variables
    N_Vis_Comp = 10; % # of components to look at
    orig_GS_CCA_weights = zeros(N_Vis_Comp,size(cifti_all,1));
    orig_Beh_CCA_weights = zeros(N_Vis_Comp,size(data_final_z,2));

    for i = 1:N_Vis_Comp
        orig_GS_CCA_weights(i,:) = corr(GS_CCA_scores_regress(:,i),cifti_all');
        orig_Beh_CCA_weights(i,:) = corr(Beh_CCA_scores_regress(:,i),data_final_z,'rows','pairwise');
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
    ciftisave_mod(cifti,'CCA_GS_Beh_Results/CCA_GS_Weights_SmithExclusion_MeanGSRegress.dtseries.nii');


      


