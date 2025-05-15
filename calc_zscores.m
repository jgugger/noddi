clc
clear

%% Load data
load all_data.mat

% grpCtl: location of controls
% grpEpi: location of patients
% d_age: demeaned age
% sex_female: location of female participants

%% regress age/sex and get residuals
icvfR=zeros(size(icvf_scalars)); % preallocate space for residuals
for ii=1:size(icvf_scalars,2) % for each roi
    mdl = fitlm([d_age(grpCtl),sex_female(grpCtl)],icvf_scalars(grpCtl,ii),'VarNames',{'age','sex','scalar'},'RobustOpts','on'); % 'RobustOpts' 
    icvfR(:,ii)=icvf_scalars(:,ii)-((sex_female*mdl.Coefficients.Estimate(3))+(d_age*mdl.Coefficients.Estimate(2))+mdl.Coefficients.Estimate(1)); 
end

odiR=zeros(size(odi_scalars)); % preallocate space for residuals
for ii=1:size(odi_scalars,2) % for each roi
    mdl = fitlm([d_age(grpCtl),sex_female(grpCtl)],odi_scalars(grpCtl,ii),'VarNames',{'age','sex','scalar'},'RobustOpts','on'); % 'RobustOpts'
    odiR(:,ii)=odi_scalars(:,ii)-((sex_female*mdl.Coefficients.Estimate(3))+(d_age*mdl.Coefficients.Estimate(2))+mdl.Coefficients.Estimate(1)); 
end

z_matrix_icvf = zeros(size(icvfR)); % preallocate space for residuals
z_matrix_odi = zeros(size(odiR)); % preallocate space for residuals

%% calculate z-scores for controls using LOO approach
for ii=1:size(icvf_scalars,2) % for each roi
    for i = find(grpCtl)' % location of controls
    tmpicvfR = icvfR;
    tmpicvfR(i,:) = [];
    tmpodiR = odiR;
    tmpodiR(i,:) = [];
    z_matrix_icvf(i,ii) = (icvfR(i,ii)-mean(tmpicvfR(:,ii)))/std(tmpicvfR(:,ii));
    z_matrix_odi(i,ii) = (odiR(i,ii)-mean(tmpodiR(:,ii)))/std(tmpodiR(:,ii));
    end
end

%% calculate z-score for patients
for ii=1:size(icvf_scalars,2) % for each roi
    for i = find(~grpCtl)' % location of patients
    z_matrix_icvf(i,ii) = (icvfR(i,ii)-mean(icvfR(grpCtl,ii)))/std(icvfR(grpCtl,ii));
    z_matrix_odi(i,ii) = (odiR(i,ii)-mean(odiR(grpCtl,ii)))/std(odiR(grpCtl,ii));
    end
end

%%
save("zscore.mat","z_matrix_icvf","z_matrix_odi")

%% Identify location of UL TLE and UL nonlesional TLE
UL_TLE = find(metadata.clinicalHypothesis1_Lateralization == "Right" | metadata.clinicalHypothesis1_Lateralization == "Left");
UL_TLE_nl = find(metadata.clinicalHypothesis1_Lateralization ~= "Bilateral" & metadata.MRI_lesionType == "NA"); % only nonlesional

%% Calculate laterality index (LI) and AUC for ICVF for lateralization of TLE
L_z_icvf = z_matrix_icvf(:,1:2:17);
R_z_icvf = z_matrix_icvf(:,2:2:18);

L_ROIList = ROIList(1:2:17);
R_ROIList = ROIList(2:2:18);

LI_z_icvf = zeros(116,9); 
epsilon = 1e-6;

for i=1:size(L_ROIList,2)
    LI_z_icvf(:,i) = (L_z_icvf(:,i)-R_z_icvf(:,i))./(abs(L_z_icvf(:,i))+abs(R_z_icvf(:,i))+epsilon);
end

LI_z_icvf_score = sum(LI_z_icvf(:,1:9),2);
LI_z_icvf = [LI_z_icvf LI_z_icvf_score];

for i = 1:10
    [X_nl,Y_nl,T_nl,AUC_nl] = perfcurve(metadata.clinicalHypothesis1_Lateralization(UL_TLE_nl),LI_z_icvf(UL_TLE_nl,i),'Right','NBoot',1000,'XVals',[0:0.05:1]);
    auc_matrix_nl(i,1:3) = AUC_nl;
    [X,Y,T,AUC] = perfcurve(metadata.clinicalHypothesis1_Lateralization(UL_TLE),LI_z_icvf(UL_TLE,i),'Right','NBoot',1000,'XVals',[0:0.05:1]);
    auc_matrix_all(i,1:3) = AUC;
    end

x = ["Amygdala" "Cingulum Hippocampus" "Entorhinal Cortex" "Fornix CST" "Fornix CB" "Hippocampus" "Parahippocampal Gyrus" "Piriform Cortex" "Uncinate Fasciculus" "Sum"];


auc_table_all_icvf = table(x',auc_matrix_all(:,1),auc_matrix_all(:,2:3),'VariableNames',["ROI","AUC","CI"]);
auc_table_nl_icvf = table(x',auc_matrix_nl(:,1),auc_matrix_nl(:,2:3),'VariableNames',["ROI","AUC","CI"]);

%% Calculate laterality index (LI) and AUC for ODI for lateralization of TLE
L_z_odi = z_matrix_odi(:,1:2:17);
R_z_odi = z_matrix_odi(:,2:2:18);

LI_z_odi = zeros(116,9); 

for i=1:size(L_ROIList,2)
    LI_z_odi(:,i) = (L_z_odi(:,i)-R_z_odi(:,i))./(abs(L_z_odi(:,i))+abs(R_z_odi(:,i))+epsilon);
end

LI_z_odi_score = sum(LI_z_odi(:,1:9),2);
LI_z_odi = [LI_z_odi LI_z_odi_score];

for i = 1:10
    [X_nl,Y_nl,T_nl,AUC_nl] = perfcurve(metadata.clinicalHypothesis1_Lateralization(UL_TLE_nl),LI_z_odi(UL_TLE_nl,i),'Right','NBoot',1000,'XVals',[0:0.05:1]);
    auc_matrix_nl(i,1:3) = AUC_nl;
    [X,Y,T,AUC] = perfcurve(metadata.clinicalHypothesis1_Lateralization(UL_TLE),LI_z_odi(UL_TLE,i),'Right','NBoot',1000,'XVals',[0:0.05:1]);
    auc_matrix_all(i,1:3) = AUC;
    end

auc_table_all_odi = table(x',auc_matrix_all(:,1),auc_matrix_all(:,2:3),'VariableNames',["ROI","AUC","CI"]);
auc_table_nl_odi = table(x',auc_matrix_nl(:,1),auc_matrix_nl(:,2:3),'VariableNames',["ROI","AUC","CI"]);