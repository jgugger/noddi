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

M_matrix = zeros(size(icvfR)); % preallocate space for residuals

%% get m-distance for controls using LOO approach
for ii=1:size(icvf_scalars,2) % for each roi
    for i = find(grpCtl)' % location of controls
    tmpicvfR = icvfR;
    tmpicvfR(i,:) = [];
    tmpodiR = odiR;
    tmpodiR(i,:) = [];
    M_matrix(i,ii)= mahal([icvfR(i,ii) odiR(i,ii)],[tmpicvfR(:,ii) tmpodiR(:,ii)]);
    end
end

%% get m-distance for patients
for ii=1:size(icvf_scalars,2) % for each roi
    for i = find(~grpCtl)' % location of patients
    M_matrix(i,ii) = mahal([icvfR(i,ii) odiR(i,ii)],[icvfR(grpCtl,ii) odiR(grpCtl,ii)]);
    end
end

%%
save("mdist.mat","M_matrix","icvfR","odiR","grpEpi","grpCtl","ROIList")

%% Identify location of UL TLE and UL nonlesional TLE
UL_TLE = find(metadata.clinicalHypothesis1_Lateralization == "Right" | metadata.clinicalHypothesis1_Lateralization == "Left");
UL_TLE_nl = find(metadata.clinicalHypothesis1_Lateralization ~= "Bilateral" & metadata.MRI_lesionType == "NA"); % only nonlesional

%% Calculate laterality index (LI) and AUC for M for lateralization of TLE
L_m = M_matrix(:,1:2:17);
R_m = M_matrix(:,2:2:18);

L_ROIList = ROIList(1:2:17);
R_ROIList = ROIList(2:2:18);

LI_m = zeros(116,9); 
epsilon = 1e-6;

for i=1:size(L_ROIList,2)
    LI_m(:,i) = (L_m(:,i)-R_m(:,i))./(L_m(:,i)+R_m(:,i)+epsilon);
end

LI_m_score = sum(LI_m(:,1:9),2);
LI_m = [LI_m LI_m_score];

for i = 1:10
    [X_nl,Y_nl,T_nl,AUC_nl] = perfcurve(metadata.clinicalHypothesis1_Lateralization(UL_TLE_nl),LI_m(UL_TLE_nl,i),'Left','NBoot',1000,'XVals',[0:0.05:1]);
    auc_matrix_nl(i,1:3) = AUC_nl;
    [X,Y,T,AUC] = perfcurve(metadata.clinicalHypothesis1_Lateralization(UL_TLE),LI_m(UL_TLE,i),'Left','NBoot',1000,'XVals',[0:0.05:1]);
    auc_matrix_all(i,1:3) = AUC;
    end

x = ["Amygdala" "Cingulum Hippocampus" "Entorhinal Cortex" "Fornix CST" "Fornix CB" "Hippocampus" "Parahippocampal Gyrus" "Piriform Cortex" "Uncinate Fasciculus" "Sum"];

auc_table_all = table(x',auc_matrix_all(:,1),auc_matrix_all(:,2:3),'VariableNames',["ROI","AUC","CI"]);
auc_table_nl = table(x',auc_matrix_nl(:,1),auc_matrix_nl(:,2:3),'VariableNames',["ROI","AUC","CI"]);

