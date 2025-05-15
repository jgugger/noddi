clc
clear

%% Import data --> Need to run calc_zscores.m and calc_mdist.m first to generate this data
load zscore.mat
load mdist.mat
load all_data.mat

%%
metadata.clinicalHypothesis1_Lateralization = fillmissing(metadata.clinicalHypothesis1_Lateralization,'constant',"Control");

x = ["L Amygdala" "L Cingulum Hippocampus" "L Entorhinal Cortex" "L Fornix CST" "L Fornix CB" "L Hippocampus" "L Parahippocampal Gyrus" "L Piriform Cortex" "L Uncinate Fasciculus"];
y = ["R Amygdala" "R Cingulum Hippocampus" "R Entorhinal Cortex" "R Fornix CST" "R Fornix CB" "R Hippocampus" "R Parahippocampal Gyrus" "R Piriform Cortex" "R Uncinate Fasciculus"];

L_icvf_z = z_matrix_icvf(:,1:2:17);
R_icvf_z = z_matrix_icvf(:,2:2:18);

L_odi_z = z_matrix_odi(:,1:2:17);
R_odi_z = z_matrix_odi(:,2:2:18);

L_m = M_matrix(:,1:2:17);
R_m = M_matrix(:,2:2:18);

f1 = figure;
f2 = figure;

%% 
figure(f1);

tiledlayout(3,9)
nexttile
    boxplot(L_icvf_z(:,1),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    title(x(1))
    ylabel('NDI z-score')
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-4.25 3.5])
for i=2:9
    nexttile
    boxplot(L_icvf_z(:,i),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    title(x(i))
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-4.25 3.5])
end
nexttile
    boxplot(L_odi_z(:,1),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    ylabel('ODI z-score')
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-3 3])
for i=2:9
    nexttile
    boxplot(L_odi_z(:,i),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-3 3])
end
nexttile
    boxplot(L_m(:,1),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    ylabel('M-distance')
    set(gca,"FontSize",15,'YLim',[0 35])
for i=2:9
    nexttile
    boxplot(L_m(:,i),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    set(gca,"FontSize",15,'YLim',[0 35])
end

%%
figure(f2);

tiledlayout(3,9)
nexttile
    boxplot(R_icvf_z(:,1),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    title(y(1))
    ylabel('NDI z-score')
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-4.25 3.5])
for i=2:9
    nexttile
    boxplot(R_icvf_z(:,i),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    title(y(i))
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-4.25 3.5])
end
nexttile
    boxplot(R_odi_z(:,1),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    ylabel('ODI z-score')
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-3 3])
for i=2:9
    nexttile
    boxplot(R_odi_z(:,i),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    set(gca,"FontSize",15,'XTickLabel',{' '},'YLim',[-3 3])
end
nexttile
    boxplot(R_m(:,1),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    ylabel('M-distance')
    set(gca,"FontSize",15,'YLim',[0 35])
for i=2:9
    nexttile
    boxplot(R_m(:,i),metadata.clinicalHypothesis1_Lateralization,"Notch","on","GroupOrder",["Control","Left","Right","Bilateral"],"Colors",'k','Symbol','+k')
    set(gca,"FontSize",15,'YLim',[0 35])
end

