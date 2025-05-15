clc
clear

%%
load mdist.mat
load all_data.mat

%%
sum_L = sum(M_matrix(:,1:2:17),2);
sum_R = sum(M_matrix(:,2:2:18),2);
sum_all = sum(M_matrix,2);

med_sum_L_L_deficits = median(sum_L(metadata.NPLateralization=="Left hemisphere deficits"))
iqr_sum_L_L_deficits = prctile(sum_L(metadata.NPLateralization=="Left hemisphere deficits"),[25 75])

med_sum_R_L_deficits = median(sum_R(metadata.NPLateralization=="Left hemisphere deficits"))
iqr_sum_R_L_deficits = prctile(sum_R(metadata.NPLateralization=="Left hemisphere deficits"),[25 75])


p_L_deficits = ranksum(sum_L(metadata.NPLateralization=="Left hemisphere deficits"),...
    sum_R(metadata.NPLateralization=="Left hemisphere deficits"),"tail","right") % alt hypothesis: median x > y

med_sum_L_R_deficits = median(sum_L(metadata.NPLateralization=="Right hemisphere deficits"))
iqr_sum_L_R_deficits = prctile(sum_L(metadata.NPLateralization=="Right hemisphere deficits"),[25 75])

med_sum_R_R_deficits = median(sum_R(metadata.NPLateralization=="Right hemisphere deficits"))
iqr_sum_R_R_deficits = prctile(sum_R(metadata.NPLateralization=="Right hemisphere deficits"),[25 75])

p_R_deficits = ranksum(sum_L(metadata.NPLateralization=="Right hemisphere deficits"),...
    sum_R(metadata.NPLateralization=="Right hemisphere deficits"),"tail","left") % alt hypothesis: median x < y