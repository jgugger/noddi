clc
clear

%% Import data --> Need to run calc_mdist.m first to generate this data
load mdist.mat

%%
tiledlayout(1,2)
nexttile
data = [icvfR(grpCtl,11) odiR(grpCtl,11)];  % location of HIPPOCAMPUS_left ROI
% Calculate the mean (centroid) and covariance matrix
mu = mean(data);                % 1x2 mean vector
sigma = cov(data);             % 2x2 covariance matrix
% Plot the data points
% figure;
scatter(data(:,1), data(:,2), 25,[0 0.4627 0.7529], 'filled'); hold on;
% hold on
plot(mu(1), mu(2),'Color',[0.6392 0.0078 0.2039], 'Marker','x', 'MarkerSize', 12, 'LineWidth', 2);  % Plot centroid
% Generate ellipse points
theta = linspace(0, 2*pi, 100);
circle = [cos(theta); sin(theta)];           % Unit circle
[U, S, ~] = svd(sigma);                      % Decompose covariance
k = 2.4477;  % 95% confidence interval for 2D normal (change as needed)
ellipse = (U * sqrt(S) * k * circle) + mu';  % Transform and shift
% Plot the ellipse
plot(ellipse(1,:), ellipse(2,:),'Color',[0 0.4627 0.7529], 'LineWidth', 2);
plot(icvfR(22,11),odiR(22,11),'Color',[0.6392 0.0078 0.2039],'Marker','o', 'MarkerFaceColor', [0.6392 0.0078 0.2039]); % Location of TLE example data
% Labels and aesthetics
title('Distribution of NODDI Parameters for Left Hippocampus in Controls');
legend('Control Data', '', '95% Confidence Ellipse','','Location','northwest');
txt = 'Left TLE Exemplar';
text(-0.057,odiR(22,11),txt,"FontSize",25)
set(gca,'FontSize',25,'LineWidth',2)
xlabel('Neurite Density Index Residual');
ylabel('Orientation Dispersion Index Residual');
nexttile
histogram(M_matrix(grpEpi,11),0:5:60,"FaceColor",[0.6392 0.0078 0.2039])
hold on
histogram(M_matrix(grpCtl,11),0:5:60,"FaceColor",[0 0.4627 0.7529])
ll = xline(M_matrix(22,11),'-',{'Left TLE Exemplar'},'LineWidth',5);
ll.FontSize = 25;
title('Mahalanobis Distance for Left Hippocampus');
set(gca,'FontSize',25,'LineWidth',2)
legend('TLE','Controls')
hold off


