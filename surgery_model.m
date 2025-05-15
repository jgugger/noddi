clc
clear

%% Load data
load surgery_data.mat

%% Create clinical features for alternative model
lesional = surgery_metadata.MRI_lesionType == "MTS";
Concord = surgery_metadata.concordance_Lateralization == "Yes" & surgery_metadata.concordance_Localization == "Yes";
duration = surgery_metadata.durationepi;

%% Perform classification with elastic net regression and nested cross-validation
X = [surgery_metadata.M_ipsi,surgery_metadata.M_contra]; % only M distances
% X = [lesional,Concord,duration]; % only clinical features
yBinom = (surgery_metadata.seizurefree_12m == "SF");

rng(1); % For reproducibility
nSamples = 30;
nFeatures = 2; % adjust

% Define number of folds for outer and inner cross-validation
nOuterFolds = 5;
nInnerFolds = 3;

% Initialize performance metrics
outerCV = cvpartition(yBinom, 'KFold', nOuterFolds);
accuracy = zeros(nOuterFolds, 1);

% Outer loop (model evaluation)
for i = 1:nOuterFolds
    fprintf('Outer fold %d/%d\n', i, nOuterFolds);
    % Split data into training and testing sets
    trainIdx = outerCV.training(i);
    testIdx = outerCV.test(i);
    XTrain = X(trainIdx, :);
    YTrain = yBinom(trainIdx);
    XTest = X(testIdx, :);
    YTest = yBinom(testIdx);

    % Inner loop (hyperparameter tuning)
    innerCV = cvpartition(YTrain, 'KFold', nInnerFolds);
    alphas = [0.1, 0.5, 0.9]; % Elastic net mixing parameters
    lambdas = logspace(-4, 1, 10); % Regularization parameters
    bestAlpha = 0;
    bestLambda = 0;
    bestAccuracy = 0;

    for alpha = alphas
        for lambda = lambdas
            innerAcc = zeros(nInnerFolds, 1);

            for j = 1:nInnerFolds
                innerTrainIdx = innerCV.training(j);
                innerTestIdx = innerCV.test(j);

                XInnerTrain = XTrain(innerTrainIdx, :);
                YInnerTrain = YTrain(innerTrainIdx);
                XInnerTest = XTrain(innerTestIdx, :);
                YInnerTest = YTrain(innerTestIdx);

                % Fit elastic net model
                [B, FitInfo] = lassoglm(XInnerTrain, YInnerTrain, 'binomial', ...
                    'Alpha', alpha, 'Lambda', lambda, 'Standardize', true);

                % Make predictions
                probs = glmval([FitInfo.Intercept; B], XInnerTest, 'logit');
                predictions = probs > 0.5;

                % Calculate accuracy
                innerAcc(j) = mean(predictions == YInnerTest);
            end

            % Average accuracy across inner folds
            meanAcc = mean(innerAcc);
            if meanAcc > bestAccuracy
                bestAccuracy = meanAcc;
                bestAlpha = alpha;
                bestLambda = lambda;
            end
        end
    end

    % Train final model on entire training set with best hyperparameters
    [B, FitInfo] = lassoglm(XTrain, YTrain, 'binomial', ...
        'Alpha', bestAlpha, 'Lambda', bestLambda, 'Standardize', true);

    % Test on outer fold test set
    probs = glmval([FitInfo.Intercept; B], XTest, 'logit');
    predictions = probs > 0.5;
    accuracy(i) = mean(predictions == YTest);
    % Compute AUC
    [~,~,~,aucScores(i)] = perfcurve(YTest, probs, 1);
end

% Bootstrap confidence intervals (95%)
alphaLevel = 0.05;
nBoot = 1000;
bootAcc = bootstrp(nBoot, @mean, accuracy);
bootAUC = bootstrp(nBoot, @mean, aucScores);

accCI = prctile(bootAcc, [100*alphaLevel/2, 100*(1-alphaLevel/2)]);
aucCI = prctile(bootAUC, [100*alphaLevel/2, 100*(1-alphaLevel/2)]);

accMean = mean(accuracy);
aucMean = mean(aucScores);

% Display results
fprintf('Mean accuracy: %.2f%%\n', accMean * 100);
fprintf('95%% CI for accuracy (bootstrap): [%.2f%%, %.2f%%]\n', accCI(1) * 100, accCI(2) * 100);

fprintf('Mean AUC: %.3f\n', aucMean);
fprintf('95%% CI for AUC (bootstrap): [%.3f, %.3f]\n', aucCI(1), aucCI(2));
