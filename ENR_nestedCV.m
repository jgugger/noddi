clc
clear

%% Load data
load data.mat

%% Classification with elastic net regression and nested cross-validation
X = data_table{:,["mTL_ipsi_PC1","mTL_contra_PC1"]}; % predictors
yBinom = (data_table.seizurefree_12m == "SF"); % outcome (SF vs NSF at 12 months)

% Generate example dataset (binary classification)
rng(1); % For reproducibility
nSamples = 30;
nFeatures = 2;

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
end

% Display results
fprintf('Mean accuracy across outer folds: %.2f%%\n', mean(accuracy) * 100);
fprintf('Standard deviation of accuracy: %.2f%%\n', std(accuracy) * 100);