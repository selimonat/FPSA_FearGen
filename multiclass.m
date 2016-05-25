clear all

load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\singletrialfixmaps\N27\kernel_defualt_nonull_N_EV_60\labels.mat')
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\singletrialfixmaps\N27\kernel_defualt_nonull_N_EV_60\trialload.mat')
label = Shuffle(labels.easy_sub);
data  = trialload;

addpath('C:/Users/user/Documents/GitHub/LIBSVM-multi-classification/src/')
cmd = '-t 0 -c 1 -q';

%make partition
kfold = 10;
P    = cvpartition(label,'k',kfold);
confusionMatrix = NaN(length(unique(label)),length(unique(label)),kfold);
totalAccuracy   = NaN(kfold,1);

for k = 1:kfold
    trainLabel = label(P.training(k))';
    trainData  = data(P.training(k),:);
    testLabel  = label(P.test(k))';
    testData   = data(P.test(k),:);
    %%
    % #######################
    % Train the SVM in one-vs-rest (OVR) mode
    % #######################
    
    model = ovrtrainBot(trainLabel, trainData, cmd);
    
    % #######################
    % Classify samples using OVR model
    % #######################
    [predict_label, accuracy, decis_values] = ovrpredictBot(testLabel, testData, model);
    [decis_value_winner, label_out] = max(decis_values,[],2);
    %%
    % #######################
    % Make confusion matrix
    % #######################
    [confusionMatrix(:,:,k),order] = confusionmat(testLabel,label_out);
    % Note: For confusionMatrix
    % column: predicted class label
    % row: ground-truth class label
    % But we need the conventional confusion matrix which has
    % column: actual
    % row: predicted
    NTest = length(testLabel);
    figure; imagesc(mean(confusionMatrix,3)');
    xlabel('actual class label');
    ylabel('predicted class label');
    totalAccuracy(k) = trace(confusionMatrix(:,:,k))/NTest;
    fprintf('Total accuracy from the SVM in fold %d: %d percent', k,totalAccuracy*100);
end
save('10foldmulticlass_random','totalAccuracy','confusionMatrix')
%%