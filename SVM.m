%% svm: subject identification
%% prepare fixmat
addpath('/home/kampermann/Documents/Code/Matlab/oop/')
clear all
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);
%prepare fixmaps
phases = 1:5;
fix = Fixmat(subjects,phases);
%% collect single trials
scale            = .1;%use .1 to have the typical 2500 pixel maps
final_size       = prod(fix.rect(end-1:end).*scale);
trialnumber      = [400 120 124 240 400];
ttrial           = length(subjects)*sum(trialnumber(phases));
datamatrix       = NaN(final_size,ttrial);
labels.sub       = NaN(1,ttrial);
labels.phase     = NaN(1,ttrial);
labels.trial     = NaN(1,ttrial);
labels.cond      = NaN(1,ttrial);
labels.pos       = NaN(1,ttrial);
v = [];
c=0;
for sub = g.ids'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr};
            fprintf('subject %d phase %d trial %d\n',sub,ph,tr);
            fix.getmaps(v);
            if ~any(isnan(fix.maps(:)))                
                c                   = c+1;
                %scale it if necessary
                if scale ~= 1
                    fix.maps        = imresize(fix.maps,scale,'method','bilinear');
                end                                
                datamatrix(:,c)     = fix.vectorize_maps;
                labels.sub(c)       = sub;
                labels.phase(c)     = ph;
                labels.trial(c)     = tr;
                labels.cond(c)      = unique(fix.deltacsp(fix.selection));
                if ismember(ph,[1 5])
                    labels.pos(c)   = mod(tr,2);%1 is 1st, 0 is 2nd
                else
                    labels.pos(c)   = NaN;
                end
            end
        end
    end
end
%%
% save(sprintf('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/singletrials.mat'),'datamatrix','labels');
%cut the nans
todelete = isnan(sum(datamatrix));
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];

% save('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/datamatrix.mat','datamatrix');
% save('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/labels.mat','labels');
%% PCA
% load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/datamatrix.mat')
% load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/datamatrix.mat','datamatrix')
%compute eigenvectors
fprintf('starting covariance computation\n')
covmat=cov(datamatrix');
fprintf('done\n')
fprintf('starting eigenvector computation\n')
[e dv] = eig(covmat);
fprintf('done\n')
dv = sort(diag(dv),'descend');plot(cumsum(dv)./sum(dv),'o-');xlim([0 200]);
eigen = fliplr(e);
% n where explained variance is > 95%
num = 134;
%collect loadings of every trial
trialload = datamatrix'*eigen(:,1:num)*diag(dv(1:num))^-.5;%dewhitened
% save('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/trialload145dewhite.mat','trialload');
%% matlab's own SVM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% svm for subject identification (test s1 vs s2)
s1 = 34;
s2 = 15;
select = [find(labels.sub == s1),find(labels.sub == s2)];

Y = labels.sub(select)';
X = trialload(select,:);
% actual svm
% Randomly partitions observations into a training set and a test
% set using stratified holdout
P = cvpartition(Y,'Holdout',0.20);
% Use a linear support vector machine classifier
options = optimset('maxiter',1000000,'Display','iter');
svmStruct = svmtrain(X(P.training,:),Y(P.training),'options',options);
% svmStruct = svmtrain(X(P.training,:),Y(P.training),'quadprog_opts',options);
%
C = svmclassify(svmStruct,X(P.test,:));
errRate = sum(Y(P.test)~= C)/P.TestSize 
conMat = confusionmat(Y(P.test),C) %one-vs-one confusion matrix
%% train phase classification
p1 = 1;
p2 = 2;
select = [find(labels.phase == p1),find(labels.phase == p2)];
Y = labels.phase(select)';
X = trialload(select,:);
% actual svm
% Randomly partitions observations into a training set and a test
% set using stratified holdout
P = cvpartition(Y,'Holdout',0.20);
% Use a linear support vector machine classifier
options = optimset('maxiter',1e8,'Display','iter');
svmStruct = svmtrain(X(P.training,:),Y(P.training),'options',options);
C = svmclassify(svmStruct,X(P.test,:));
errRate = sum(Y(P.test)~= C)/P.TestSize 
conMat = confusionmat(Y(P.test),C) %one-vs-one confusion matrix
%% libSVM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load things
if ispc
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/datamatrix.mat')%loads singletrial fixationmaps
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/labels.mat')%loads label struct (sub,phase,trial)
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/alleigen.mat')%load all e
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/trialload145dewhited.mat')%loads trialload of 145 best EV
elseif isunix
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/labels.mat')
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/trialload145dewhited.mat')
end
    
%reorder labels for later convenience during confusion matrix computation
c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end
% %%
% Classified = [];
% Real       = [];
% tfold      = 2;
% for s1 = unique(labels.easy_sub);
%     for s2 = unique(labels.easy_sub);
%         if s1 ~= s2
%             fprintf('%g vs. %g...\n',s1,s2);
%             %
%             select = [find(labels.easy_sub == s1),find(labels.easy_sub == s2)];
%             Y      = labels.easy_sub(select)';
%             X      = trialload(select,:);
%             % actual svm
%             % Randomly partitions observations into a training set and a test
%             % set using stratified holdout
%             P      = cvpartition(Y,'Kfold',tfold);
%             % Use a linear support vector machine classifier
%             for nfold = 1:tfold
%                 svmStruct  = svmtrain(X(P.training(nfold),:),Y(P.training(nfold)),'autoscale',1,'options',options,'kernel_function','linear','method','ls');
%                 Classified = [Classified; svmclassify(svmStruct,X(P.test(nfold),:)) ];
%                 Real       = [Real;Y(P.test(nfold))];
%             end
%         end
%     end
% end
% conMat = confusionmat(Real,Classified);
% conMat = conMat./repmat(sum(conMat,2),1,28);
% subplot(1,2,1);imagesc(conMat);colorbar;axis image
% subplot(1,2,2);imagesc(triu(conMat,1));colorbar;axis image
% %

%% subject classification
%% k-fold subjects (all phases collapsed)
addpath('/home/kampermann/Documents/Code/Matlab/libsvm/matlab/')
if ispc
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/datamatrix.mat')%loads singletrial fixationmaps
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/labels.mat')%loads label struct (sub,phase,trial)
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/alleigen.mat')%load all e
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/trialload145dewhited.mat')%loads trialload of 145 best EV
elseif isunix
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/labels.mat')
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/trialload145dewhited.mat')
end
%reorder labels for later convenience during confusion matrix computation
c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end

kfold   = 10;

ids            = unique(labels.easy_sub);
groupconmat    = NaN(length(ids),length(ids),kfold);
diag_conmat    = NaN(length(ids),length(ids)-1,kfold);
collectcorrect = [];

for s1 = ids
    for s2 = ids
        if s1<s2
            fprintf('Classifying %d vs %d...\n',s1,s2)
            
            select = [find(labels.easy_sub == s1),find(labels.easy_sub == s2)];
            Y = labels.easy_sub(select)';
            X = trialload(select,:);
            
            indices = crossvalind('Kfold',Y,kfold);
            for k = 1:kfold
                test = (indices == k); train = ~test;
                cmd = '-t 0 -c 1 -q';
                model = svmtrain(Y(train), X(train,:), cmd);
                [predicted_label, accuracy, decision_values]  = svmpredict(Y(test), X(test,:), model);
                conMat    = confusionmat(Y(test),predicted_label,'order',[s1,s2]);
                %store in in one-vs-one-confmatrix (rows = class,columns=classhat)
                groupconmat(s1,s2,k) = conMat(1,2)/sum(conMat(1,:),2);%s1 classified as s2
                groupconmat(s2,s1,k) = conMat(2,1)/sum(conMat(2,:),2);%s2 classified as s1
                diag_conmat(s1,min(find(isnan(diag_conmat(s1,:,k)))),k) = conMat(1,1)/sum(conMat(1,:),2);
                diag_conmat(s2,min(find(isnan(diag_conmat(s2,:,k)))),k) = conMat(2,2)/sum(conMat(2,:),2);
                collectcorrect = [collectcorrect accuracy(1)];%over aaaall comparisons
            end
        end
    end
end

%% k-fold subjects (respecting phases)
addpath('/home/kampermann/Documents/Code/Matlab/libsvm/matlab/')
if ispc
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/datamatrix.mat')%loads singletrial fixationmaps
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/labels.mat')%loads label struct (sub,phase,trial)
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/alleigen.mat')%load all e
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/trialload145dewhited.mat')%loads trialload of 145 best EV
elseif isunix
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/labels.mat')
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/trialload145dewhited.mat')
end
%reorder labels for later convenience during confusion matrix computation
c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end

kfold   = 10;

phases         = 1:5;
ids            = unique(labels.easy_sub);
groupconmat    = NaN(length(ids),length(ids),length(phases),kfold);
diag_conmat    = NaN(length(ids),length(ids)-1,length(phases),kfold);
collectcorrect = [];

for ph = phases
    for s1 = ids
        for s2 = ids
            if s1 < s2;
                fprintf('Phase %d - classifying %d vs %d...',ph,s1,s2)
                select = [find(labels.easy_sub == s1),find(labels.easy_sub == s2)];
                select = intersect(select,find(labels.phase == ph));
                Y = labels.easy_sub(select)';
                X = trialload(select,:);
                 indices = crossvalind('Kfold',Y,kfold);
                for k = 1:kfold
                    test = (indices == k); train = ~test;
                    cmd = '-t 0 -c 1 -q';
                    model = svmtrain(Y(train), X(train,:), cmd);
                    [predicted_label, accuracy, decision_values]  = svmpredict(Y(test), X(test,:), model);
                    conMat    = confusionmat(Y(test),predicted_label,'order',[s1,s2]);
                    %store in in one-vs-one-confmatrix (rows = class,columns=classhat)
                    groupconmat(s1,s2,ph,k) = conMat(1,2)/sum(conMat(1,:),2);%s1 classified as s2
                    groupconmat(s2,s1,ph,k) = conMat(2,1)/sum(conMat(2,:),2);%s2 classified as s1
                    diag_conmat(s1,min(find(isnan(diag_conmat(s1,:,ph,k)))),ph,k) = conMat(1,1)/sum(conMat(1,:),2);
                    diag_conmat(s2,min(find(isnan(diag_conmat(s2,:,ph,k)))),ph,k) = conMat(2,2)/sum(conMat(2,:),2);
                    collectcorrect = [collectcorrect accuracy(1)];%over aaaall comparisons
                end
            end
        end
    end
end

%% phases k-fold, all subjects collapsed
addpath('/home/kampermann/Documents/Code/Matlab/libsvm/matlab/')

if ispc
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/trialload145dewhited.mat','trialload');
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/labels.mat');
elseif isunix
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/trialload145dewhited.mat','trialload');
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/labels.mat')
end

kfold = 10;

groupconmat    = NaN(5,5,kfold);
diag_conmat    = NaN(5,4,kfold);
collectcorrect = [];

for p1 = 1:5
    for p2 = 1:5
        if p1<p2;
            fprintf('Classifying phase %d vs %d...',p1,p2)
            select = [find(labels.phase == p1),find(labels.phase == p2)];
            Y = labels.phase(select)';
            X = trialload(select,:);
            indices = crossvalind('Kfold',Y,kfold);
            for k = 1:kfold
                test = (indices == k); train = ~test;
                cmd = '-t 0 -c 1 -q';
                model = svmtrain(Y(train), X(train,:), cmd);
                [predicted_label, accuracy, decision_values]  = svmpredict(Y(test), X(test,:), model);
                conMat    = confusionmat(Y(test),predicted_label,'order',[p1,p2]);
                %one-vs-one confusion matrix
                %store in in one-vs-one-confm (rows = class,columns=classhat)
                groupconmat(p1,p2,k) = conMat(1,2)/sum(conMat(1,:),2);%s1 classified as s2
                groupconmat(p2,p1,k) = conMat(2,1)/sum(conMat(2,:),2);%s2 classified as s1
                diag_conmat(p1,min(find(isnan(diag_conmat(p1,:,k)))),k) = conMat(1,1)/sum(conMat(1,:),2);
                diag_conmat(p2,min(find(isnan(diag_conmat(p2,:,k)))),k) = conMat(2,2)/sum(conMat(2,:),2);
                collectcorrect = [collectcorrect accuracy(1)];%over aaaall comparisons
            end
        end
    end
end

%% phases k-fold, respecting subjects
addpath('/home/kampermann/Documents/Code/Matlab/libsvm/matlab/')

if ispc
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/trialload145dewhited.mat','trialload');
    load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/labels.mat');
elseif isunix
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/trialload145dewhited.mat','trialload');
    load('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/labels.mat')
end

kfold = 10;
ids            = unique(labels.easy_sub);
nsub           = length(ids);
groupconmat    = NaN(5,5,nsub,kfold);
diag_conmat    = NaN(5,5,nsub,kfold);
collectcorrect = [];

for s = 1:nsub
    for i = 1:5
        for j = 1:5
            if i<j;
                fprintf('Sub %d - classifying phase %d vs %d...',s,i,j)
                %get data
                sub = ids(s);
                p1  = i;
for p1 = 1:5
    for p2 = 1:5
        if p1<p2;
            fprintf('Classifying phase %d vs %d...',p1,p2)
            select = [find(labels.phase == p1),find(labels.phase == p2)];
            Y = labels.phase(select)';
            X = trialload(select,:);
            indices = crossvalind('Kfold',Y,kfold);
            for k = 1:kfold
                test = (indices == k); train = ~test;
                cmd = '-t 0 -c 1 -q';
                model = svmtrain(Y(train), X(train,:), cmd);
                [predicted_label, accuracy, decision_values]  = svmpredict(Y(test), X(test,:), model);
                conMat    = confusionmat(Y(test),predicted_label,'order',[p1,p2]);
                %one-vs-one confusion matrix
                %store in in one-vs-one-confm (rows = class,columns=classhat)
                groupconmat(p1,p2,k) = conMat(1,2)/sum(conMat(1,:),2);%s1 classified as s2
                groupconmat(p2,p1,k) = conMat(2,1)/sum(conMat(2,:),2);%s2 classified as s1
                diag_conmat(p1,min(find(isnan(diag_conmat(i,:,k)))),k) = conMat(1,1)/sum(conMat(1,:),2);
                diag_conmat(p2,min(find(isnan(diag_conmat(j,:,k)))),k) = conMat(2,2)/sum(conMat(2,:),2);
                collectcorrect = [collectcorrect accuracy(1)];%over aaaall comparisons
            end
        end
    end
end
                p2  = j;
                select = [find(labels.phase == p1),find(labels.phase == p2)];
                select = intersect(select,find(labels.sub==sub));
                Y = labels.phase(select)';
                X = trialload(select,:);
                indices = crossvalind('Kfold',Y,kfold);
                for k = 1:kfold
                    test = (indices == k); train = ~test;
                    cmd = '-t 0 -c 1 -q';
                    model = svmtrain(Y(train), X(train,:), cmd);
                    [predicted_label, accuracy, decision_values]  = svmpredict(Y(test), X(test,:), model);
                    conMat    = confusionmat(Y(test),predicted_label,'order',[p1,p2]);
                    %one-vs-one confusion matrix
                    %store in in one-vs-one-confm (rows = class,columns=classhat)
                    groupconmat(i,j,s,k) = conMat(1,2)/sum(conMat(1,:),2);%s1 classified as s2
                    groupconmat(j,i,s,k) = conMat(2,1)/sum(conMat(2,:),2);%s2 classified as s1
                    diag_conmat(i,min(find(isnan(diag_conmat(i,:,s,k)))),s,k) = conMat(1,1)/sum(conMat(1,:),2);
                    diag_conmat(j,min(find(isnan(diag_conmat(j,:,s,k)))),s,k) = conMat(2,2)/sum(conMat(2,:),2);
                    collectcorrect = [collectcorrect accuracy(1)];%over aaaall comparisons
                end
            end
        end
    end
end