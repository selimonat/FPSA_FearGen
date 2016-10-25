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
dv = sort(diag(dv),'descend');
plot(cumsum(dv)./sum(dv),'o-');xlim([0 200]);
eigen = fliplr(e);
% n where explained variance is > 95%
num = min(find((cumsum(dv)./sum(dv))>.95));
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% classifying SI types of observers, alphas etc
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),subjects);
g               = Group(subjects);
g.getSI(3);
[mat,tag] = g.parameterMat;

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
labels.stim      = NaN(1,ttrial);
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
                labels.stim(c)      = unique(fix.file(fix.selection));
                if ismember(ph,[1 5])
                    labels.pos(c)   = mod(tr,2);%1 is 1st, 0 is 2nd
                else
                    labels.pos(c)   = NaN;
                end
            end
        end
    end
end

%cut the nans
todelete = isnan(sum(datamatrix));
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];
labels.stim(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end
todelete = find(ismember(labels.cond,[500 1000 3000]));
fprintf('Will delete another %g trials...\n',length(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];
labels.stim(:,todelete)=[];
labels.easy_sub(:,todelete)=[];


% save('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/N23/datamatrix.mat','datamatrix');
% save('/home/kampermann/Documents/fearcloud/data/midlevel/singletrialfixmaps/N23/labels.mat','labels');
%% PCA
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
num = 108;
%collect loadings of every trial
trialload = datamatrix'*eigen(:,1:num)*diag(dv(1:num))^-.5;%dewhitened

%%
SI1 = find(g.SI<prctile(g.SI,33));                         %first third
SI2 = find(g.SI>prctile(g.SI,33) & g.SI<prctile(g.SI,66)); %second third
SI3 = find(g.SI>prctile(g.SI,66));                         %third third

labels.SI = nan(1,length(labels.easy_sub));
labels.SI(ismember(labels.easy_sub,SI1))     = 1;
labels.SI(ismember(labels.easy_sub,SI2))     = 2;
labels.SI(ismember(labels.easy_sub,SI3))     = 3;

a =unique([labels.sub' labels.SI'],'rows');
boxplot(g.SI,a(:,2))
t = title('Classes SI');set(t,'FontSize',14);
ylabel('SI')
xlabel('Classes')
%%
alpha_bef = mean(mat(:,[1 3]),2);
alpha_bef1 = find(alpha_bef<prctile(alpha_bef,33));
alpha_bef2= find(alpha_bef>prctile(alpha_bef,33)&alpha_bef<prctile(alpha_bef,66));
alpha_bef3 = find(alpha_bef>prctile(alpha_bef,66));

labels.alpha_bef = nan(1,length(labels.easy_sub));
labels.alpha_bef(ismember(labels.easy_sub,alpha_bef1))  = 1;
labels.alpha_bef(ismember(labels.easy_sub,alpha_bef2))  = 2;
labels.alpha_bef(ismember(labels.easy_sub,alpha_bef3))  = 3;
%%
g = Group(unique(labels.sub));
g.getSI(3)
[mat tags] = g.parameterMat;
alpha_bef = mean(mat(:,[1 3]),2);
alpha_bef_good =  find(alpha_bef<median(alpha_bef));
alpha_bef_bad =  find(alpha_bef>=median(alpha_bef));
labels.alpha_bef2 = nan(1,length(labels.easy_sub));
labels.alpha_bef2(ismember(labels.easy_sub,alpha_bef_good))  = 1;
labels.alpha_bef2(ismember(labels.easy_sub,alpha_bef_bad))  = 0;
a = unique([labels.easy_sub' labels.alpha_bef2'],'rows')
gscatter(a(:,1),alpha_bef,a(:,2));
%
g = Group(unique(labels.sub));
g.getSI(3)
SI_good =  find(g.SI>median(g.SI));
SI_bad  =  find(g.SI<=median(g.SI));
labels.SI2 = nan(1,length(labels.easy_sub));
labels.SI2(ismember(labels.easy_sub,SI_good))  = 1;
labels.SI2(ismember(labels.easy_sub,SI_bad))  = 0;%-1
a = unique([labels.easy_sub' labels.SI2'],'rows')
gscatter(a(:,1),g.SI,a(:,2));

sigmatest_good = find(g.sigma_test<median(g.sigma_test));
sigmatest_bad = find(g.sigma_test>=median(g.sigma_test));
labels.sigmatest2 = nan(1,length(labels.easy_sub));
labels.sigmatest2(ismember(labels.easy_sub,sigmatest_good))  = 1;
labels.sigmatest2(ismember(labels.easy_sub,sigmatest_bad))  = 0;
a = unique([labels.easy_sub' labels.sigmatest2'],'rows')
gscatter(a(:,1),g.sigma_test,a(:,2));


%% plot and analyze svm_analysis confusionmatrices further
%subject classification for 5 phases
% figure
% accuracy = [];
% for i = 1:5
%     a = result(:,:,i);
%     scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,27]);
%     accuracy(i) = mean(diag(scaled));
%     accuracystd(i) = std(diag(scaled));
%     subplot(1,5,i)
%     imagesc(scaled)
%     colorbar
%     caxis([0 1])
%     axis image
%     set(gca,'XTick',[],'YTick',[])
%     colorbar('off')
%     title(['ph' num2str(i)])
% end
% colorbar
% caxis([0 1])
%%
% plot 8cond "classified as csp" (svm_analysis method 7) 
figure;
subplot(1,2,1); [h(1),e(1)] = barwitherr(std(mean(result(:,:,:,1),2),0,3)./sqrt(size(result,3)),mean(mean(result(:,:,:,1),2),3));
title('baseline');ylabel('Classified as CSP')
box off
axis square
set(gca,'XTick',[4 8],'XTickLabel',{'CS+','CS-'})
subplot(1,2,2); [h(2) e(2)] = barwitherr(std(mean(result(:,:,:,2),2),0,3)./sqrt(size(result,3)),mean(mean(result(:,:,:,2),2),3));
title('testphase');ylabel('Classified as CSP')
axH = findall(gcf,'type','axes');
% set(axH,'ylim',[0 .25])
SetFearGenBarColors(h(2));
set(e,'LineWidth',2);
box off
axis square
set(gca,'XTick',[4 8],'XTickLabel',{'CS+','CS-'})
%%

%subjects:
[row col] = GetSubplotNumber(size(result,3));
nsp       = sort([row col]);
figure(1)%baseline
for i = 1:size(result,3)
    
    subplot(nsp(1),nsp(2),i);
    bar(mean(result(:,:,i,1),2));
    axis off
end
figure(2)%testphase
for i = 1:size(result,3)
    subplot(nsp(1),nsp(2),i);
    h = bar(mean(result(:,:,i,2),2));
    axis off
    SetFearGenBarColors(h)
end
%fit Gauss to that
data = squeeze(mean(result(:,:,:,2),2));
for i = 1:size(data,2)
[o(i)]   = FitGauss(-135:45:180,data(:,i),3);
params(i,:) = o(i).Est;
LL(i)    = o(i).Likelihood;
pval(i)  = o(i).pval;
end

figure%testphase gaussians
x_HD = linspace(-135,180,500);
for i = 1:size(result,3)
    subplot(nsp(1),nsp(2),i);
    h = barwitherr(std(result(:,:,i,2),0,2),-135:45:180,mean(result(:,:,i,2),2));
    set(gca,'XTick',[])
    SetFearGenBarColors(h)
    hold on;
    plot(x_HD,o(i).fitfun(x_HD,o(i).Est)+ mean(data(:,i)),'k-','linewidth',2)
    title(num2str(o(i).pval))
    ylim([0 0.3])
    line(xlim,[0.125 0.125],'LineStyle',':')
end


%% %% collect single trials by fixations
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
g = Group(subjects);g.getSI(3);
phases = [2 4];
fix = Fixmat(subjects,phases);
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
labels.stim      = NaN(1,ttrial);
v = [];
c=0;
for nfix = [1 2 3 4];
    for sub = g.ids'
        for ph = [2 4]
            for tr = 1:max(fix.trialid(fix.phase == ph))
                v = {'subject' sub, 'phase' ph 'trialid' tr 'fix' nfix};
                fprintf('Fix %d subject %d phase %d trial %d\n',nfix,sub,ph,tr);
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
                    labels.stim(c)      = unique(fix.file(fix.selection));
                    labels.fix(c)       = nfix;
                    if ismember(ph,[1 5])
                        labels.pos(c)   = mod(tr,2);%1 is 1st, 0 is 2nd
                    else
                        labels.pos(c)   = NaN;
                    end
                end
            end
        end
    end
end
%cut the nans
todelete = isnan(sum(datamatrix));
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];
labels.stim(:,todelete)=[];
labels.fix(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end
todelete = find(labels.cond>180);
fprintf('Will delete another %g trials...\n',length(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];
labels.stim(:,todelete)=[];
labels.easy_sub(:,todelete)=[];
labels.fix(:,todelete)=[];


%% 
a = squeeze(mean(result,3));%bootstraps
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,5]);%scale by rowsums
%%
%insubject:
w0 = w;
w = squeeze(mean(w0(:,:,:,4),3));
num=60;
% fix.getsubmaps;
% fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
% meanmap = Scale(mean(fix.maps,3));

for n = 1:size(w,2);
    hp(:,:,n) = squeeze(eigen(:,1:num)*w(:,n));
end
% for n = 1:27;fix.maps(:,:,n) = reshape(hp(:,n),[50 50]);end
for n = 1:27;fix.maps(:,:,n) = reshape(hp(:,n),[50 50]).*meanmap;end
for n = 1:27;fix.maps(:,:,n) = reshape(hp(:,n),[50 50]);end


%% analyses
% load all the results, from ph4 as well as ph2, 
% both random = 0 and random = 1
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\findingparams\SVM_simulations_ph4_NEV18_fwhm5.mat')
result4_r0 = result;
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\findingparams\SVM_simulations_ph4_NEV18_fwhm5_random.mat')
result4_r1 = result;
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\findingparams\SVM_simulations_ph2_NEV18_fwhm5.mat')
result2_r0 = result;
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\findingparams\SVM_simulations_ph2_NEV18_fwhm5_random.mat')
result2_r1 = result;

%statistics
% csp phase 4 better than chance?
a = squeeze(mean(result4_r0(4,:,:),2));
b = squeeze(mean(result4_r1(4,:,:),2));
[h,p] = ttest(a,b)

%csp phase 2 better than chance?
a = squeeze(mean(result2_r0(4,:,:),2));
b = squeeze(mean(result2_r1(4,:,:),2));
[h,p] = ttest(a,b)

%fit tuning to ph4
%average
data.x = -135:45:180;
data.y = squeeze(mean(mean(result4_r0,2),3))';
data.ids = NaN;
t4 = Tuning(data);
t4.SingleSubjectFit(8)
10.^-t4.fit_results.pval
%don't average
data.y = squeeze(mean(result4_r0,2))'
data.x = repmat(-135:45:180,[26 1]);
data.ids = 1:26;
t4 = Tuning(data);
t4.GroupFit(8)
10.^t4.groupfit.pval
%fit tuning to ph2
data.x = -135:45:180;
data.y = squeeze(mean(mean(result2_r0,2),3))';
data.ids = NaN;
t2 = Tuning(data);
t2.SingleSubjectFit(8)
10.^-t2.fit_results.pval
%don't average
data.y = squeeze(mean(result2_r0,2))'
data.x = repmat(-135:45:180,[26 1]);
data.ids = 1:26;
t2 = Tuning(data);
t2.GroupFit(8)
10.^-t2.groupfit.pval

%because this is also fittable, we check CS+ vs CS- here
a = squeeze(mean(result2_r0(4,:,:),2));
b = squeeze(mean(result2_r0(8,:,:),2));
[h,p] = ttest(a,b)


%% plot that thing @figures_JoV
r2 = result2_r0;
r4 = result4_r0;
r20 = result2_r1;
r40 = result4_r1;