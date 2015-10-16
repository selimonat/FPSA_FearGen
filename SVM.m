%% svm: subject identification
%% prepare fixmat
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
trialnumber = [400 120 124 240 400];
datamatrix       = NaN(2500,length(subjects)*sum(trialnumber(phases)));
labels.sub       = NaN(1,length(subjects)*sum(trialnumber(phases)));
labels.phase     = NaN(1,length(subjects)*sum(trialnumber(phases)));
labels.trial     = NaN(1,length(subjects)*sum(trialnumber(phases)));
v = [];
c=0;
nancount = 0;
for sub = g.ids'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            c=c+1;
            v = {'subject' sub, 'phase' ph 'trialid' tr};
            fprintf('subject %d phase %d trial %d\n',sub,ph,tr)
            fix.kernel_fwhm = 50;
            fix.getmaps(v)
            fix.maps = imresize(fix.maps,[50 50]);
            dummy = fix.vectorize_maps;
            if any(isnan(dummy))
                nancount = nancount+1; %mark this b***
                c = c-1;%forget this trial
            else
%             datamatrix(:,c) = dummy;
            labels.sub(c)       = sub;
            labels.phase(c)     = ph;
            labels.trial(c)     = tr;
            end
        end
    end
end
% save(sprintf('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/singletrials.mat'),'datamatrix','labels');
%cut the nans
datamatrix(:,(end-(nancount-1)):end)=[];
labels.sub(:,(end-(nancount-1)):end)=[];
labels.phase(:,(end-(nancount-1)):end)=[];
labels.trial(:,(end-(nancount-1)):end)=[];
% save(sprintf('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/singletrials.mat'),'datamatrix','labels');
%% PCA
load(sprintf('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/singletrials.mat'),'datamatrix','labels');
%compute eigenvectors
fprintf('starting eigenvector computation\n')
covmat=cov(datamatrix');
fprintf('done\n')
fprintf('starting eigenvector computation\n')
[e dv] = eig(covmat);
fprintf('done\n')
dv = sort(diag(dv),'descend');plot(cumsum(dv)./sum(dv),'o-');xlim([0 100]);
eigen = fliplr(e);
%take only some of them
num = 55;
%collect loadings of every trial
trialload = datamatrix'*eigen(:,1:num);
%% load things
load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/datamatrix.mat')%loads singletrial fixationmaps
load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/labels1500.mat')%loads label struct (sub,phase,trial)
load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/alleigen1500.mat')%load all e
num = 55; eigen = e(:,1:num);
trialload = datamatrix'*eigen(:,1:num);
%% svm for subject identification (test s1 vs s2)
s1 = 7;
s2 = 10;
select = [find(labels.sub == s1),find(labels.sub == s2)];

Y = labels.sub(select)';
Y(find(Y == s1))=1;
Y(find(Y == s2))=-1;
X = trialload(select,:);
% actual svm
% Randomly partitions observations into a training set and a test
% set using stratified holdout
P = cvpartition(Y,'Holdout',0.20);
% Use a linear support vector machine classifier
options = optimset('maxiter',100000);
svmStruct = svmtrain(X(P.training,:),Y(P.training),...
'Kernel_Function','rbf','quadprog_opts',options);
%
C = svmclassify(svmStruct,X(P.test,:));
errRate = sum(Y(P.test)~= C)/P.TestSize  %correct classification rate
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
svmStruct = svmtrain(X(P.training,:),Y(P.training));
C = svmclassify(svmStruct,X(P.test,:));
errRate = 1-sum(Y(P.test)~= C)/P.TestSize  %correct classification rate
conMat = confusionmat(Y(P.test),C) %one-vs-one confusion matrix


%% later - one-vs-one (for all)
groupconmat = NaN(length(g.ids),length(g.ids));
for i = 1:length(g.ids)
    for j = 1:length(g.ids);
        if i<j;
            %get data
            s1 = g.ids(i);
            s2 = g.ids(j);
            select = [find(labels.sub == s1),find(labels.sub == s2)];
            Y = labels.sub(select)';
            X = trialload(select,:);
            P = cvpartition(Y,'Holdout',0.20);%prepares trainings vs testset
            
            options = optimset('maxiter',100000);
            svmStruct = svmtrain(X(P.training,:),Y(P.training),...
                'Kernel_Function','rbf','Method','QP','quadprog_opts',options);
            
            
            C         = svmclassify(svmStruct,X(P.test,:));
            correct   = sum(Y(P.test)== C)/P.TestSize;  %correct classification rate
            conMat    = confusionmat(Y(P.test),C,'order',[s1,s2]);
            %one-vs-one confusion matrix
            %store in in one-vs-one-confm (rows = class,columns=classhat)
            groupconmat(i,i) = conMat(1,1);%s1 classified as s1
            groupconmat(i,j) = conMat(1,2);%s1 classified as s2
            groupconmat(j,i) = conMat(2,1);%s2 classified as s1
            groupconmat(j,j) = conMat(2,2);%s2 classified as s2
        end
    end
end
 

%% actual svm
% Randomly partitions observations into a training set and a test
% set using stratified holdout
P = cvpartition(classes,'Holdout',0.20);
% Use a linear support vector machine classifier
svmStruct = svmtrain(loadings(P.training,:),classes(P.training));
C = svmclassify(svmStruct,loadings(P.test,:));
errRate = 1-sum(classes(P.test)~= C)/P.TestSize  %correct classification rate
conMat = confusionmat(classes(P.test),C) %one-vs-one confusion matrix
