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
scale            = .1;%use .1 to have the typical 2500 pixel maps
final_size       = prod(fix.rect(end-1:end).*scale);
fix.kernel_fwhm  = 30;
trialnumber      = [400 120 124 240 400];
ttrial           = length(subjects)*sum(trialnumber(phases));
datamatrix       = NaN(final_size,ttrial);
labels.sub       = NaN(1,ttrial);
labels.phase     = NaN(1,ttrial);
labels.trial     = NaN(1,ttrial);
v = [];
c=0;
for sub = g.ids'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr};
            fprintf('subject %d phase %d trial %d\n',sub,ph,tr);
            fix.getmaps(v);
            if ~isempty(fix.maps)                
                c                   = c+1;
                %scale it if necessary
                if scale ~= 1
                    fix.maps        = imresize(fix.maps,scale);
                end                                
                datamatrix(:,c)     = fix.vectorize_maps;
                labels.sub(c)       = sub;
                labels.phase(c)     = ph;
                labels.trial(c)     = tr;
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

% save(sprintf('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/singletrials.mat'),'datamatrix','labels');
%% PCA
% load(sprintf('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/singletrials.mat'),'datamatrix','labels');
%compute eigenvectors
fprintf('starting covariance computation\n')
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
options    = optimset('maxiter',100000);
%reorder labels for later convenience during confusion matrix computation
c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end
%%
Classified = [];
Real       = [];
tfold      = 2;
for s1 = unique(labels.easy_sub);
    for s2 = unique(labels.easy_sub);
        if s1 ~= s2
            fprintf('%g vs. %g...\n',s1,s2);
            %
            select = [find(labels.easy_sub == s1),find(labels.easy_sub == s2)];
            Y      = labels.easy_sub(select)';
            X      = trialload(select,:);
            % actual svm
            % Randomly partitions observations into a training set and a test
            % set using stratified holdout
            P      = cvpartition(Y,'Kfold',tfold);
            % Use a linear support vector machine classifier
            for nfold = 1:tfold
                svmStruct  = svmtrain(X(P.training(nfold),:),Y(P.training(nfold)),'autoscale',1,'options',options,'kernel_function','linear','method','ls');
                Classified = [Classified; svmclassify(svmStruct,X(P.test(nfold),:)) ];
                Real       = [Real;Y(P.test(nfold))];
            end
        end
    end
end
conMat = confusionmat(Real,Classified);
conMat = conMat./repmat(sum(conMat,2),1,28);
subplot(1,2,1);imagesc(conMat);colorbar;axis image
subplot(1,2,2);imagesc(triu(conMat,1));colorbar;axis image
%
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
