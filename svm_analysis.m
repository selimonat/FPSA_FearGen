function [joint] = svm_analysis(analysis_type,trialload,labels)
%SVM_ANALYSIS allows different types of analysis classifying subjects,
%phases, conditions,...
%
%  [joint] = SVM_ANALYSIS(analysis_type,trialload,labels) performs the analysis specified by ANALYSIS_TYPE.
%  You can choose from the following options:
%  1 - Classifies subjects without differentiating phases.
%  2 - Classifies subjects, within single phases
%      number of trials is reduced to 95 per phase, to have balanced data (discr>>baseline/cond)     
%  3 - Classifies phases without differentiating subjects.
%      number of trials is not controlled, as ntrials very high anyway 
%  4 - Classifies phases, within single subjects
%      number of trials is reduced to 95 per phase, to have balanced data (discr>>baseline/cond)   
%  5 - Classifies 1st vs. 2nd face within Discr. Tasks (1 and 5 seperately)
%  6 - Classifies Conditions for baseline and testphase (8x8)
%
% SVM_ANALYSIS needs a datamatrix, e.g. containing trialloads, 
% which can be found here:
% .../data/midlevel/singletrialfixmaps/trialload134dewhite.mat', 
% as well as labels, found here:
% .../data/midlevel/singletrialfixmaps/labels.mat'

path = setPath;

r=0; %so far no randomization implemented

nbootstrap    = 2;
cmd           = '-t 0 -c 1 -q'; %t 0: linear, -c 1: criterion, -q: quiet
ids           = unique(labels.easy_sub);
nsub          = length(ids);nsub=3;
numtrials     = [400 120 124 240 400];
phases        = 1:5;

if analysis_type == 1
    name_analysis = 'subjects'; %classify subjects, collapse phases
    savepath      = fullfile(path,[name_analysis '_rand' num2str(r),filesep]);
    mkdir(savepath)
    
    joint = [];
    Classified = uint8([]);
    Real       = uint8([]);
    for n = 1:nbootstrap
        for s1 = 1:nsub
            for s2 = 1:nsub
                if s1 < s2;
                    select                          = ismember(labels.easy_sub,[s1 s2]);
                    Y                               = labels.easy_sub(select)';
                    X                               = trialload(select,:);
                    P                               = cvpartition(Y,'Holdout',.2);%prepares trainings vs testset
                    tic
                    model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
                    [predicted_label, accuracy, ~]  = svmpredict(Y(P.test), X(P.test,:), model);
                    Classified                      = [Classified; uint8(predicted_label) ];
                    Real                            = [Real;uint8(Y(P.test))];
                    fprintf('Run %d - Classifying %d vs %d... in %g seconds ...\n',n,s1,s2,toc)
                end
            end
        end
    end
    dummy         = confusionmat(Real,Classified);
    joint = dummy./sum(dummy(:));
    joint = joint./repmat(sum(joint,2),[1 size(joint,2)]);
    
    save(fullfile(savepath,'joint.mat'),'joint')
    
elseif analysis_type == 2
    name_analysis = 'subjects_inphase'; %classify subjects, respect phases
    savepath      = fullfile(path,[name_analysis '_rand' num2str(r),filesep]);
    mkdir(savepath)
    
    joint = [];
    for ph = phases
        Classified = uint8([]);
        Real       = uint8([]);
        for n = 1:nbootstrap
            ind = labels.phase==ph;
            for s1 = 1:nsub
                for s2 = 1:nsub
                    if s1 < s2;
                        select    = logical(ismember(labels.easy_sub,[s1 s2]).*ind);
                        Y         = labels.easy_sub(select)';
                        X         = trialload(select,:);
                        Y1        = randsample(find(Y == s1),95);
                        Y2        = randsample(find(Y == s2),95);
                        Y         = Y([Y1;Y2]);
                        X         = X([Y1;Y2]);
                        P         = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                        tic
                        model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
                        [predicted_label, accuracy, ~]  = svmpredict(Y(P.test), X(P.test,:), model);
                        Classified                      = [Classified; predicted_label];
                        Real                            = [Real;Y(P.test)];
                        fprintf('Classifying %d vs %d... Fold %d in %g seconds- ...\n',s1,s2,n,toc)
                    end
                end
            end
        end
        dummy         = confusionmat(Real,Classified);
        joint(:,:,ph) = dummy./sum(dummy(:));
        joint(:,:,ph) = joint(:,:,ph)./repmat(sum(joint(:,:,ph),2),[1 size(joint(:,:,ph),2)]);
    end
    
    save(fullfile(savepath,'joint.mat'),'joint')
    
elseif analysis_type == 3
    name_analysis = 'phases';  %classify phases, collapse subjects
    savepath      = fullfile(path,[name_analysis '_rand' num2str(r),filesep]);
    mkdir(savepath)
    
    joint      = [];
    Classified = uint8([]);
    Real       = uint8([]);
    for n = 1:nbootstrap
        for p1 = phases
            for p2 = phases
                if p1 < p2;
                    select                          = ismember(labels.phase,[p1 p2]);
                    Y                               = labels.phase(select)';
                    X                               = trialload(select,:);
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    tic
                    model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
                    [predicted_label, accuracy, ~]  = svmpredict(Y(P.test), X(P.test,:), model);
                    Classified                      = [Classified; uint8(predicted_label)];
                    Real                            = [Real;uint8(Y(P.test))];
                    fprintf('Run %d - Classifying %d vs %d... in %g seconds ...\n',n,p1,p2,toc)
                end
            end
        end
    end
    dummy         = confusionmat(Real,Classified);
    joint = dummy./sum(dummy(:));
    joint = joint./repmat(sum(joint,2),[1 size(joint,2)]);
    
    save(fullfile(savepath,'joint.mat'),'joint')
    
elseif analysis_type == 4
    name_analysis = 'phases_insubject'; %classify subjects, across phases
    savepath      = fullfile(path,[name_analysis '_rand' num2str(r),filesep]);
    mkdir(savepath)
    
    joint   = [];
    for sub = 1:nsub
        Classified = uint8([]);
        Real       = uint8([]);
        for n = 1:nbootstrap
            ind      = labels.easy_sub == sub;
            for p1 = phases
                for p2 = phases
                    if p1 < p2;
                        select    = logical(ismember(labels.phase,[p1 p2]).*ind);
                        Y         = labels.phase(select)';
                        X         = trialload(select,:);
                        Y1        = randsample(find(Y == p1),95);
                        Y2        = randsample(find(Y == p2),95);
                        Y         = Y([Y1;Y2]);
                        X         = X([Y1;Y2]);
                        P         = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                        tic
                        model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
                        [predicted_label, accuracy, ~]  = svmpredict(Y(P.test), X(P.test,:), model);
                        Classified                      = [Classified; predicted_label];
                        Real                            = [Real;Y(P.test)];
                        fprintf('Subject %d - Classifying %d vs %d... in %g seconds- ...\n',sub,p1,p2,toc)
                    end
                end
            end
        end
        dummy         = confusionmat(Real,Classified);
        joint(:,:,sub) = dummy./sum(dummy(:));
        joint(:,:,sub) = joint(:,:,sub)./repmat(sum(joint(:,:,sub),2),[1 size(joint(:,:,sub),2)]);
    end
    save(fullfile(savepath,'joint.mat'),'joint')
    
elseif analysis_type == 5
    name_analysis = '1st_2nd_face'; %classify 1st vs 2nd face in Discr
    savepath      = fullfile(path,[name_analysis '_rand' num2str(r),filesep]);
    mkdir(savepath)
    
    joint   = [];
    pc=0;
    for ph = [1 5]
        Classified = uint8([]);
        Real       = uint8([]);
        pc = pc+1;
        for n = 1:nbootstrap
            ind = labels.phase == ph;
            select    = logical(ismember(labels.pos,[0 1]).*ind);
            Y                               = labels.pos(select)';
            X                               = trialload(select,:);
            P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
            tic
            model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
            [predicted_label, accuracy, ~]  = svmpredict(Y(P.test), X(P.test,:), model);
            Classified                      = [Classified; uint8(predicted_label)];
            Real                            = [Real;uint8(Y(P.test))];
            fprintf('Phase %d - Run %d - Classifying 1st vs 2nd... in %g seconds ...\n',ph,n,toc)
        end
        dummy         = confusionmat(Real,Classified);
        joint(:,:,pc) = dummy./sum(dummy(:));
        joint(:,:,pc) = joint(:,:,pc)./repmat(sum(joint(:,:,pc),2),[1 size(joint(:,:,pc),2)]);
    end
    
    save(fullfile(savepath,'joint.mat'),'joint')
    
elseif analysis_type == 6
    name_analysis = 'conditions'; %classify conditions in FG
    savepath      = fullfile(path,[name_analysis '_rand' num2str(r),filesep]);
    mkdir(savepath)
    
    joint=[];
    pc=0;
    for ph = [2 4]
        ind = labels.phase == ph;
        pc=pc+1;
        for sub = 1:nsub
            Classified = uint8([]);
            Real       = uint8([]);
            ind = logical(ismember(labels.easy_sub,sub).*ind);
            for c1 = -135:45:180
                for c2 = -135:45:180
                    if c1 < c2
                        for n = 1:nbootstrap
                            select    = logical(ismember(labels.cond,[c1 c2]).*ind);
                            Y                               = labels.cond(select)';
                            X                               = trialload(select,:);
                            P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                            tic
                            model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
                            [predicted_label, accuracy, ~]  = svmpredict(Y(P.test), X(P.test,:), model);
                            Classified                      = [Classified; uint8(predicted_label)];
                            Real                            = [Real;uint8(Y(P.test))];
                            fprintf('Phase %d - Run %d - Classifying 1st vs 2nd... in %g seconds ...\n',ph,n,toc)
                        end
                    end
                end
            end
            dummy         = confusionmat(Real,Classified);
            joint(:,:,sub,pc) = dummy./sum(dummy(:));
            joint(:,:,sub,pc) = joint(:,:,pc)./repmat(sum(joint(:,:,pc),2),[1 size(joint(:,:,pc),2)]);
        end
    end
    
    save(fullfile(savepath,'joint.mat'),'joint')
end



    function [path] = setPath
        if ispc || ismac
            [~,version] = GetGit(fullfile(homedir,'Documents','GitHub','fearcloud'));
            path = fullfile(homedir,'Google Drive','EthnoMaster','data','midlevel','svm_analysis',['version' version]);
            mkdir(path)
            addpath([homedir '/Documents/Code/Matlab/libsvm/matlab'])
        elseif isunix
            [~,version] = GetGit(fullfile(homedir,'Documents','Code','Matlab','fearcloud'));
            path = fullfile(homedir,'Documents','fearcloud','data','midlevel','svm_analysis',['version' version]);
            mkdir(path)
            addpath([homedir '/Documents/Code/Matlab/libsvm/matlab'])
        end
    end


end

