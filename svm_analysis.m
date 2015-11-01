function [result] = svm_analysis(analysis_type,data,labels)
%SVM_ANALYSIS allows different types of analysis classifying subjects,
%phases, conditions,...
%
%  [result] = SVM_ANALYSIS(analysis_type,data,labels) performs the analysis specified by ANALYSIS_TYPE.
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
nsub          = length(ids);nsub=5;;
numtrials     = [400 120 124 240 400];
phases        = 1:5;
start_time    = [];
savepath      = [];

if analysis_type == 1
    name_analysis = 'subjects'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result        = [];
    for n = 1:nbootstrap
        Init;
        for s1 = 1:nsub
            for s2 = 1:nsub
                if s1 < s2;
                    select                          = ismember(labels.easy_sub,[s1 s2]);
                    Y                               = labels.easy_sub(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2);%prepares trainings vs testset
                    %
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,s1,s2,toc,toc(start_time)/60);
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified);
    end
    
elseif analysis_type == 2
    name_analysis = 'subjects_inphase'; %classify subjects, respect phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result = [];
    for ph = phases
        ind = labels.phase==ph;
        for n = 1:nbootstrap
            Init;
            %====
            %actually this inner part of forloop could go to the Classify
            %function. for that would need to homogenize variables across all the
            %analysis cases...
            for s1 = 1:nsub
                for s2 = 1:nsub
                    if s1 < s2;
                        select    = logical(ismember(labels.easy_sub,[s1 s2]).*ind);
                        Y         = labels.easy_sub(select)';
                        X         = data(select,:);
                        Y1        = randsample(find(Y == s1),95);
                        Y2        = randsample(find(Y == s2),95);
                        Y         = Y([Y1;Y2]);
                        X         = X([Y1;Y2]);
                        P         = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                        tic;
                        Classify;
                        fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,s1,s2,toc,toc(start_time)/60);
                    end
                end
            end
            result(:,:,n,ph)  = confusionmat(Real,Classified);
        end
    end
    
elseif analysis_type == 3
    name_analysis = 'phases';  %classify phases, collapse subjects
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result      = [];
    for n = 1:nbootstrap
        Init;
        for p1 = phases
            for p2 = phases
                if p1 < p2;
                    select                          = ismember(labels.phase,[p1 p2]);
                    Y                               = labels.phase(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,p1,p2,toc,toc(start_time)/60);
                end
            end
        end
        result(:,:,n) = confusionmat(Real,Classified);
    end
    
elseif analysis_type == 4
    name_analysis = 'phases_insubject'; %classify subjects, across phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    result   = [];
    for sub = 1:nsub
        ind        = labels.easy_sub == sub;
        for n = 1:nbootstrap
            Init;
            for p1 = phases
                for p2 = phases
                    if p1 < p2;
                        select    = logical(ismember(labels.phase,[p1 p2]).*ind);
                        Y         = labels.phase(select)';
                        X         = data(select,:);
                        Y1        = randsample(find(Y == p1),95);
                        Y2        = randsample(find(Y == p2),95);
                        Y         = Y([Y1;Y2]);
                        X         = X([Y1;Y2]);
                        P         = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                        tic
                        Classify;
                        fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,p1,p2,toc,toc(start_time)/60);
                    end
                end
            end
            result(:,:,n,sub)= confusionmat(Real,Classified);
        end
    end
    
elseif analysis_type == 5
    name_analysis = 'FirstVsSecondFace_insubject'; %classify 1st vs 2nd face during discrimination within subjects
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    result   = [];
    pc=0;
    for ph = [1 5]
        pc = pc+1;
        ind = labels.phase == ph;
        for sub = 1:nsub
            ind  = logical(ismember(labels.easy_sub,sub).*ind);
            for n = 1:nbootstrap
                Init;
                ind                 = labels.phase == ph;
                select              = logical(ismember(labels.pos,[0 1]).*ind);
                Y                   = labels.pos(select)';
                X                   = data(select,:);
                P                   = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                tic
                Classify;
                fprintf('Analysis: %s, Run %d, %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,1,0,toc,toc(start_time)/60);
                result(:,:,n,sub,pc)    = confusionmat(Real,Classified);
            end
        end
    end
    
elseif analysis_type == 6
    name_analysis = 'conditions'; %classify conditions in FG
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    result=[];
    pc=0;
    for ph = [2 4]
        ind = labels.phase == ph;
        pc=pc+1;
        for sub = 1:nsub
            ind  = logical(ismember(labels.easy_sub,sub).*ind);
            for n = 1:nbootstrap
                Init;
                for c1 = -135:45:180
                    for c2 = -135:45:180
                        if c1 < c2
                            select              = logical(ismember(labels.cond,[c1 c2]).*ind);
                            if sum(select) ~= 0
                                Y                   = labels.cond(select)';
                                X                   = data(select,:);
                                P                   = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                                tic
                                Classify;
                                fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,c1,c2,toc,toc(start_time)/60);
                            else
                                fprintf('SKIPPING THIS CLASSIFICATION due to lack of data!!!!!\n');
                            end
                        end
                    end
                end
                if ~any([isempty(Real),isempty(Classified)])
                    result(:,:,n,sub,pc) = confusionmat(Real,Classified);
                else
                    result(:,:,n,sub,pc) = NaN;
                end                
            end
        end
    end
end
    
    save(fullfile(savepath,'result.mat'),'result')
    
    function Classify
    model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
    [predicted_label]               = svmpredict(Y(P.test), X(P.test,:), model);
    Classified                      = [Classified; predicted_label];
    Real                            = [Real;Y(P.test)];
    end
    function Init
    Classified = uint8([]);
    Real       = uint8([]);
    start_time = tic;
    end
    function PrepareSavePath
    savepath      = fullfile(path,[name_analysis '_rand' num2str(r),filesep]);
    if exist(savepath) == 0;mkdir(savepath);end
    fprintf('Created save path: %s\n',savepath);
    end
    
    function [path] = setPath
    if ispc || ismac
        [~,version] = GetGit(fileparts(which(mfilename)));
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

