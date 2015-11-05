function [result] = svm_analysis_behavioral(analysis_type,data,labels)
%SVM_ANALYSIS allows different types of analysis classifying subjects,
%phases, conditions,...
%
%  [result] = SVM_ANALYSIS(analysis_type,data,labels) performs the analysis specified by ANALYSIS_TYPE.
%  You can choose from the following options:
%  1 - Classifies by SI (1 / 2 / 3, 3 is high SI)
%  2 - Classifies by alpha before, phase 1 (1, 2, 3 - 3 is high alpha)
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
nsub          = length(ids);
phases        = 1:5;
start_time    = [];%init variables here so that they are global
savepath      = [];

if analysis_type == 1
    name_analysis = 'subjects_by_SI'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    ind = labels.phase == 4;
    
    result        = [];
    for n = 1:nbootstrap
        Init;
        for c1 = unique(labels.SI)
            for c2 = unique(labels.SI)
                if c1 < c2;
                    select    = logical(ismember(labels.SI,[c1 c2]).*ind);
                    Y                               = labels.SI(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    %
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,c1,c2,toc,toc(start_time)/60);
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified);
    end
    
elseif analysis_type == 2
    name_analysis = 'subjects_by_alpha_ave'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    ind = ismember(labels.phase,[1 5]);
    result        = [];
    for n = 1:nbootstrap
        Init;
        for c1 = unique(labels.alpha_ave)
            for c2 = unique(labels.alpha_ave)
                if c1 < c2;
                    select    = logical(ismember(labels.alpha_ave,[c1 c2]).*ind);
                    Y                               = labels.alpha_ave(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    %
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,c1,c2,toc,toc(start_time)/60);
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified);
    end
elseif analysis_type == 3
    name_analysis = 'subjects_by_alpha_bef'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    ind = labels.phase == 1;
    result        = [];
    for n = 1:nbootstrap
        Init;
        for c1 = unique(labels.alpha_bef)
            for c2 = unique(labels.alpha_bef)
                if c1 < c2;
                    select    = logical(ismember(labels.alpha_bef,[c1 c2]).*ind);
                    Y                               = labels.alpha_bef(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    %
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,c1,c2,toc,toc(start_time)/60);
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified);
    end
elseif analysis_type == 3
    name_analysis = 'subjects_by_alpha_aft'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    ind = labels.phase == 5;
    result        = [];
    for n = 1:nbootstrap
        Init;
        for c1 = unique(labels.alpha_aft)
            for c2 = unique(labels.alpha_aft)
                if c1 < c2;
                    select    = logical(ismember(labels.alpha_aft,[c1 c2]).*ind);
                    Y                               = labels.alpha_aft(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    %
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,c1,c2,toc,toc(start_time)/60);
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified);
    end
elseif analysis_type == 3
    name_analysis = 'subjects_by_alpha_aft'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    ind = labels.phase == 5;
    result        = [];
    for n = 1:nbootstrap
        Init;
        for c1 = unique(labels.alpha_aft)
            for c2 = unique(labels.alpha_aft)
                if c1 < c2;
                    select    = logical(ismember(labels.alpha_aft,[c1 c2]).*ind);
                    Y                               = labels.alpha_aft(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    %
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,c1,c2,toc,toc(start_time)/60);
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified);
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

