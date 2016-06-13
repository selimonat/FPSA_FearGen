function [result,model,w] = svm_analysis(analysis_type,data,labels,r,nbootstrap,neig)
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
%  7 - Trains on binary data (1 = csp, 0 = rest), tests each cond to be csp
%      or not (takes different fixmaps as input!)
%  8 - same as 7, but subjects collapsed
%
% SVM_ANALYSIS needs a datamatrix, e.g. containing trialloads,
% which can be found here:
% .../data/midlevel/singletrialfixmaps/trialload134dewhite.mat',
% as well as labels, found here:
% .../data/midlevel/singletrialfixmaps/labels.mat'

path = setPath;


% nbootstrap    = 100;
cmd           = '-t 0 -c 1 -q -w1 1 -w0 1'; %t 0: linear, -c 1: criterion, -q: quiet
ids           = unique(labels.easy_sub);
nsub          = length(ids);
numtrials     = [400 120 124 240 400];
phases        = 1:5;
start_time    = [];%init variables here so that they are global
savepath      = [];

if analysis_type == 1
    name_analysis = 'subjects'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result        = [];
    for n = 1:nbootstrap
        Init;
        if r
            labels.easy_sub = Shuffle(labels.easy_sub);
        end
        for s1 = 1:nsub
            for s2 = 1:nsub
                if s1 < s2;
                    select                          = ismember(labels.easy_sub,[s1 s2]);
                    Y                               = labels.easy_sub(select)';
                    X                               = data(select,:);
                    P                               = cvpartition(Y,'Holdout',.2);%prepares trainings vs testset
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,s1,s2,toc,toc(start_time)/60);
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified);
    end
elseif analysis_type == 111
    name_analysis = 'subjects_balanced'; %classify subjects, collapse phases
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result        = [];
    for n = 1:nbootstrap
        Init;
        for s1 = 1:nsub
            for s2 = 1:nsub
                if s1 < s2;
                    select   = ismember(labels.easy_sub,[s1 s2]);
                    Y        = labels.easy_sub(select)';
                    X        = data(select,:); 
                    Y1        = randsample(find(Y == s1),47*5);
                    Y2        = randsample(find(Y == s2),47*5);
                    Y         = Y([Y1;Y2]);
                    X         = X([Y1;Y2],:);
                    randomizeifwanted;
                    P         = cvpartition(Y,'Holdout',.2);%prepares trainings vs testset
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
                        Y1        = randsample(find(Y == s1),47);
                        Y2        = randsample(find(Y == s2),47);
                        Y         = Y([Y1;Y2]);
                        X         = X([Y1;Y2],:);
                        randomizeifwanted;
                        P         = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                        tic;
                        Classify;
                        fprintf('Analysis: %s, Phase %d - Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,ph,n,s1,s2,toc,toc(start_time)/60);
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
                    randomizeifwanted;
                    P                               = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
                    tic
                    Classify;
                    fprintf('Analysis: %s, Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,p1,p2,toc,toc(start_time)/60);
                end
            end
        end
        result(:,:,n) = confusionmat(Real,Classified);
    end
elseif analysis_type == 333
    name_analysis = 'phases_balanced'; %classify subjects, across phases
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
                        Y1        = randsample(find(Y == p1),47);
                        Y2        = randsample(find(Y == p2),47);
                        Y         = Y([Y1;Y2]);
                        X         = X([Y1;Y2],:);
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
                        Y1        = randsample(find(Y == p1),47);
                        Y2        = randsample(find(Y == p2),47);
                        Y         = Y([Y1;Y2]);
                        X         = X([Y1;Y2],:);
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
            
            Init;
            ind                 = labels.phase == ph;
            select              = logical(ismember(labels.pos,[0 1]).*ind);
            Y                   = labels.pos(select)';
            X                   = data(select,:);
            P                   = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
            tic
            Classify;
            fprintf('Analysis: %s, Run %d, %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,1,1,0,toc,toc(start_time)/60);
            result(:,:,1,sub,pc)    = confusionmat(Real,Classified);
            
        end
    end
    
elseif analysis_type == 6
    name_analysis = 'conditions'; %classify conditions in FG
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result=[];
    pc=0;
    for ph = [2 4]
        indph = labels.phase == ph;
        pc=pc+1;
        for sub = 1:nsub
            inds  = logical(ismember(labels.easy_sub,sub).*indph);
            for n = 1:nbootstrap
                Classified = int16([]);
                Real       = int16([]);
                start_time = tic;
                for c1 = -135:45:180
                    for c2 = -135:45:180
                        if c1 < c2
                            select              = logical(ismember(labels.cond,[c1 c2]).*inds);
                            if sum(select) ~= 0
                                Y                   = labels.cond(select)';
                                X                   = data(select,:);
                                P                   = cvpartition(Y,'Holdout',.5); %prepares trainings vs testset
                                tic
                                Classify;
                                fprintf('Analysis: %s,Phase %d Sub %d Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,ph,sub,n,c1,c2,toc,toc(start_time)/60);
                            else
                                fprintf('SKIPPING THIS CLASSIFICATION due to lack of data!!!!!\n');
                                fprintf('Phase %d Sub %d Run %d - Classifying %d vs %d \n',ph,sub,n,c1,c2)
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
elseif analysis_type == 7
    name_analysis = 'conditions_1vsrest_weighted_test'; %classify conditions in FG
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result = nan(8,nbootstrap,nsub,2);
    w = nan(size(data,2),nsub,nbootstrap,2);
    pc = 0;
    for ph = [2 4]
        indph = labels.phase == ph;
        pc = pc+1;
        for sub = 1%:nsub;
             ind  = logical(ismember(labels.easy_sub,sub).*indph);
             Ybin = double(labels.csp(ind))';%binary labelling here
             Ybin0 = Ybin;
             Ycond = double(labels.cond(ind))';
             X = data(ind,:);
             for n = 1:nbootstrap
                 keyboard
                 if r == 1
                     Ybin  = Shuffle(Ybin);%randi(2,[1 length(labels.csp(ind))])-1;
                 end
                 P = cvpartition(Ycond,'Holdout',.5); % respecting conditions
                 model   = svmtrain(Ybin(P.training), X(P.training,:), cmd);
                 try
                     w(:,sub,n,pc)          = model.SVs'*model.sv_coef;
                 end
                 cc=0;
                 for cond = unique(Ycond)'
                     cc=cc+1;
                     if sum(Ycond == cond) ~= 0;
                         fprintf('Analysis %s, Phase %d, Sub %d Run %d - Classifying cond %d... \n',name_analysis,ph,sub,n,cond);
                         i              = logical(P.test.*(Ycond==cond));
                         [predicted]    = svmpredict(Ybin0(i), X(i,:), model);
                         result(cc,n,sub,pc)  = sum(predicted)./length(predicted);%predicted as csp
                     else
                         
                         fprintf('SKIPPING THIS CLASSIFICATION due to lack of data!!!!!\n');
                     end
                 end
             end
        end
    end
elseif analysis_type == 77
    
    name_analysis = 'conditions_CSPvsCSN_test'; %classify conditions in FG
    warning('No randomization implemented yet - press key to run anyway.')    
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result      = nan(8,nbootstrap,nsub,1);
    w           = nan(neig,nsub,nbootstrap,1);
    pc = 0;
    for ph = 4
        indph = labels.phase == ph;%indices of this phase.
        pc = pc+1;%increament the phase counter
        for sub = 1:27;
            ind_all  = logical(ismember(labels.easy_sub,sub).*indph);%this subject, this phase.
            Ycond    = double(labels.cond(ind_all))';%labels of the fixation maps for this subject in this phase.
            X        = data(ind_all,1:neig);%fixation maps of this subject in this phase.
            for n = 1:nbootstrap%
                P       = cvpartition(Ycond,'Holdout',.5); % divide training and test datasets respecting conditions                         
                i       = logical(P.training.*ismember(Ycond,[0 180]));%train using only the CS+ and CS? conditions.
                model   = svmtrain(Ycond(i), X(i,:), cmd);%TRAIN!
                %%
                try
                    w(:,sub,n,pc)          = model.SVs'*model.sv_coef;
                catch
                    keyboard
                end
                %%
                cc=0;
                for cond = unique(Ycond)'
                    cc                          = cc+1;
                    i                           = logical(P.test.*ismember(Ycond,cond));%find all indices that were not used for training belonging to COND.
                    dummy                       = svmpredict(Ycond(i), X(i,:), model);%get the predictions
                    dummy                       = dummy == 0;%binarize: 1=CS+, 0=Not CS+
                    result(cc,n,sub)            = sum(dummy)/length(dummy);
                end                  
            end
        end
    end
    
elseif analysis_type == 8
    warning('this is a test version, probably needs debugging...Press any key to run anyway.')
    pause
    name_analysis = 'conditions_1vsrest_subjcollapsed'; %classify conditions in FG
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result = nan(8,nbootstrap,2);
    pc = 0;
    for ph = [2 4]
        ind = labels.phase == ph;
        pc = pc+1;
        Ybin = double(labels.csp(ind))';%binary labelling here
        Ycond = double(labels.cond(ind))';
        X = data(ind,:);
        for n = 1:nbootstrap
            P = cvpartition(Ycond,'Holdout',.5); % respecting conditions
            model   = svmtrain(Ybin(P.training), X(P.training,:), cmd);
            cc=0;
            for cond = unique(Ycond)'
                cc=cc+1;
                if sum(Ycond == cond) ~= 0;
                    fprintf('Analysis %s, Phase %d, Run %d - Classifying cond %d... \n',name_analysis,ph,n,cond);
                    i                = logical(P.test.*(Ycond==cond));
                    [predicted]      = svmpredict(Ybin(i), X(i,:), model);
                    result(cc,n,pc)  = sum(predicted)./length(predicted);%predicted as csp
                else
                    fprintf('SKIPPING THIS CLASSIFICATION due to lack of data!!!!!\n');
                end
            end
        end
    end
elseif analysis_type == 9
    name_analysis = 'stimID_1vsrest'; %classify conditions in FG
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result = nan(8,nbootstrap,nsub,2);
    pc = 0;
    for ph = [2 4]
        indph = labels.phase == ph;
        pc = pc+1;
        for sub = 1:nsub;
            ind  = logical(ismember(labels.easy_sub,sub).*indph);
            Ybin  = double(labels.stim(ind)==1)';%binary labelling here
            Ycond = double(labels.stim(ind))';
            X = data(ind,:);
            for n = 1:nbootstrap
                P = cvpartition(Ycond,'Holdout',.5); % respecting conditions
                model   = svmtrain(Ybin(P.training), X(P.training,:), cmd);
                cc=0;
                for cond = unique(Ycond)'
                    cc=cc+1;
                    if sum(Ycond == cond) ~= 0;
                        fprintf('Analysis %s, Phase %d, Sub %d Run %d - Classifying cond %d... \n',name_analysis,ph,sub,n,cond);
                        i              = logical(P.test.*(Ycond==cond));
                        [predicted]    = svmpredict(Ybin(i), X(i,:), model);
                        result(cc,n,sub,pc)  = sum(predicted)./length(predicted);%predicted as csp
                    else
                        fprintf('SKIPPING THIS CLASSIFICATION due to lack of data!!!!!\n');
                    end
                end
            end
        end
    end
elseif analysis_type == 10
    name_analysis = 'stimID_1vs1_subjcollapsed'; %classify stimIDs in FG
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result = nan(8,nbootstrap,nsub,2);
    
    for n = 1:nbootstrap
        Init;
        pc = 0;
        for ph = [2 4]
            indph = labels.phase == ph;
            pc = pc+1;
            for s1 = 1:8
                for s2 = 1:8
                    if s1 < s2;
                        select = logical(ismember(labels.stim,[s1 s2]).*indph);
                        Y      = labels.easy_sub(select)';
                        X      = data(select,:);
                        P      = cvpartition(Y,'Holdout',.2);%prepares trainings vs testset
                        %
                        tic
                        Classify;
                        fprintf('Analysis: %s, Phase %d - Run %d - Classifying %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,ph,n,s1,s2,toc,toc(start_time)/60);
                    end
                end
            end
        end
        fprintf('===============\nFinished run %d in %g minutes...\n===============\n',n,toc(start_time)/60);
        result(:,:,n,pc) = confusionmat(Real,Classified);
    end
   
elseif analysis_type == 11
    name_analysis = 'conditions_1vsrest_singlefixations'; %classify conditions in FG
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result = nan(8,nbootstrap,nsub,length(unique(labels.fix)),length(unique(labels.phase)));
    for f = unique(labels.fix)
        indf = labels.fix==f;
        pc = 0;
        for ph = [2 4]
            indph = (labels.phase == ph).*indf;
            pc = pc+1;
            for sub = 1:nsub;
                ind  = logical(ismember(labels.easy_sub,sub).*indph);
                Ybin = double(labels.cond(ind)==0)';%binary labelling here
                Ycond = double(labels.cond(ind))';
                X = data(ind,:);
                for n = 1:nbootstrap
                    P = cvpartition(Ycond,'Holdout',.5); % respecting conditions
                    model   = svmtrain(Ybin(P.training), X(P.training,:), cmd);
                    cc=0;
                    for cond = unique(Ycond)'
                        cc=cc+1;
                        if sum(Ycond == cond) ~= 0;
                            %fprintf('Analysis %s, Phase %d, Sub %d Run %d - Classifying cond %d... \n',name_analysis,ph,sub,n,cond);
                            i              = logical(P.test.*(Ycond==cond));
                            [predicted]    = svmpredict(Ybin(i), X(i,:), model);
                            result(cc,n,sub,f,pc)  = sum(predicted)./length(predicted);%predicted as csp
                        else
                            fprintf('SKIPPING THIS CLASSIFICATION due to lack of data!!!!!\n');
                        end
                    end
                end
            end
        end
    end
elseif analysis_type == 12
    name_analysis = 'CSP_CSN_insubj_test'; %classify CSP vs CSN
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result      = [];
    for n = 1:nbootstrap
        Init;
        pc = 0;
        for ph = [2 4]
            indph = labels.phase == ph;
            pc = pc+1;
            for sub = 1:nsub;
                ind  = logical(ismember(labels.easy_sub,sub).*indph);
                select = logical(ismember(labels.cond,[0 180]).*ind);
                Y = double(labels.cond(select)==0)';%binary labelling here
                X = data(select,:);
                P = cvpartition(Y,'Holdout',.5); % respecting conditions
                tic
                Classify;
                fprintf('Analysis: %s, Phase %d - Run %d - Sub %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,ph,n,sub,toc,toc(start_time)/60);
                result(:,:,n,sub,pc) = confusionmat(Real,Classified,'order',[1 0]);
                w(:,n,sub,pc) = model.SVs'*model.sv_coef;
            end
        end
    end
elseif analysis_type == 13
    name_analysis = 'CSP_CSN_acrosssubj'; %classify CSP vs CSN
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    result      = [];
    for n = 1:nbootstrap
        Init;
        pc = 0;
        for ph = [2 4]
            indph = labels.phase == ph;
            pc = pc+1;
            select = logical(ismember(labels.cond,[0 180]).*indph);
            Y = double(labels.cond(select)==0)';%binary labelling here
            X = data(select,:);
            P = cvpartition(Y,'Holdout',.5); % respecting conditions
            tic
            Classify;
            fprintf('Analysis: %s, Phase %d - Run %d ... in %g seconds, cumulative time %g minutes...\n',name_analysis,ph,n,toc,toc(start_time)/60);
            result(:,:,n,pc) = confusionmat(Real,Classified,'order',[1 0]);
            w(:,n) = model.SVs'*model.sv_coef;
        end
    end
elseif analysis_type ==14
    name_analysis = 'CSPvCSN_insubject'; %classify CSP vs CSN within the subject
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    result   = [];
    for n = 1:nbootstrap
        for sub = 1:nsub
            ind  = logical(ismember(labels.easy_sub,sub));
            Init;
            Y                   = labels.csp(ind)';
            X                   = data(ind,:);
            P                   = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
            tic
            Classify;
            fprintf('Analysis: %s, Run %d, Sub %d %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,sub,1,0,toc,toc(start_time)/60);
            result(:,:,sub,n) = confusionmat(Real,Classified,'order',[1 0]);
            w(:,sub,n)          = model.SVs'*model.sv_coef;
        end
    end
elseif analysis_type ==15
    name_analysis = 'CSPvCSN_collapsesubject'; %classify CSP vs CSN within the subject
    fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
    PrepareSavePath;
    
    result   = [];
    for n = 1:nbootstrap
        Init;
        Y                   = labels.csp';
        X                   = data;
        P                   = cvpartition(Y,'Holdout',.2); %prepares trainings vs testset
        tic
        Classify;
        fprintf('Analysis: %s, Run %d,  %d vs %d... in %g seconds, cumulative time %g minutes...\n',name_analysis,n,1,0,toc,toc(start_time)/60);
        result(:,:,n) = confusionmat(Real,Classified,'order',[1 0]);
        w(:,n)          = model.SVs'*model.sv_coef;
    end
end
try
    save(fullfile(savepath,'result.mat'),'result','model','w')
catch
    save(fullfile(savepath,'result.mat'),'result')
end
    

   function Classify
        model                           = svmtrain(Y(P.training), X(P.training,:), cmd);
        [predicted_label]               = svmpredict(Y(P.test), X(P.test,:), model);
        Classified                      = [Classified; predicted_label];
        Real                            = [Real;Y(P.test)];
    end
    function randomizeifwanted
        if r == 1
           warning('Randomizing labels happening...!')
           Y = Shuffle(Y);
        end
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
        if ismac
            [~,version] = GetGit(fileparts(which(mfilename)));
            path = fullfile('/Volumes/feargen2/project_fearcloud/','data','midlevel','svm_analysis',['version' version]);
            mkdir(path)
            addpath([homedir '/Documents/Code/Matlab/libsvm/matlab'])
         elseif ispc
            t = datestr(now,30);
            path = 'C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\20160531';
            mkdir(path)
            addpath('C:\Users\user\Documents\GitHub\libsvm\matlab\')
        elseif isunix
            [~,version] = GetGit(fullfile(homedir,'Documents','Code','Matlab','fearcloud'));
            path = fullfile(homedir,'Documents','fearcloud','data','midlevel','svm_analysis',['version' version]);
            mkdir(path)
            addpath([homedir '/Documents/Code/Matlab/libsvm/matlab'])
        end
    end

end

