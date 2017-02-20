function [R AVEHP result HP]=SVM_simulations_BDNF(phase,varargin)
%%
%This script is about SVM training CS+ vs. CS- for phases 2 and 4 based on
%our discussions. It collects the data and computes the eigenvalues on the
%fly for different sorts of parameters (kernel_fwhm, number of
%eigenvalues). As LK pointed out the number of trials are lower in the
%baseline, all the test-training sessions should use the same number of
%trials. For example to keep the comparisons comparable, in baseline 11
%trials in the baseline, with .5 hold-out, should correspond to hold-out of
%.25 in the test. This is taken into account by the HOLDOUT_RATIO. This
%script will run train-test cycles for different smootheness levels and
%with different numbers of eigenvectors (up to TEIG).
random = 0;

tbootstrap       = 1000;%number of bootstraps
% phase            = 2;%baseline or test
holdout_ratio    = [NaN .5 NaN .5];%see above %1/6 is half of 1/3 (11 vs 33 trials)
teig             = 100;%how many eigenvalues do you want to include for each run.
fwhm             = 0;%counter
R                = [];%result storage
AVE              = [];%result storate
HP               = [];
eigval           = [];
sample120first = 1;

subjects = FearCloud_RSA('get_subjects');
% subjects        = Project.subjects_bdnf(Project.subjects_ET);
% subjects        = setdiff(subjects,58);
%%
for kernel_fwhm = 29.6%10:10:100;
    fwhm            = fwhm  + 1;
    %
    fix             = Fixmat(subjects,phase);%get the data
    fix.kernel_fwhm = kernel_fwhm;%set the smoothening factor
    fix.unitize     = 1;%unitize or not, this is important and I will make some tests with this.
    %% get number of trials per condition and subject: Sanity check...
    M               = [];%this will contain number of trials per subject and per condition. some few people have 10 trials (not 11) but that is ok. This is how I formed the subject_exlusion variable.
    sub_c           = 0;
    for ns = subjects(:)'
        sub_c = sub_c + 1;
        nc = 0;
        for cond = unique(fix.deltacsp)
            nc          = nc +1;
            i           = ismember(fix.subject,ns).*ismember(fix.deltacsp,cond);
            i           = logical(i);
            M(sub_c,nc) = length(unique(fix.trialid(i)));
        end
    end
    %% get all the single trials in a huge matrix D together with labels.
    global_counter = 0;
    clear D;%the giant data matrix
    clear labels;%and associated labels.
    for ns = subjects(:)'
        for deltacsp = -135:45:180;
            i              = (fix.subject == ns).*(fix.deltacsp == deltacsp);
            trials         = unique(fix.trialid(i == 1));
            if nargin >1
                switch varargin{1}
                    case '1'
                        trials = trials(ismember(trials,1:120));
                    case '2'
                        trials = trials(ismember(trials,121:240));
                    case '3'
                        trials = trials(ismember(trials,241:360));
                end
            end
            trial_counter  = 0;
            for trialid = trials
                trial_counter       = trial_counter + 1;
                global_counter      = global_counter +1;
                c                   = global_counter;
                v                   = {'subject' ns 'phase' phase 'deltacsp' deltacsp 'trialid' trialid};
                fix.getmaps(v);
                D(:,c)              = Vectorize(imresize(fix.maps,.1));
                labels.sub(c)       = ns;
                labels.phase(c)     = phase;
                labels.trial(c)     = trial_counter;%some people have less trials check it out with plot(labels.trial)
                labels.cond(c)      = deltacsp;
            end
        end
    end
    %% DATA2LOAD get the eigen decomposition: D is transformed to TRIALLOAD
    fprintf('starting covariance computation\n')
    covmat    = cov(D');
    fprintf('done\n')
    fprintf('starting eigenvector computation\n')
    [e dv]    = eig(covmat);
    fprintf('done\n')
    dv        = sort(diag(dv),'descend');
    eigval(fwhm,:) = dv;
%     plot(cumsum(dv)./sum(dv),'o-');xlim([0 200]);drawnow
    eigen     = fliplr(e);
    %collect loadings of every trial
    trialload = D'*eigen(:,1:teig)*diag(dv(1:teig))^-.5;%dewhitened
    %% LIBSVM business
    cev = 0;
    for neig    = 1:teig%test different numbers of eigenvectors
        cev = cev+1;
        sub_counter = 0;
        result      = [];
        w           = [];
        for sub = unique(labels.sub)%go subject by subject
            fprintf('Smooth:%d-Eig:%d-Sub:%d\n',fwhm,neig,sub);
            sub_counter = sub_counter + 1;
            ind_all     = ismember(labels.sub,sub);%this subject, this phase.
            %
            for n = 1:tbootstrap%
                Ycond   = double(labels.cond(ind_all))';%labels of the fixation maps for this subject in this phase.
                X       = trialload(ind_all,1:neig);%fixation maps of this subject in this phase.
                % this takes into account, that ntrials_base = 120 but
                % ntrials_test = 360
                % so first we chose a subset of 120 trials, then do the
                % usual holdout procedure on that.
                if sample120first == 1
                    P120 = cvpartition(Ycond,'Holdout',2/3); % we just take this to get 1/3 (11 rep) of the data, respecting conditions
                    %now replace Ycond and X, we just need this subsample
                    %of trials
                    Ycond   = Ycond(P120.training);%labels of the fixation maps for this subject in this phase.
                    X       = X(P120.training,:);%fixation maps of this subject in this phase.
                end
                % now normal Holdout for every phase (which should all have the
                % same number of trials now)
                P       = cvpartition(Ycond,'Holdout',holdout_ratio(phase)); % divide training and test datasets respecting conditions
                i       = logical(P.training.*ismember(Ycond,[0 180]));%train using only the CS+ and CS? conditions.
                if random ==1
                    warning('Randomizing labels as wanted! \n')
                    model   = svmtrain(Shuffle(Ycond(i)), X(i,1:neig), '-t 0 -c 1 -q'); %t 0: linear, -c 1: criterion, -q: quiet
                else
                    model   = svmtrain(Ycond(i), X(i,1:neig), '-t 0 -c 1 -q'); %t 0: linear, -c 1: criterion, -q: quiet
                end
                % get the hyperplane
                try
                    w(:,sub_counter,n)          = model.SVs'*model.sv_coef;
                catch
                    keyboard%sanity check: stop if something is wrong
                end
                %%
                cc=0;
                for cond = unique(Ycond)'
                    cc                          = cc+1;
                    i                           = logical(P.test.*ismember(Ycond,cond));%find all indices that were not used for training belonging to COND.
                    [~, dummy]                  = evalc('svmpredict(Ycond(i), X(i,:), model);');%doing it like this supresses outputs.
                    dummy                       = dummy == 0;%binarize: 1=CS+, 0=Not CS+
                    result(cc,n,sub_counter)    = sum(dummy)/length(dummy);%get the percentage of CS+ responses for each CONDITION,BOOTSTR,SUBJECT
                end
            end
        end
        %once the data is there compute 2 important metrics:
        R(:,cev,fwhm)      = mean(mean(result,2),3);%average across bootstaps the classification results
        AVEHP(:,:,cev,fwhm)  = reshape(mean(eigen(:,1:neig)*mean(w,3),2),[50 50 1]);%average hyperplane across subjects
        HP(:,:,:,cev,fwhm) = reshape(eigen(:,1:neig)*mean(w,3),[50 50 size(eigen(:,1:neig)*mean(w,3),2)]); %single hyperplanes
    end
%     figure(1000);imagesc(R(:,:,fwhm));
%     figure(1001);plot(R(4,:,fwhm)-R(end,:,fwhm),'o-');
%     drawnow
end
try
    save(sprintf('C:/Users/Lea/Documents/Experiments/project_bdnf/data/midlevel/svm_analysis/findingparams/SVM_simulations_phase%d_NEV1to%d_FWHM30_alltrials.mat',phase,neig),'R','AVEHP','HP','eigval')
    disp('saved successfully')
    keyboard
catch
    fprintf('couldnt save, do sth!')
    keyboard
end
