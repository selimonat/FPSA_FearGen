function [R AVE result HP]=SVM_simulations(phase)
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

subject_exlusion = [20 22 7];%these people do not have much trials so kick them out.
tbootstrap       = 100;%number of bootstraps
nfix              = 2;
% phase            = 2;%baseline or test
holdout_ratio    = [NaN .5 NaN .25];%see above
teig             = 100;%how many eigenvalues do you want to include for each run.
fwhm             = 0;%counter
R                = [];%result storage
AVE              = [];%result storate
HP               = [];
%%
for kernel_fwhm = 10:10:100;
    fwhm            = fwhm  + 1;
    %
    fix             = Fixmat(setdiff(Project.subjects_1500,7),phase);%get the data
    fix.kernel_fwhm = kernel_fwhm;%set the smoothening factor
    fix.unitize     = 1;%unitize or not, this is important and I will make some tests with this.
    %% get number of trials per condition and subject: Sanity check...
    M               = [];%this will contain number of trials per subject and per condition. some few people have 10 trials (not 11) but that is ok. This is how I formed the subject_exlusion variable. 
    sub_c           = 0;
    for ns =setdiff(Project.subjects_1500,subject_exlusion)
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
    for ns = setdiff(Project.subjects_1500,subject_exlusion)
        for deltacsp = -135:45:180;
            i              = (fix.subject == ns).*(fix.deltacsp == deltacsp);
            trials         = unique(fix.trialid(i == 1));
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
            fprintf('Eig:%d-Smooth:%d-Sub:%d\n',neig,fwhm,sub);
            sub_counter = sub_counter + 1;            
            ind_all     = ismember(labels.sub,sub);%this subject, this phase.                        
            %
            for n = 1:tbootstrap%
                Ycond   = double(labels.cond(ind_all))';%labels of the fixation maps for this subject in this phase.            
                X       = trialload(ind_all,1:neig);%fixation maps of this subject in this phase.
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
        AVE(:,:,cev,fwhm)  = reshape(mean(eigen(:,1:neig)*mean(w,3),2),[50 50 1]);%average the hyperplans
        HP = reshape(eigen(:,1:neig)*mean(w,3),[50 50 26]);
    end
    figure(1000);imagesc(R(:,:,fwhm));
    figure(1001);plot(R(4,:,fwhm)-R(end,:,fwhm),'o-');
    drawnow
end
try
    save(sprintf('C:/Users/Lea/Documents/Experiments/FearCloud_Eyelab/data/midlevel/svm_analysis/findingparams/SVM_simulations_phase%d_NEV%d_FWHM%d.mat',phase,neig,kernel_fwhm),'R','AVE','HP')
catch
    fprintf('couldnt save, do sth!')
    keyboard
end
