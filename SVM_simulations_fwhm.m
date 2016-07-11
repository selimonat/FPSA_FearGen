function [R,Rbin, AVE,result,confmats]=SVM_simulations_fwhm(subjects,param)
%%
%This script is about SVM for sharp vs wide fear generalization. It collects the data and computes the eigenvalues on the
%fly for different sorts of parameters (kernel_fwhm, number of
%eigenvalues). This script will run train-test cycles for different smootheness levels and
%with different numbers of eigenvectors (up to TEIG).
% subject_exlusion = [20 22 7];%these people do not have much trials so kick them out.

addpath('/home/onat/Documents/Code/Matlab/edfread/current');
tbootstrap       = 5;%number of bootstraps
% phase            = 2;%baseline or test
holdout_ratio    = .2;%see above
teig             = 100;%how many eigenvalues do you want to include for each run.
fwhm             = 0;%counter
scale            = .1;
R                = [];%result storage
Rbin             = [];
AVE              = [];%result storate
p = Project;
%%
for kernel_fwhm = 10:10:80;
    fwhm            = fwhm  + 1;
    %
    fprintf('Time: %s \n',datestr(now,'HH:MM:SS'))
    fprintf('Collecting fixmat for kernel_fwhm = %g... \n',kernel_fwhm)
    [~,fix]             = evalc('Fixmat(subjects,4)');%get the data
    fix.kernel_fwhm = kernel_fwhm;%set the smoothening factor
    fix.unitize     = 1;%unitize or not, this is important and I will make some tests with this.
    final_size      = 50*50;
    
    %% get all the single trials in a huge matrix D together with labels.
    ttrial           = length(subjects)*240;
    D                = NaN(final_size,ttrial);
    labels.sub       = NaN(1,ttrial);
    labels.trial     = NaN(1,ttrial);
    labels.fwhm     = NaN(1,ttrial);
    v = [];
    c=0;
    fprintf('single trials for sub ')
    for sub = subjects(:)'
         fprintf('%d - ',sub);
        for tr = 1:max(fix.trialid)
            v = {'subject' sub 'trialid' tr 'deltacsp' fix.realcond};
            fix.getmaps(v);
            if ~any(isnan(fix.maps(:)))
                c                   = c+1;
                %scale it if necessary
                if scale ~= 1
                    fix.maps        = imresize(fix.maps,scale,'method','bilinear');
                end
                D(:,c)              = fix.vectorize_maps;
                labels.sub(c)       = sub;
                labels.trial(c)     = tr;
                labels.fwhm(c)     = param(subjects==sub) < median(param);
            end
        end
    end
   
%cut the nans
todelete = isnan(labels.sub);
fprintf('\nWill delete %g trials...\n',sum(todelete));
D(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.fwhm(:,todelete)=[]; 

%% DATA2LOAD get the eigen decomposition: D is transformed to TRIALLOAD
fprintf('starting covariance computation... ')
covmat    = cov(D');
fprintf('done. ')
fprintf('starting eigenvector computation... ')
[e dv]    = eig(covmat);
fprintf('done.\n')
dv        = sort(diag(dv),'descend');
%     plot(cumsum(dv)./sum(dv),'o-');xlim([0 200]);drawnow
eigen     = fliplr(e);
%collect loadings of every trial
trialload = D'*eigen(:,1:teig)*diag(dv(1:teig))^-.5;%dewhitened
%% LIBSVM business

for neig    = 1:teig%test different numbers of eigenvectors
    fprintf('Time: %s Eig:%d - Smooth:%d \n',datestr(now,'HH:MM:SS'),neig,fwhm);%this subject, this phase.
    %
    result      = [];
    w           = [];
    confmats    = [];
    for n = 1:tbootstrap%
        
        Y       = double(labels.fwhm)';%labels of the fixation maps for this subject in this phase.
        X       = trialload(:,1:neig);%fixation maps of this subject in this phase.
        P       = cvpartition(Y,'Holdout',holdout_ratio); %divide training and test datasets respecting conditions
        i       = logical(P.training);%train using only the CS+ and CS? conditions.
        nclass  = hist(Y(i),[0 1]);
        cmq     = sprintf('-t 0 -c 1 -q -w0 1 -w1 %g',nclass(1)/nclass(2));
        model   = svmtrain(Y(i), X(i,1:neig), cmq); %t 0: linear, -c 1: criterion, -q: quiet
        % get the hyperplane
        try
            w(:,n)          = model.SVs'*model.sv_coef;
        catch
            keyboard%sanity check: stop if something is wrong
        end
        %%
        
        i                           = logical(P.test);%.*ismember(Ycond,cond));%find all indices that were not used for training belonging to COND.
        [~,predicted]               = evalc('svmpredict(Y(i), X(i,:), model);');%doing it like this supresses outputs.
        confmats(:,:,n,neig,fwhm)   = confusionmat(Y(i),predicted,'order',[1 0]);
        result(n,neig,fwhm)         = sum(Y(i)==predicted)./length(predicted);
    end
    a = squeeze(mean(confmats(:,:,:,neig,fwhm),3));%mean confusionmat for this configuration (across bootstrps)
    %once the data is there compute 2 important metrics:
    R(:,neig,fwhm)      = mean(result(:,neig,fwhm),1);%average across bootstaps the classification results
    Rbin(:,neig,fwhm)   = diag((a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,2]));
    AVE(:,:,neig,fwhm)  = reshape(eigen(:,1:neig)*mean(w,2),[50 50]);%average the hyperplanes
end
end
try
     save(sprintf('%shome%skampermann%sDocuments%sfearcloud%sdata%smidlevel%ssvm_analysis%sfindingparams%sSVM_simulations_fwhm_NEV_%g_FWHM_%g_%s.mat',filesep,filesep,filesep,filesep,filesep,filesep,filesep,filesep,filesep,neig,fwhm,datestr(now,'ddmmyy'),'result','confmats','R','Rbin','AVE','subjects','param'));
catch
    keyboard
end
