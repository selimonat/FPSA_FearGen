p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
fix             = Fixmat(subjects,[1 2 3 4 5]);
fix.getsubmaps
fix.maps   = imresize(fix.maps,.1,'method','bilinear');



%% all phases, for subject identification
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);

phases = 1:5;
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
v = [];
c=0;
for sub = subjects(:)'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr 'deltacsp' fix.realcond};
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

%cut the nans
todelete = isnan(labels.sub);
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end

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


%% Discrimination
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),subjects);
g               = Group(subjects);
[mat,tag] = g.parameterMat;
clear g

%good discriminator vs bad
crit = median(mean(mat(:,[1 3]),2));
sub_good = subjects(mean(mat(:,[1 3]),2)<=crit);

%prepare fixmaps
phases =1;
fix = Fixmat(subjects,phases);
% collect single trials
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
labels.discr     = NaN(1,ttrial);
v = [];
c=0;
deletecount = 0;
for sub = subjects(:)'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr 'deltacsp' fix.realcond};
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
                labels.discr(c)     = ismember(sub,sub_good);
                if ismember(ph,[1 5])
                    labels.pos(c)   = mod(tr,2);%1 is 1st, 0 is 2nd
                else
                    labels.pos(c)   = NaN;
                end
            else
                deletecount = deletecount +1;
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
labels.discr(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end

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


% save ~/Desktop/phase1.mat datamatrix labels
%% 
%% feargen 
clear all
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);
g.getSI(8);
[mat,tag] = g.parameterMat;
clear g;
subs_goodsigma = subjects(mat(:,13)  >=median(mat(:,13)));
subs_sharpener = subjects(mat(:,14)  >=median(mat(:,14)));
% 
%crit = log(mat(:,13)+1)-log(mat(:,12)+1);
% subs_sharpener = subjects(crit >= median(crit));
% subs_goodsigma = subjects(log(mat(:,13)+1) >=median(log(mat(:,13)+1)));

%prepare fixmaps
phases = 4;
fix = Fixmat(subjects,phases);
% fix.kernel_fwhm = 36;
% collect single trials
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
labels.SI        = NaN(1,ttrial);
labels.sigma     = NaN(1,ttrial);
v = [];
c=0;
for sub = subjects(:)'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr 'deltacsp' fix.realcond};
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
                labels.SI(c)        = ismember(sub,subs_sharpener);
                labels.sigma(c)     = ismember(sub,subs_goodsigma);
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
todelete = isnan(labels.sub);
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];
labels.SI(:,todelete)=[];
labels.sigma(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end


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
num =60;%min(find((cumsum(dv)./sum(dv))>.95));
%collect loadings of every trial
trialload = datamatrix'*eigen(:,1:num)*diag(dv(1:num))^-.5;%dewhitened



%% all phases, for subject identification
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);

phases = 1:5;
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
v = [];
c=0;
for sub = subjects(:)'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr 'deltacsp' fix.realcond};
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

%cut the nans
todelete = isnan(labels.sub);
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end

%% 
%% feargen4 - but connected to RATE AND PMF
clear all
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
% mask            = p.getMask('RATE');
% subjects        = intersect(find(mask),subjects);
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),subjects);
g               = Group(subjects);
g.getSI(3);
[mat,tag] = g.parameterMat;

subs_goodsigma = subjects(mat(:,13)  <= median(mat(:,13)));
subs_sharpener = subjects(mat(:,14)  >= median(mat(:,14)));
subs_goodalpha = subjects(mean(mat(:,[1 3]),2) <= median(mean(mat(:,[1 3]),2)));
%test: mean(mean(mat(find(~ismember(g.ids,subs_goodalpha)),[1 3]),2))

%prepare fixmaps
phases = 4;
fix = Fixmat(subjects,phases);
% fix.kernel_fwhm = 36;
% collect single trials
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
labels.SI        = NaN(1,ttrial);
labels.sigma     = NaN(1,ttrial);
labels.discr     = NaN(1,ttrial);
v = [];
c=0;
for sub = g.ids'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr 'deltacsp' fix.realcond};
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
                labels.SI(c)        = ismember(sub,subs_sharpener);
                labels.sigma(c)     = ismember(sub,subs_goodsigma);
                labels.discr(c)     = ismember(sub,subs_goodalpha);
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
todelete = isnan(labels.sub);
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];
labels.discr(:,todelete)=[];
labels.SI(:,todelete)=[];
labels.sigma(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end


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




%% this finds the number of trials per subjects in phase
for n=1:max(labels.easy_sub);trials(n,:)=hist(labels.phase(labels.easy_sub==n),unique(labels.phase(labels.easy_sub==n)));end
min(trials(:))
%% results SUBJECTS N=27
size(result,3)
a = mean(result,3);
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,27]);
mean(diag(scaled))
std(diag(scaled))
%% 2classes (alpha,kappa,SI)
size(result,3)
a = mean(result,3);
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,2]);
mean(diag(scaled))
std(diag(scaled))

%generalization/discrimination
a = result;
for n=1:100;scaled(:,:,n) = (a(:,:,n)./sum(sum(a(:,:,n))))./repmat(sum((a(:,:,n)./sum(sum(a(:,:,n)))),2),[1,2]);end
for n=1:100;correct(n,:) = diag(scaled(:,:,n));end
performance = mean(mean(correct,2));
sd_perf     = std(mean(correct,2));

%% hyperplanes
w0 = w;
w = mean(w0,2);
num = size(w,1);
hp = eigen(:,1:num)*w;
fix = Fixmat(6,3);
fix.getsubmaps;
fix.maps = reshape(hp,[50 50]);
fix.plot;


g = Group(unique(labels.sub));g.getSI(3);[mat tags] = g.parameterMat; clear g;
phase = 4;
param = mat(:,13);%mean(mat(:,[1 3]),2);

fix = Fixmat(unique(labels.sub),phase);fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmap    = fix.vectorize_maps;
hpload     = hp(:)'*subjmap;
[r,pval]=corr(hpload',param)

p = Project; mask = p.getMask('PMF');
valid = ismember(unique(labels.sub),find(sum(mask,2)==4));
[r,pval]=corr(hpload(valid)',mean(mat(valid,[1 3]),2))


%% classify csp vs csn
clear all
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
fix = Fixmat(subjects,4);
phases = 4;

deltacsp = [0 180];
% fix.kernel_fwhm = 36;
% collect single trials
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
labels.csp       = NaN(1,ttrial);

v = [];
c=0;
for cspcsn = deltacsp(:)'
    for sub = subjects(:)'
        for ph = phases
            for tr = 1:max(fix.trialid(fix.phase == ph))
                v = {'subject' sub, 'phase' ph 'trialid' tr 'deltacsp' cspcsn};
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
                    labels.csp(c)       = cspcsn==0;
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
todelete = isnan(labels.sub);
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];
labels.csp(:,todelete)=[];


c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end


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
% num = min(find((cumsum(dv)./sum(dv))>.95));
num=60;
%collect loadings of every trial
trialload = datamatrix'*eigen(:,1:num)*diag(dv(1:num))^-.5;%dewhitened


%% subject identificatino
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\versionABCDEFG\subjects_rand0\result.mat')
a = result;
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,27]);
subplot(1,2,1);imagesc(mean(scaled,3))
axis square
subplot(1,2,2);a = CancelDiagonals(mean(scaled,3),0);imagesc(a);
axis square

