%% Discrimination
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),subjects);
g               = Group(subjects);
g.getSI(3);
[mat,tag] = g.parameterMat;

%good discriminator vs bad
crit = median(mean(mat(:,[1 3]),2));
sub_good = subjects(mean(mat(:,[1 3]),2)<=crit);

%prepare fixmaps
phases =1;
fix = Fixmat(subjects,phases);
fix.kernel_fwhm = 36;
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
                labels.discr(c)     = ismember(sub,sub_good);
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
g.getSI(3);
[mat,tag] = g.parameterMat;

subs_goodsigma = subjects(mat(:,13)  <= median(mat(:,13)));
subs_sharpener = subjects(mat(:,14)  >= median(mat(:,14)));


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
num = min(find((cumsum(dv)./sum(dv))>.95));
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
%% results
size(result,3)
a = mean(result,3);
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,27]);
mean(diag(scaled))
std(diag(scaled))
% this finds the number of trials per subjects in phase
for n=1:max(labels.easy_sub);trials(n,:)=hist(labels.phase(labels.easy_sub==n),unique(labels.phase(labels.easy_sub==n)));end
min(trials(:))

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

g = Group(unique(labels.sub));g.getSI(3);[mat tags] = g.parameterMat; fix = Fixmat(g.ids,4);fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmap    = fix.vectorize_maps;
hpload     = hp(:)'*subjmap;
[r,pval]=corr(hpload',g.SI)
%correlating ph04 with something of ph01
p = Project;
mask = p.getMask('ET_discr');ET_people = find(mask);mask = p.getMask('PMF');discr_people = find(sum(mask,2)==4);
ETPMF = intersect(ET_people,discr_people);
ETPMF = intersect(ETPMF,Project.subjects_1500);
valid = ismember(g.ids,ETPMF);
sum(valid)
%or
p = Project;
mask = p.getMask('ET_feargen');ET_people = find(mask);mask = p.getMask('RATE');rate_people = find(mask);
ETRATE = intersect(ET_people,rate_people);
ETRATE = intersect(ETRATE,Project.subjects_1500);
valid = ismember(g.ids,ETRATE);
sum(valid)


