%%BDNF Stack
addpath([homedir '/Documents/GitHub/fancycarp/'])
p = Project;
% mask = p.getMask('ET_feargen');
%create groups. 1 = GG, 2 = AG,AA
g1 = Group(Project.subjects_bdnf(Project.BDNF ==1));
g2 = Group(Project.subjects_bdnf(Project.BDNF ==2));
N1 = length(g1.ids);
N2 = length(g2.ids);
%get SI
g1.getSI(8);
g2.getSI(8);
%
threshold = .05;
valid2 = prod([g2.tunings.rate{3}.pval;g2.tunings.rate{4}.pval] > -log10(threshold));
valid1 = prod([g1.tunings.rate{3}.pval;g1.tunings.rate{4}.pval] > -log10(threshold));
g1.sigma_cond(~valid1) = NaN;
g2.sigma_cond(~valid2) = NaN;
g1.sigma_test(~valid1) = NaN;
g2.sigma_test(~valid2) = NaN;
g1.SI(~valid1) = NaN;
g2.SI(~valid2) = NaN;

%%
g1.tunings.rate{3}.GroupFit(8)
g1.tunings.rate{4}.GroupFit(8)
g2.tunings.rate{3}.GroupFit(8)
g2.tunings.rate{4}.GroupFit(8)
%% plot the results as bars with fitted Gaussian
g1.PlotRatingResults
t = supertitle('Group 1');set(t,'FontSize',14);
g2.PlotRatingResults
t = supertitle('Group 2');set(t,'FontSize',14);

%% barplot sigma_cond|sigma_test|SI
c=0;
for n = [g1 g2]
    c=c+1;
    subplot(1,2,c)
    bar(1,nanmean(n.sigma_cond));
    hold on;
    bar(2,nanmean(n.sigma_test));
    bar(3,nanmean(n.SI));
    errorbar(1,nanmean(n.sigma_cond),nanstd(n.sigma_cond)./sqrt(sum(~isnan(n.sigma_cond))),'k.','linewidth',2)
    errorbar(2,nanmean(n.sigma_test),nanstd(n.sigma_test)./sqrt(sum(~isnan(n.sigma_test))),'k.','linewidth',2)
    errorbar(3,nanmean(n.SI),nanstd(n.SI)./sqrt(sum(~isnan(n.sigma_cond))),'k.','linewidth',2)
    t = title(['Group ' num2str(c)]);set(t,'FontSize',14);
    set(gca,'xtick',[1 2 3],'xticklabel',{'Cond' 'Test' 'SI'})
    ylabel(sprintf('Fear tuning parameter (deg)'))
    set(gca,'fontsize',14);
    hold off;
    box off;
end
EqualizeSubPlotYlim(gcf)
%% get mu
for s = 1:length(g1.ids)
    mu1(s,1) = g1.tunings.rate{3}.singlesubject{s}.Est(3);
    mu1(s,2) = g1.tunings.rate{4}.singlesubject{s}.Est(3);
end
for s = 1:length(g2.ids)
    mu2(s,1) = g2.tunings.rate{3}.singlesubject{s}.Est(3);
    mu2(s,2) = g2.tunings.rate{4}.singlesubject{s}.Est(3);
end
mat1 = [g1.sigma_cond g1.sigma_test g1.SI mu1];
mat2 = [g2.sigma_cond g2.sigma_test g2.SI mu2];
tags = {'kappa_cond','kappa_test','SI','Mu_cond','mu_test'};

%% ttests
[H,P,CI,STATS] = ttest2(g1.sigma_cond,g2.sigma_cond)
[H,P,CI,STATS] = ttest2(g1.sigma_test,g2.sigma_test)
[H,P,CI,STATS] = ttest2([g1.sigma_test ;g1.sigma_cond],[g2.sigma_test ;g2.sigma_cond]);
[H,P,CI,STATS] = ttest2(g1.SI,g2.SI)

[H,P,CI,STATS] = ttest2(abs(mat1(:,end-1)),abs(mat2(:,end-1)))%mu cond
[H,P,CI,STATS] = ttest2(abs(mat1(:,end)),abs(mat2(:,end)))%mu test

%% ANOVA
[p,tbl,stats] = anova1([g1.sigma_cond; g2.sigma_cond],[ones(length(g1.ids),1); ones(length(g2.ids),1)*2])
[p,tbl,stats] = anova1([g1.sigma_test; g2.sigma_test],[ones(length(g1.ids),1); ones(length(g2.ids),1)*2])
[p,tbl,stats] = anova1([g1.SI; g2.SI],[ones(length(g1.ids),1); ones(length(g2.ids),1)*2])
%% ANOVAN
data  = [g1.sigma_cond; g2.sigma_cond; g1.sigma_test; g2.sigma_test];
group = [[ones(length(g1.ids),1); ones(length(g2.ids),1)*2   ;ones(length(g1.ids),1); ones(length(g2.ids),1)*2 ] [ones(length(g1.ids)*2,1); ones(length(g2.ids)*2,1)*2  ]];
[p,tbl,stats] = anovan(data,group,'model','interaction','varnames',char('group', 'timepoint'));

%% fixation maps
fix  = Fixmat(10,2);%this is a dummy to be filled later.
fix1 = Fixmat(g1.ids,4);
fix2 = Fixmat(g2.ids,4);
%get maps for g1 and g2, male and female
%g1
v = [];
c = 0;
for n = 1:2;
    c    = c+1;
    v{c} = {'subject' intersect(g1.ids,find(Project.gender == n))};
end
fix1.getmaps(v{:});
%g2
v = [];
c = 0;
for n = 1:2;
    c    = c+1;
    v{c} = {'subject' intersect(g2.ids,find(Project.gender == n))};
end
fix2.getmaps(v{:});
fix.maps = cat(3,fix1.maps,fix2.maps);
fix.plot;
%cocktailblanked
fix.maps = fix.maps - repmat(mean(fix.maps(:,:,1:4),3),[1 1 4]);
fix.plot;

%% corr matrix 8x8 conditions
subplot(2,2,1);imagesc(median(fix1.condcorrmat(2,0),3));ylabel('Group 1');title('Baseline')
subplot(2,2,2);imagesc(median(fix1.condcorrmat(4,0),3));title('Testphase')
subplot(2,2,3);imagesc(median(fix2.condcorrmat(2,0),3));ylabel('Group 2');title('Baseline')
subplot(2,2,4);imagesc(median(fix2.condcorrmat(4,0),3));title('Testphase')
%%
%gender differences
g = Group(Project.subjects_BDNF);
g.getSI(3);
Y = g.SI;
genes  = Project.BDNF(g.ids);
gender = Project.gender(subjects);
[p,tbl,stats,terms] = anovan(Y,{genes gender}, 'model','interaction','varnames',{'BDNF','Gender'});
%% classify groups in testphase
% prepare data
clear all
p               = Project;
subjects        = find(p.getMask('ET_feargen'));
fix = Fixmat(subjects,4);
phases = 4;

deltacsp = [0 180];
% fix.kernel_fwhm = 36;
% collect single trials
scale            = .1;%use .1 to have the typical 2500 pixel maps
final_size       = prod(fix.rect(end-1:end).*scale);
trialnumber      = ['fill me in'];
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
num = min(find((cumsum(dv)./sum(dv))>.95));
%collect loadings of every trial
trialload = datamatrix'*eigen(:,1:num)*diag(dv(1:num))^-.5;%dewhitened

