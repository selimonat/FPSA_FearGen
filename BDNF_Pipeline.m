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
[mat1 tags] = g1.misesMat;
mat2 = g2.misesMat;
threshold = .05;
valid2 = prod([g2.tunings.rate{3}.pval;g2.tunings.rate{4}.pval] > -log10(threshold));
valid1 = prod([g1.tunings.rate{3}.pval;g1.tunings.rate{4}.pval] > -log10(threshold));
g1.sigma_cond(~valid1) = NaN;
g2.sigma_cond(~valid2) = NaN;
g1.sigma_test(~valid1) = NaN;
g2.sigma_test(~valid2) = NaN;
g1.SI(~valid1) = NaN;
g2.SI(~valid2) = NaN;
mat1(~valid1,:) = NaN;
mat2(~valid2,:) = NaN;

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
%% barplot Kappa|SI|Mu
figure;
subplot(1,2,1)
bar([1 2 4 5 7 8],[nanmean(mat1(:,1)),nanmean(mat2(:,1)),nanmean(mat1(:,2)),nanmean(mat2(:,2)),nanmean(mat1(:,3)),nanmean(mat2(:,3))])
box off
hold on;errorbar([1 2 4 5 7 8],...
    [nanmean(mat1(:,1)),nanmean(mat2(:,1)),nanmean(mat1(:,2)),nanmean(mat2(:,2)),nanmean(mat1(:,3)),nanmean(mat2(:,3))],...
    [nanstd(mat1(:,1))./sqrt(sum(valid1)),nanstd(mat2(:,1))./sqrt(sum(valid2)),nanstd(mat1(:,2))./sqrt(sum(valid1)),nanstd(mat2(:,2))./sqrt(sum(valid2)),nanstd(mat1(:,3))./sqrt(sum(valid1)),nanstd(mat2(:,3))./sqrt(sum(valid2))],'k.','LineWidth',1.5)
legend('Group 1','Group 2')
set(gca,'XTick',[ 1.5 4.5 7.5],'XTickLabel',{'kappa_cond' 'kappa_test' 'SI'})
legend('Group 1')
subplot(1,2,2)
bar([1 2 4 5],[nanmean(abs(mat1(:,4))),nanmean(abs(mat2(:,4))),nanmean(abs(mat1(:,5))),nanmean(abs(mat2(:,5)))])
box off
hold on;errorbar([1 2 4 5],...
    [nanmean(abs(mat1(:,4))),nanmean(abs(mat2(:,4))),nanmean(abs(mat1(:,5))),nanmean(abs(mat2(:,5)))],...
    [nanstd(abs(mat1(:,4)))./sqrt(sum(valid1)),nanstd(abs(mat2(:,4)))./sqrt(sum(valid2)),nanstd(abs(mat1(:,5)))./sqrt(sum(valid1)),nanstd(abs(mat2(:,5)))./sqrt(sum(valid2))],'k.','LineWidth',1.5)
legend('Group 1','Group 2')
set(gca,'XTick',[ 1.5 4.5],'XTickLabel',{'Mu_cond' 'Mu_test'},'YTick',0:10:45)
legend('Group 1')

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
triallnoad = datamatrix'*eigen(:,1:num)*diag(dv(1:num))^-.5;%dewhitened

%% SCR
p = Project;
phasenames = {'base$' 'cond$' 'test$'};
phaseconds = [9 3 9];
bdnf_problems = [2 13 18 31 32 39 50];
subjects = Project.subjects_bdnf(~ismember(Project.subjects_bdnf,bdnf_problems));
for ph = 1:3
    sc = 0;
    for sub = subjects(:)'
        fprintf('Working on phase %d, subject %d .. \n',ph,sub)
        sc=sc+1;
        s = Subject(sub);
        ids(sc) = s.id;
        s.scr.cut(s.scr.findphase(phasenames{ph}));
        s.scr.run_ledalab;
        s.scr.plot_tuning_ledalab(1:phaseconds(ph))
        if ismember(ph,[1 3])
            scr_bars(sc,:,ph) = s.scr.fear_tuning;
        else
            scr_bars(sc,[4 8 9],ph) = s.scr.fear_tuning;
        end
        close all
        clear s
    end
end
subs1 = Project.subjects_bdnf(Project.BDNF ==1);
subs2 = Project.subjects_bdnf(Project.BDNF ==2);
members(:,1) = ismember(subjects,subs1);
members(:,2) = ismember(subjects,subs2);
%% OR
load('C:\Users\user\Google Drive\EthnoMaster\BDNF\midlevel\scr_N70_BCT.mat');
%% plot
ylims = [0 1];
Yticks = 0:0.2:1;
for g = [1 2]
    for n = [1 3]
        subpl((-1+g)*3+n) = subplot(2,3,(-1+g)*3+n);
        b = bar(1:8,mean(scr_bars(members(:,g),1:8,n)));
        hold on;
        ylim(ylims)
        xlim([0 9])
        e = errorbar(1:8,mean(scr_bars(members(:,g),1:8,n)),std(scr_bars(members(:,g),1:8,n))./sqrt(size(scr_bars,1)),'k.');
        set(gca,'XTick',1:8,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',Yticks)
        SetFearGenBarColors(b)
        set(e,'LineWidth',2,'Color','k')
        axis square
        box off
    end
    subplot(2,3,(-1+g)*3+1);ylabel('SCR [muS]')
    %cond special treatment
    n=2;
    subpl((-1+g)*3+n) = subplot(2,3,(-1+g)*3+n);
    b(1) = bar(4,mean(scr_bars(members(:,g),4,2)));
    hold on;
    e(1) = errorbar(4,mean(scr_bars(members(:,g),4,2)),std(scr_bars(members(:,g),4,2))./sqrt(size(scr_bars,1)),'k.');
    ylim(ylims)
    xlim([0 9])
    set(gca,'XTick',1:8,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',Yticks);
    b(2) =  bar(8,mean(scr_bars(members(:,g),8,2)));
    e(2) = errorbar(8,mean(scr_bars(members(:,g),8,2)),std(scr_bars(members(:,g),8,2))./sqrt(size(scr_bars,1)));
    set(b(1),'FaceColor',[1 0 0],'EdgeColor','w');
    set(b(2),'FaceColor','c','EdgeColor','w');
    set(e,'LineWidth',2,'Color','k');
    axis square
    box off
    try
        for n = 1:size(scr_bars,3)
            subplot(2,3,(-1+g)*3+n)
            line([0 9],repmat(mean(scr_bars(members(:,g),9,n)),[2 1]),'Color','k','LineWidth',1.3,'LineStyle',':')
        end
    end
end
subplot(2,3,2);title('Group 1')
subplot(2,3,5);title('Group 2')