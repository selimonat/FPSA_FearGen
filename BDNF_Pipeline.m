%%BDNF Stack
addpath([homedir '/Documents/GitHub/fancycarp/'])
p = Project;
mask = p.getMask('ET_feargen');
%create groups. 1 = GG, 2 = AG,AA
g1 = Group(intersect(find(mask),Project.subjects_bdnf(Project.BDNF ==1)));
g2 = Group(intersect(find(mask),Project.subjects_bdnf(Project.BDNF ==2)));
N1 = length(g1.ids);
N2 = length(g2.ids);
%get SI
g1.getSI(3);
g2.getSI(3);
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
    bar(1,mean(n.sigma_cond));
    hold on;
    bar(2,mean(n.sigma_test));
    bar(3,mean(n.SI));
    errorbar(1,mean(n.sigma_cond),std(n.sigma_cond)./sqrt(length(n.ids)),'k.','linewidth',2)
    errorbar(2,mean(n.sigma_test),std(n.sigma_test)./sqrt(length(n.ids)),'k.','linewidth',2)
    errorbar(3,mean(n.SI),std(n.SI)./sqrt(length(n.ids)),'k.','linewidth',2)
    ylabel(sprintf('Fear tuning parameter (deg)'))
    t = title(['Group ' num2str(c)]);set(t,'FontSize',14);
    set(gca,'xtick',[1 2 3],'xticklabel',{'Cond' 'Test' 'SI'})
    hold off;
    box off;
end
EqualizeSubPlotYlim(gcf)
%% ttests
[H,P,CI,STATS] = ttest2(g1.sigma_cond,g2.sigma_cond)
[H,P,CI,STATS] = ttest2(g1.sigma_test,g2.sigma_test)
[H,P,CI,STATS] = ttest2([g1.sigma_test ;g1.sigma_cond],[g2.sigma_test ;g2.sigma_cond]);
[H,P,CI,STATS] = ttest2(g1.SI,g2.SI)
%% ANOVA
[p,tbl,stats] = anova1([g1.sigma_cond; g2.sigma_cond],[ones(length(g1.ids),1); ones(length(g2.ids),1)*2])
[p,tbl,stats] = anova1([g1.sigma_test; g2.sigma_test],[ones(length(g1.ids),1); ones(length(g2.ids),1)*2])
[p,tbl,stats] = anova1([g1.SI; g2.SI],[ones(length(g1.ids),1); ones(length(g2.ids),1)*2])
%% ANOVAN
data  = [g1.sigma_cond         ; g2.sigma_cond;g1.sigma_test; g2.sigma_test]
group = [[ones(length(g1.ids),1); ones(length(g2.ids),1)*2   ;ones(length(g1.ids),1); ones(length(g2.ids),1)*2 ] [ones(length(g1.ids)*2,1); ones(length(g2.ids)*2,1)*2  ]]
[p,tbl,stats] = anovan(data,group,'model','interaction');

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
