%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p    = Project;
mask = p.getMask('PMF');
Subjects = intersect(find(sum(mask,2)==4),[Project.subjects_1500]);
g15= Group(Subjects);
mat15 = g15.parameterMat;
clear g15%for memory reasons
Subjects = intersect(find(sum(mask,2)==4),Project.subjects_600);
g6 = Group(Subjects);
mat6 = g6.parameterMat;
clear g6
mat = [mat15;mat6];
%%
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\pmf_data_N51.mat','mat','tags','mat15','mat6')
%% plot all alphas for N=51, 1500 and 600 seperately
bar([1 2 4 5 8 9 11 12 15 16 18 19],[mean(mat(:,1)) mean(mat(:,2)) mean(mat(:,3)) mean(mat(:,4)) mean(mat15(:,1)) mean(mat15(:,2)) mean(mat15(:,3)) mean(mat15(:,4)) mean(mat6(:,1)) mean(mat6(:,2)) mean(mat6(:,3)) mean(mat6(:,4))],'FaceColor','white')
hold on;
errorbar([1 2 4 5 8 9 11 12 15 16 18 19],[mean(mat(:,1)) mean(mat(:,2)) mean(mat(:,3)) mean(mat(:,4)) mean(mat15(:,1)) mean(mat15(:,2)) mean(mat15(:,3)) mean(mat15(:,4)) mean(mat6(:,1)) mean(mat6(:,2)) mean(mat6(:,3)) mean(mat6(:,4))],[std(mat(:,1))./sqrt(51) std(mat(:,2))./sqrt(51) std(mat(:,3))./sqrt(51) std(mat(:,4))./sqrt(51) std(mat15(:,1))./sqrt(25) std(mat15(:,2))./sqrt(25) std(mat15(:,3))./sqrt(25) std(mat15(:,4))./sqrt(25) std(mat6(:,1))./sqrt(26) std(mat6(:,2))./sqrt(26) std(mat6(:,3))./sqrt(26) std(mat6(:,4))./sqrt(26)],'k.')
set(gca,'xtick',[1.5 4.5 8.5 11.5 15.5 18.5],'xticklabel',{'CS+' 'CS-' 'CS+' 'CS-' 'CS+' 'CS-'})

%% plot differences for CSP and CSN
clf
bar([1 2 4 5 7 8],[mean(mat(:,1)-mat(:,2)) mean(mat(:,3)-mat(:,4)) mean(mat15(:,1)-mat15(:,2)) mean(mat15(:,3)-mat15(:,4)) mean(mat6(:,1)-mat6(:,2)) mean(mat6(:,3)-mat6(:,4))],'FaceColor','white')
hold on;
errorbar([1 2 4 5 7 8],[mean(mat(:,1)-mat(:,2)) mean(mat(:,3)-mat(:,4)) mean(mat15(:,1)-mat15(:,2)) mean(mat15(:,3)-mat15(:,4)) mean(mat6(:,1)-mat6(:,2)) mean(mat6(:,3)-mat6(:,4))],...
   [std(mat(:,1)-mat(:,2))./sqrt(length(mat)) std(mat(:,3)-mat(:,4))./sqrt(length(mat)) std(mat15(:,1)-mat15(:,2))./sqrt(length(mat15)) std(mat15(:,3)-mat15(:,4))./sqrt(length(mat15)) std(mat6(:,1)-mat6(:,2))./sqrt(length(mat6)) std(mat6(:,3)-mat6(:,4))./sqrt(length(mat6))],'k.')
set(gca,'xtick',[1.5 4.5 7.5],'xticklabel',{'pooled' '1500ms' '600ms'})

%% plot resnik figure
%better worse
csp_better = length(find(mat(:,1)-mat(:,3)>0))./length(mat(:,1));
csn_better = length(find(mat(:,2)-mat(:,4)>0))./length(mat(:,2));
bar([csp_better 1-csp_better; csn_better 1-csn_better],0.4,'stacked')
set(gca,'XTicklabel',{'CS+','CS-'})
xlim([0.5 2.5])
n1 = length(find(mat(:,1)-mat(:,3)>0)); N1 = length(mat(:,1));
n2 = length(find(mat(:,2)-mat(:,4)>0)); N2 = length(mat(:,2));
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [ones(n1,1); repmat(2,N1-n1,1); ones(n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

%%
[h,p,ci,stats] = ttest(mat(:,1),mat(:,2))%csp before to after
[h,p,ci,stats] = ttest(mat(:,3),mat(:,4))%csn before to after
[h,p,ci,stats] = ttest(mat(:,1)-mat(:,2),mat(:,3)-mat(:,4))%csp improvement different than csn improvement?

[h,p,ci,stats] = ttest(mat15(:,1),mat15(:,2))%csp before to after
[h,p,ci,stats] = ttest(mat15(:,3),mat15(:,4))%csn before to after
[h,p,ci,stats] = ttest(mat15(:,1)-mat15(:,2),mat15(:,3)-mat15(:,4))%csp improvement diff than csn improvement

[h,p,ci,stats] = ttest(mat6(:,1),mat6(:,2))%csp before to after
[h,p,ci,stats] = ttest(mat6(:,3),mat6(:,4))%csn before to after
[h,p,ci,stats] = ttest(mat6(:,1)-mat6(:,2),mat6(:,3)-mat6(:,4))%csp improvement diff than csn improvement
%% 
clf
[r,pval]=corr(mat15);
hold off
imagesc(r)
set(gca,'XTick',1:length(mat15),'YTick',1:length(mat15),'YTickLabel',tags)
axis image
colorbar
hold on;
[y x]=ind2sub(size(r),find(pval<0.05));
text(x,y,'x','FontSize',8)
[y x]=ind2sub(size(r),find(pval<0.01));
text(x,y,'x','FontSize',11)
[y x]=ind2sub(size(r),find(pval<0.001));
text(x,y,'x','FontSize',15)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p    = Project;
mask = p.getMask('RATEgauss');
Subjects = intersect(find(mask),Project.subjects_1500);
g15 = Group(Subjects);
g15.getSI(3);
[mat15 tags] = g15.parameterMat;
clear g15
Subjects = intersect(find(mask),Project.subjects_600);
g6 = Group(Subjects);
g6.getSI(3);
[mat6 tags] = g6.parameterMat;
clear g6
mat = [mat15;mat6];



figure
bar([1 2 3 6 7 8 11 12 13],[mean(mat(:,12)) mean(mat(:,13)) mean(mat(:,14)) mean(mat15(:,12)) mean(mat15(:,13)) mean(mat15(:,14)) mean(mat6(:,12)) mean(mat6(:,13)) mean(mat6(:,14))],'facecolor','white')
hold on;
errorbar([1 2 3 6 7 8 11 12 13],[mean(mat(:,12)) mean(mat(:,13)) mean(mat(:,14)) mean(mat15(:,12)) mean(mat15(:,13)) mean(mat15(:,14)) mean(mat6(:,12)) mean(mat6(:,13)) mean(mat6(:,14))],...
    [std(mat(:,12))./sqrt(length(mat)) std(mat(:,13))./sqrt(length(mat)) std(mat(:,14))./sqrt(length(mat)) std(mat15(:,12))./sqrt(length(mat15)) std(mat15(:,13))./sqrt(length(mat15)) std(mat15(:,14))./sqrt(length(mat15))...
    mean(mat6(:,12))./sqrt(length(mat6)) mean(mat6(:,13))./sqrt(length(mat6)) mean(mat6(:,14))./sqrt(length(mat6))],'k.')
set(gca,'xticklabel',{'Cond' 'Test' 'SI' 'Cond' 'Test' 'SI' 'Cond' 'Test' 'SI'})

[h,p,ci,stats] = ttest(mat(:,12),mat(:,13))
[h,p,ci,stats] = ttest(mat15(:,12),mat15(:,13))
[h,p,ci,stats] = ttest(mat6(:,12),mat6(:,13))

% plot bars and Gaussians in Feargen Colors
%
g.tunings.rate{3}.GroupFit(3)
g.tunings.rate{4}.GroupFit(3)
clf
subplot(1,2,1);h = bar(unique(g.tunings.rate{3}.x(1,:)),g.tunings.rate{3}.y_mean);SetFearGenBarColors(h);
axis square
hold on;
errorbar(unique(g.tunings.rate{3}.x(1,:)),g.tunings.rate{3}.y_mean,g.tunings.rate{3}.y_std./sqrt(27),'k.','LineWidth',2)
xlim([-160 200])
box off
set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'})
x = linspace(g.tunings.rate{3}.x(1,1),g.tunings.rate{3}.x(1,end),100);
% plot(x ,  g.tunings.rate{3}.singlesubject{1}.fitfun( x,mean(g.tunings.rate{3}.params(:,1:2))) ,'k--','linewidth',1);
plot(x ,  g.tunings.rate{3}.groupfit.fitfun( x,g.tunings.rate{3}.groupfit.Est(1:2)) ,'k--','linewidth',2);
hold off
set(gca,'fontsize',15)
%
subplot(1,2,2);h = bar(unique(g.tunings.rate{4}.x(1,:)),g.tunings.rate{4}.y_mean);SetFearGenBarColors(h);hold on;
errorbar(unique(g.tunings.rate{4}.x(1,:)),g.tunings.rate{4}.y_mean,g.tunings.rate{4}.y_std./sqrt(27),'k.','LineWidth',2)
EqualizeSubPlotYlim(gcf)
axis square
box off
xlim([-160 200])
set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'})
x = linspace(g.tunings.rate{4}.x(1,1),g.tunings.rate{4}.x(1,end),100);
plot(x ,  g.tunings.rate{4}.groupfit.fitfun( x,g.tunings.rate{4}.groupfit.Est(1:2)) ,'k-','linewidth',2);
% plot(x ,  g.tunings.rate{3}.singlesubject{1}.fitfun( x,mean(g.tunings.rate{4}.params(:,1:2))) ,'k--','linewidth',1);

set(gca,'fontsize',15)
hold off
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha x feargen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p               = Project;
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),Project.subjects_1500);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g15               = Group(subjects);
g15.getSI(8)
[mat15, tags] = g15.parameterMat;
clear g15;
p               = Project;
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),Project.subjects_600);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g6               = Group(subjects);
g6.getSI(8)
mat6 = g6.parameterMat;
clear g6;
mat = [mat15;mat6];

mat(:,end+1) = mean(mat(:,[1 3]),2);
mat15(:,end+1) = mean(mat15(:,[1 3]),2);
mat6(:,end+1) = mean(mat6(:,[1 3]),2);
tags{end+1} = 'initial_alpha'
%%
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat','mat','tags','mat15','mat6')
%% corr alpha x feargen
[r,pval]=corrcoef(mean(mat(:,[1 3]),2),mat(:,12))
[r,pval]=corrcoef(mean(mat(:,[1 3]),2),mat(:,13))
[r,pval]=corrcoef(mean(mat(:,[1 3]),2),mat(:,14))


% plot correlations between threshold and sigma/SI
subplot(1,3,1)
plot(mat(:,12),mean(mat(:,[1 3]),2),'k.','Markersize',25);box off;lsline
title(sprintf('r = %1.2g',corr2(mean(mat(:,[1 3]),2),mat(:,12))))
xlabel(sprintf('sigma_{cond} (in deg)'))
axis square
set(gca,'fontsize',20)
ylabel(sprintf('Perceptual Discrimination \n (alpha in deg)'))

subplot(1,3,2)
plot(mat(:,13),mean(mat(:,[1 3]),2),'k.','Markersize',25);box off;lsline
title(sprintf('r = %1.2g',corr2(mean(mat(:,[1 3]),2),mat(:,13))))
xlabel(sprintf('sigma_{test} (in deg)'));
axis square
set(gca,'fontsize',20)

subplot(1,3,3)
plot(mat(:,14),mean(mat(:,[1 3]),2),'k.','Markersize',25);box off;lsline
title(sprintf('r = %1.2g',corr2(mean(mat(:,[1 3]),2),mat(:,14))))
xlabel(sprintf(sprintf('Sharpening Index')))
axis square
set(gca,'fontsize',20)

%% 
clf
[r,pval]=corr(matbig);
hold off
imagesc(r)
set(gca,'XTick',1:length(mat),'YTick',1:length(mat),'YTickLabel',tags)
axis image
colorbar
hold on;
[y x]=ind2sub(size(r),find(pval<0.05));
text(x,y,'x','FontSize',8)
[y x]=ind2sub(size(r),find(pval<0.01));
text(x,y,'x','FontSize',11)
[y x]=ind2sub(size(r),find(pval<0.001));
text(x,y,'x','FontSize',15)
VerticalXlabel(tags,'interpreter','none')



