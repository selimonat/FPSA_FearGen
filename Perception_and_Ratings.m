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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot all alphas for N=51, 1500 and 600 seperately
bar([1 2 4 5 7 8 10 11 13 14 16 17],[mean(mat(:,1)) mean(mat(:,2)) mean(mat(:,3)) mean(mat(:,4)) mean(mat15(:,1)) mean(mat15(:,2)) mean(mat15(:,3)) mean(mat15(:,4)) mean(mat6(:,1)) mean(mat6(:,2)) mean(mat6(:,3)) mean(mat6(:,4))],'FaceColor','white')
hold on;
errorbar([1 2 4 5 7 8 10 11 13 14 16 17],[mean(mat(:,1)) mean(mat(:,2)) mean(mat(:,3)) mean(mat(:,4)) mean(mat15(:,1)) mean(mat15(:,2)) mean(mat15(:,3)) mean(mat15(:,4)) mean(mat6(:,1)) mean(mat6(:,2)) mean(mat6(:,3)) mean(mat6(:,4))],[std(mat(:,1))./sqrt(51) std(mat(:,2))./sqrt(51) std(mat(:,3))./sqrt(51) std(mat(:,4))./sqrt(51) std(mat15(:,1))./sqrt(25) std(mat15(:,2))./sqrt(25) std(mat15(:,3))./sqrt(25) std(mat15(:,4))./sqrt(25) std(mat6(:,1))./sqrt(26) std(mat6(:,2))./sqrt(26) std(mat6(:,3))./sqrt(26) std(mat6(:,4))./sqrt(26)],'k.')
set(gca,'xtick',[1.5 4.5 8.5 11.5 15.5 18.5],'xticklabel',{'CS+' 'CS-' 'CS+' 'CS-' 'CS+' 'CS-'})

%% plot actual figure
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
load('C:\Users\onat\Desktop\Lea\behavioral_data.mat','data')
[h,p,ci,stats] = ttest(data(:,1),data(:,3))%before the experiment, CSP vs CSN

[h,p,ci,stats] = ttest(data(:,1),data(:,2))%CSP before to after
[h,p,ci,stats] = ttest(data(:,3),data(:,4))%CSN before to after

[h,p,ci,stats] = ttest(data(:,9),data(:,10),'tail','right') % csp_imprv vs csn_imprv
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p    = Project;
mask = p.getMask('RATE');
Subjects = intersect(find(mask),[Project.subjects_1500 Project.subjects_600]);
g = Group(Subjects);
g.getSI(3);
[mat tags] = g.parameterMat;

figure
bar([1 2 3 6 7 8 10 11 12],[mean(mat(:,12)) mean(mat(:,13)) mean(mat(:,14)) mean(mat(1:18,12)) mean(mat(1:18,13)) mean(mat(1:18,14)) mean(mat(19:end,12)) mean(mat(19:end,13)) mean(mat(19:end,14))],'facecolor','white')
hold on;
errorbar([1 2 3 6 7 8 10 11 12],[mean(mat(:,12)) mean(mat(:,13)) mean(mat(:,14)) mean(mat(1:18,12)) mean(mat(1:18,13)) mean(mat(1:18,14)) mean(mat(19:end,12)) mean(mat(19:end,13)) mean(mat(19:end,14))],...
    [std(mat(:,12))./sqrt(37) std(mat(:,13))./sqrt(37) std(mat(:,14))./sqrt(37) std(mat(1:18,12))./sqrt(18) std(mat(1:18,13))./sqrt(18) std(mat(1:18,14))./sqrt(18) std(mat(19:end,12))./sqrt(19) std(mat(19:end,13))./sqrt(19) std(mat(19:end,14))./sqrt(19)],'k.')
set(gca,'xticklabel',{'Cond_all' 'Test_all' 'SI_all' 'Cond_1500' 'Test_1500' 'SI_1500' 'Cond_600' 'Test_600' 'SI_600'})

figure;
bar([1 2 3],mean(mat(:,12:14)));hold on;
errorbar([1 2 3],mean(mat(:,12:14)),std(mat(:,12:14))./sqrt(length(mat)),'k.','LineWidth',2)

[h,p,ci,stats] = ttest(data(:,12),data(:,13))

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
g15.getSI(3)
mat15 = g15.parameterMat;
clear g15;
p               = Project;
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),Project.subjects_600);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g6               = Group(subjects);
g6.getSI(3)
mat6 = g6.parameterMat;
clear g6;
mat = [mat15;mat6];
%% corr alpha x feargen
[r,p]=corrcoef(mean(mat(:,[1 3]),2),mat(:,12))
[r,p]=corrcoef(mean(mat(:,[1 3]),2),mat(:,13))
[r,p]=corrcoef(mean(mat(:,[1 3]),2),mat(:,14))


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
xlabel(sprintf(sprintf('sigma_{cond} - sigma_{test} (in deg) \n Sharpening Index')))
axis square
set(gca,'fontsize',20)

%% all correlations
clf
mat(:,15) = mean(mat(:,[1 3]),2);
mat(:,16) = mean(mat(:,[2 4]),2);
tags{15} = 'alpha_before';
tags{16} = 'alpha_after';
[r,pval]=corr(mat);
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

