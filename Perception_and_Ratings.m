%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conditioning check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = Project;
mask = p.getMask('RATE');
subs = intersect(find(mask),Project.subjects_1500);
g = Group(subs);
[ratings sds]= g.getRatings(2:4);
%%
clf
c = 0;
for n = 2:4
    c = c+1;
    subpl(c) =  subplot(3,3,c);
    b = bar(-135:45:180,mean(ratings(:,:,n)));
    hold on;
    e = errorbar(-135:45:180,mean(ratings(:,:,n)),std(ratings(:,:,n))./sqrt(size(ratings,1)),'k.');
    set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'})
    SetFearGenBarColors(b)
    set(e,'LineWidth',2,'Color','k')
    ylim([0 10])
    xlim([-180 225])
    axis square
    box off
end
subplot(3,3,1);ylabel('Rating of p(shock)')
%add Groupfit line

x = -150:0.1:195;
subplot(3,3,1);
plot(x,mean(mean(ratings(:,:,2))),'k-','LineWidth',1.5)
subplot(3,3,2);
plot(x,VonMises(deg2rad(x),5.4473,1.3072,deg2rad(-8.0356),1.7061),'k-','LineWidth',1.5)
subplot(3,3,3);
plot(x,VonMises(deg2rad(x),4.6438,2.2963,deg2rad(-1.4574),1.9557),'k-','LineWidth',1.5)

%%scr 
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_1500BCT.mat','subjects','scr_bars')
ylims = [0 0.7];%600 0.3 1.1  %1500: 0.7
Yticks = [0 0.3 0.6];%600 [0.3 0.6 0.9], %1500: [0 0.3 0.6]

% plot
for n = [1 3]
subpl(3+n) = subplot(3,3,3+n);
b = bar(-135:45:180,mean(scr_bars(:,1:8,n)));
hold on;
ylim(ylims)
xlim([-180 225])
e = errorbar(-135:45:180,mean(scr_bars(:,1:8,n)),std(scr_bars(:,1:8,n))./sqrt(size(scr_bars,1)),'k.');
set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',Yticks)
SetFearGenBarColors(b)
set(e,'LineWidth',2,'Color','k')
axis square
box off
end

subplot(3,3,4);ylabel('SCR [muS]')
%cond special treatment
subpl(5) = subplot(3,3,5);

b(1) = bar(4,mean(scr_bars(:,4,2)));
hold on;
e(1) = errorbar(4,mean(scr_bars(:,4,2)),std(scr_bars(:,4,2))./sqrt(size(scr_bars,1)),'k.');
ylim(ylims)
set(gca,'XTick',1:8,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',Yticks);
b(2) =  bar(8,mean(scr_bars(:,8,2)));
e(2) = errorbar(8,mean(scr_bars(:,8,2)),std(scr_bars(:,8,2))./sqrt(size(scr_bars,1)));
set(b(1),'FaceColor',[1 0 0],'EdgeColor','w');
set(b(2),'FaceColor','c','EdgeColor','w');
set(e,'LineWidth',2,'Color','k');
axis square
box off
xlim([0 9])
try
    for n = 1:size(scr_bars,3)
        subplot(3,3,3+n)
        line([-180 225],repmat(mean(scr_bars(:,9,n)),[2 1]),'Color','k','LineWidth',1.3,'LineStyle',':')
    end
end

% add von mises groupfit curve to SCR bars @testphase
subplot(3,3,6);
plot(x,VonMises(deg2rad(x),0.1678,2.2181,deg2rad(-5.0893),0.2204),'k-','LineWidth',1.5)
subplot(3,3,4)
plot(x,mean(mean(scr_bars(:,1:8,1))),'k-','LineWidth',1.5)
% single subject example
% sub = 7;
% g2 = Group(sub); g2.ModelRatings(2,8);g2.ModelRatings(3,8);g2.ModelRatings(4,8);
% [rate sds] = g2.getRatings(2:4);

for sub = 1:length(g.ids)
    for n = 7:9
        subplot(3,3,n)
        b=bar(-135:45:180,rate(:,:,n-5));
        SetFearGenBarColors(b);
        hold on;
        errorbar(-135:45:180,rate(:,:,n-5),sds(:,:,n-5)./sqrt(2),'k.','LineWidth',2)
        set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',[0,5,10])
        params = g2.tunings.rate{n-5}.singlesubject{1}.Est(1:4);
        plot(x,VonMises(deg2rad(x),params(1),params(2),deg2rad(params(3)),params(4)),'k-','LineWidth',1.5)
        ylim([0 11])
        axis square
        box off
    end
end

%% SEARCH
% single subject example
clf
for sub = 1:length(g.ids)
    close all
    for n = 7:9
        subplot(3,3,n)
        b=bar(-135:45:180,ratings(sub,:,n-5));
        SetFearGenBarColors(b);
        hold on;
        errorbar(-135:45:180,ratings(sub,:,n-5),sds(sub,:,n-5)./sqrt(2),'k.','LineWidth',2)
        set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',[0,5,10])
        ylim([0 10])
        title(num2str(sub))
    end
    pause
end
%% add generalization in fixations
p =Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
fix = Fixmat(subjects,2:4);
counts0 = fix.histogram;close;
%cond: 11=null, 10=oddball, 9=ucs
for n = 1:3
    counts(:,:,n) = nanzscore(counts0(:,1:8,n)')';
end
c = 0;
for n = 1:3
    c = c+1;
    subpl(c) =  subplot(3,3,6+c);
    b = bar(1:8,mean(counts(:,1:8,n)));
    hold on;
    e = errorbar(1:8,mean(counts(:,1:8,n)),std(counts(:,1:8,n))./sqrt(size(counts,1)),'k.');
    set(gca,'XTick',1:8,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'})
    SetFearGenBarColors(b)
    set(e,'LineWidth',2,'Color','k')
    xlim([0 9])
    axis square
    box off
end
subplot(3,3,7);ylabel('fixation count')
%%
% ttest
[h,pval,ci,stats] = ttest(scr_bars(:,4,2),scr_bars(:,8,2));%conditioning check

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
clf
subplot(1,2,1)
b = bar([1 2 4 5 8 9 11 12 15 16 18 19],[mean(mat15(:,1)) mean(mat15(:,2)) mean(mat15(:,3)) mean(mat15(:,4))...
    mean(mat6(:,1)) mean(mat6(:,2)) mean(mat6(:,3)) mean(mat6(:,4)) mean(mat(:,1)) mean(mat(:,2)) mean(mat(:,3)) mean(mat(:,4))]);
hold on;
e = errorbar([1 2 4 5 8 9 11 12 15 16 18 19],[mean(mat15(:,1)) mean(mat15(:,2)) mean(mat15(:,3)) mean(mat15(:,4))...
    mean(mat6(:,1)) mean(mat6(:,2)) mean(mat6(:,3)) mean(mat6(:,4)) mean(mat(:,1)) mean(mat(:,2)) mean(mat(:,3)) mean(mat(:,4)) ],...
    [std(mat15(:,1))./sqrt(25) std(mat15(:,2))./sqrt(25) std(mat15(:,3))./sqrt(25) std(mat15(:,4))./sqrt(25) std(mat6(:,1))./sqrt(26) std(mat6(:,2))./sqrt(26) std(mat6(:,3))./sqrt(26) std(mat6(:,4))./sqrt(26)...
    std(mat(:,1))./sqrt(51) std(mat(:,2))./sqrt(51) std(mat(:,3))./sqrt(51) std(mat(:,4))./sqrt(51) ],'k.');
set(e,'LineWidth',1.5)
set(b,'FaceColor','white')
set(gca,'xtick',[1.5 4.5 8.5 11.5 15.5 18.5],'xticklabel',{'CS+' 'CS-' 'CS+' 'CS-' 'CS+' 'CS-'})
ylim([35 75])
xlim([-1 21])
box off


%% plot differences for CSP and CSN
subplot(1,2,2)
b = bar([1 2 4 5 7 8],[mean(mat15(:,1)-mat15(:,2)) mean(mat15(:,3)-mat15(:,4)) mean(mat6(:,1)-mat6(:,2)) mean(mat6(:,3)-mat6(:,4)) mean(mat(:,1)-mat(:,2)) mean(mat(:,3)-mat(:,4))]);
hold on;
e = errorbar([1 2 4 5 7 8],[mean(mat15(:,1)-mat15(:,2)) mean(mat15(:,3)-mat15(:,4)) mean(mat6(:,1)-mat6(:,2)) mean(mat6(:,3)-mat6(:,4)) mean(mat(:,1)-mat(:,2)) mean(mat(:,3)-mat(:,4))],...
   [ std(mat15(:,1)-mat15(:,2))./sqrt(length(mat15)) std(mat15(:,3)-mat15(:,4))./sqrt(length(mat15)) ...
   std(mat6(:,1)-mat6(:,2))./sqrt(length(mat6)) std(mat6(:,3)-mat6(:,4))./sqrt(length(mat6)) std(mat(:,1)-mat(:,2))./sqrt(length(mat)) std(mat(:,3)-mat(:,4))./sqrt(length(mat))],'k.');
set(e,'LineWidth',2)
set(gca,'xtick',[1.5 4.5 7.5],'xticklabel',{'1500ms' '600ms' 'pooled' })
set(e,'LineWidth',2)
set(b,'FaceColor','white')
box off
ylim([-5 35])
%% plot mean(CSP - CSN)
subplot(1,2,2)
b = bar([mean(mat15(:,11)) mean(mat6(:,11)) mean(mat(:,11))]);
hold on;
e = errorbar([mean(mat15(:,11)) mean(mat6(:,11)) mean(mat(:,11))],...
  [std(mat15(:,11))./sqrt(length(mat15)) std(mat6(:,11))./sqrt(length(mat6)) std(mat(:,11))./sqrt(length(mat6))],'k.');
set(e,'LineWidth',1.5)
set(b,'FaceColor','white','BarWidth',0.3)
box off
ylim([-3 20])
set(gca,'xticklabel',{'1500ms' '600ms' 'pooled' },'YTick',[0:5:20])


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

%% is improvement connected with SCR parameters?
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_1500BCT.mat','subjects','scr_bars')
subscr = subjects;clear subjects;
p = Project;
mask = p.getMask('PMF');
subpmf = intersect(find(sum(mask,2)==4),Project.subjects_1500);
valid = intersect(subpmf,subscr);
g = Group(valid);
mat = g.parameterMat;
[r,pval]=corr(mat(:,11),scr_bars(ismember(subscr,valid),4,2)-scr_bars(ismember(subscr,valid),8,2))%ph3 = test
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
%% check only conditioned people
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_1500BCT.mat')
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\pmf_data_N51.mat','mat','tags','mat15','mat6')
p    = Project;
mask = p.getMask('PMF');
Subjects = intersect(find(sum(mask,2)==4),[Project.subjects_600]);
subjects_pmf = Subjects;
subjects_scr = intersect(find(sum(p.getMask('SCR'),2)==3),Project.subjects_600);
valid_pmf = ismember(subjects_pmf,subjects_scr);

valid_scr = ismember(subjects_scr,subjects_pmf);
scr15pmfscr = scr_bars(valid_scr,:,:);
%take ratings also now
subjects_rate = intersect(find(p.getMask('RATE')),Project.subjects_600);
valid_subs = subjects_rate(ismember(subjects_rate,intersect(subjects_pmf,subjects_scr)));

%1500
g = Group(valid_subs);g.getSI(8);dummy = g.getSCRmeans('cond$');SCR_crit15 = dummy(:,1:2);
[mat15 tags] = g.parameterMat;
mat15 = [mat15 SCR_crit15];
tags{end+1} = 'SCR_csp';
tags{end+1} = 'SCR_csn';
%improvement CSP>CSN for SCR-cond people and SCR_notcond people
bar(1:2,[mean(mat15(SCR_crit15(:,1)>SCR_crit15(:,2),9)-mat15(SCR_crit15(:,1)>SCR_crit15(:,2),10)),mean(mat15(SCR_crit15(:,1)<=SCR_crit15(:,2),9)-mat15(SCR_crit15(:,1)<=SCR_crit15(:,2),10))])
hist(mat15(SCR_crit15(:,1)>SCR_crit15(:,2),9)-mat15(SCR_crit15(:,1)>SCR_crit15(:,2),10),20)
hold on
hist(mat15(SCR_crit15(:,1)<=SCR_crit15(:,2),9)-mat15(SCR_crit15(:,1)<=SCR_crit15(:,2),10),20,'Color','r')
[h,p,ci,stats]= ttest2(mat15(SCR_crit15(:,1)>SCR_crit15(:,2),9)-mat15(SCR_crit15(:,1)>SCR_crit15(:,2),10),mat15(SCR_crit15(:,1)<=SCR_crit15(:,2),9)-mat15(SCR_crit15(:,1)<=SCR_crit15(:,2),10))

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
%%
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\ratings_mises.mat')

figure
bar([1 2 3 6 7 8 11 12 13],[mean(mat15(:,12)) mean(mat15(:,13)) mean(mat15(:,14)) mean(mat6(:,12)) mean(mat6(:,13)) mean(mat6(:,14)) mean(mat(:,12)) mean(mat(:,13)) mean(mat(:,14))],'facecolor','white')
hold on;
e = errorbar([1 2 3 6 7 8 11 12 13],[mean(mat15(:,12)) mean(mat15(:,13)) mean(mat15(:,14)) mean(mat6(:,12)) mean(mat6(:,13)) mean(mat6(:,14)) mean(mat(:,12)) mean(mat(:,13)) mean(mat(:,14))],...
    [std(mat15(:,12))./sqrt(length(mat15)) std(mat15(:,13))./sqrt(length(mat15)) std(mat15(:,14))./sqrt(length(mat15))...
    mean(mat6(:,12))./sqrt(length(mat6)) mean(mat6(:,13))./sqrt(length(mat6)) mean(mat6(:,14))./sqrt(length(mat6))...
    std(mat(:,12))./sqrt(length(mat)) std(mat(:,13))./sqrt(length(mat)) std(mat(:,14))./sqrt(length(mat))],'k.')
set(gca,'xticklabel',{'Cond' 'Test' 'SI' 'Cond' 'Test' 'SI' 'Cond' 'Test' 'SI'})
set(e,'LineWidth',2)
[h,p,ci,stats] = ttest(mat(:,12),mat(:,13))
[h,p,ci,stats] = ttest(mat15(:,12),mat15(:,13))
[h,p,ci,stats] = ttest(mat6(:,12),mat6(:,13))

% plot bars and Gaussians in Feargen Colors
%
% g.ModelRatings(3,8)
% g.ModelRatings(4,8)
g.tunings.rate{3} = Tuning(g.Ratings(3));
g.tunings.rate{4} = Tuning(g.Ratings(4));
g.tunings.rate{3}.GroupFit(8)
g.tunings.rate{4}.GroupFit(8)
g.PlotRatingResults

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
VerticalXlabel(tags,'interpreter','none')



