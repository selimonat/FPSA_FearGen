clear all;
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
fix             = Fixmat(subjects,[1 2 3 4 5]);

fix.getsubmaps;
fix.maps = mean(fix.maps,3);
fix.plot
%discr vs feargen
v = [];
c = 0;
for ph = [1 4]
    c = c+1;
    v{c} = {'phase' ph};
end
fix.getmaps(v{:});
fix.plot('linear',[2 1])

%%%%%%%%%%%%%
% behavioral results 
%%%%%%%%%%%%%
%% big figure for behavioral results
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\ratings_600.mat','ratings');
sp = [2 5];
clf
c = 6;
for n = 2:4
    c = c+1;
    subpl(c) =  subplot(sp(1),sp(2),c);
    b = bar(-135:45:180,mean(ratings(:,:,n)));
    hold on;
    e = errorbar(-135:45:180,mean(ratings(:,:,n)),std(ratings(:,:,n))./sqrt(size(ratings,1)),'k.');
    set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',[0 5 10])
    SetFearGenBarColors(b)
    set(e,'LineWidth',2,'Color','k')
    ylim([0 10])
    xlim([-180 225])
    axis square
    box off
end
subplot(sp(1),sp(2),7);ylabel('Rating of p(shock)')
%add Groupfit line
params = [5.1801 0.7288 deg2rad(0.8714) 1.8312; 4.4396 1.1171 deg2rad(-0.0598) 1.9506]; %600ms cond;test
% params = [5.4473,1.3072,deg2rad(-8.0356),1.7061;4.6438,2.2963,deg2rad(-1.4574),1.9557] %1500 cond;test
x = -150:0.1:195;
subplot(sp(1),sp(2),7);
plot(x,mean(mean(ratings(:,:,2))),'k-','LineWidth',3.5)
subplot(sp(1),sp(2),8);
plot(x,VonMises(deg2rad(x),params(1,1),params(1,2),params(1,3),params(1,4)),'k-','LineWidth',1.5)
line([0 180],[9 9],'Color','k','LineWidth',1.3)
text(45,9,'***','FontSize',20)
subplot(sp(1),sp(2),9);
plot(x,VonMises(deg2rad(x),params(2,1),params(2,2),params(2,3),params(2,4)),'k-','LineWidth',1.5)
line([0 180],[9 9],'Color','k','LineWidth',1.3)
text(45,9,'***','FontSize',20)

%%scr 
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_600BCT.mat','scr_bars','subjects')
ylims = [0.3 .9];%600 0.3 1.1  %1500: 0.7
Yticks = [0.3 0.6 0.9];%600 [0.3 0.6 0.9], %1500: [0 0.3 0.6]

% plot
for n = [1 3]
subplot(sp(1),sp(2),n+1);
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
subplot(sp(1),sp(2),2);ylabel('SCR [muS]')
title('Free Viewing')
subplot(sp(1),sp(2),4)
title('Generalization')
%cond special treatment
subplot(sp(1),sp(2),3);
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
line([4 8],[.85 0.85],'Color','k','LineWidth',1.3)
text(5.7,.85,'*','FontSize',20)
title('Conditioning')

try
    for n = 1:size(scr_bars,3)
        subplot(sp(1),sp(2),1+n)
        line([-180 225],repmat(mean(scr_bars(:,9,n)),[2 1]),'Color','k','LineWidth',1.3,'LineStyle',':')
    end
end
% 
% % add von mises groupfit curve to SCR bars @testphase
subplot(sp(1),sp(2),2);
plot(x,mean(mean(scr_bars(:,1:8,1))),'k-','LineWidth',3.5)
subplot(sp(1),sp(2),4);
plot(x,mean(mean(scr_bars(:,1:8,3))),'k-','LineWidth',3.5)

% pmf at 1 and 5
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\params_pmf_600.mat','params')
pmf = mean(params,3);
pmf(:,2) = repmat(mean(params(:,2)),[4 1]);
x = linspace(0,180,1000);
subplot(sp(1),sp(2),1);
plot(x,PAL_Weibull(pmf(2,:),x),'c','LineWidth',3)%csn before
hold on;
plot(x,PAL_Weibull(pmf(1,:),x),'r','LineWidth',3)%csp before
set(gca,'YTick',0:.25:1);
axis square
box off
ylim([0 1]);
xlim([0 120]);
subplot(sp(1),sp(2),5);
plot(x,PAL_Weibull(pmf(4,:),x),'c','LineWidth',3)%csn after
hold on;
plot(x,PAL_Weibull(pmf(3,:),x),'r','LineWidth',3)%csp after
set(gca,'YTick',0:.25:1);
axis square
box off
ylim([0 1]);
xlim([0 120]);
legend('CS-','CS+','location','southeast')
legend boxoff
hold off
%% histogram of alpha and kappa
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate2pmf_N54.mat')
ind = subs6;
sp = [2 4];
crit = 22.5; % for later noshifter vs shifter improvement

clf
subplot(sp(1),sp(2),1)%FWHM
hold off
b = linspace(20,180,6);
histf(mat(ind,18),b,'FaceColor','b','EdgeColor','none','FaceAlpha',.6);
hold on;
histf(mat(ind,19),b,'FaceColor','k','EdgeColor','none','FaceAlpha',.6)
axis square;box off
hold off
ylabel('counts')
title('FWHM')
text(10,7.8e-3,'A','FontSize',20);
ylim([0 5])
subplot(sp(1),sp(2),2)%MU
hold off
b = linspace(0,40,6);
histf(abs(mat(ind,15)),b,'FaceColor','b','EdgeColor','none','FaceAlpha',.6);
hold on;
histf(abs(mat(ind,16)),b,'FaceColor','k','EdgeColor','none','FaceAlpha',.6)
axis square;box off
hold off
ylabel('counts')
title('abs(Mu)')
hold off
ylim([0 5])
legend('Conditioning','Generalization','orientation','vertical','location','best')
legend boxoff
xlim([0 40])
%
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat')
% load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\max4.mat')
ind = logical(abs(mat6(:,16))<= crit);
subplot(sp(1),sp(2),5)%alpha results
ylabel('threshold alpha (deg)')
b = bar([1 2 4 5],mean(mat6(:,1:4)));
SetFearGenBarColors(b,[1 0 0;.6 0 0;0 1 1;0 .6 .6]')
line([1 2],[72 72],'Color','k','LineWidth',2)
text(1,75,'***','FontSize',16)
hold on;
e = errorbar([1 2 4 5],mean(mat6(:,1:4)),std(mat6(:,1:4))./sqrt(length(mat6)),'k.');
set(e,'LineWidth',2)
set(gca,'xtick',[1.5 4.5],'xticklabel',{'CS+' 'CS-'})
set(e,'LineWidth',2)
axis square;box off;
xlim([0 6])
subplot(sp(1),sp(2),6)%csp improvement shifters vs no shifters
b = bar(1:2,[mean(mat6(ind,11)) mean(mat6(~ind,11))]);
set(b,'BarWidth',.5,'FaceColor','k')
hold on;box off
e = errorbar(1:2,[mean(mat6(ind,11)) mean(mat6(~ind,11))],[std(mat6(ind,11))./sqrt(sum(ind)) std(mat6(~ind,11))./sqrt(sum(~ind))],'k.','LineWidth',1.5);
e = errorbar(1,mean(mat6(ind,11)),std(mat6(ind,11))./sqrt(sum(ind)),0,'w.','LineWidth',1.5);
axis square
xlim([0 3])
set(gca,'XTick',[1 2],'XTicklabel',{'CS+ centered' 'shifted'})
ylim([-5 20])
ylabel('corrected improvement (deg)')
text(.5,78,'B','FontSize',20);
hold off;
%% correlation
subplot(sp(1),sp(2),[3 4 7 8])
text(20,180,'C','FontSize',20);
hold on;
plot(mat6(:,17),mat6(:,18),'bo','MarkerFaceColor','b','MarkerSize',10)
l = lsline;set(l,'LineWidth',2)
box off
xlabel('initial alpha')
ylabel('FWHM Cond')
% regression
X = mat6(:,17); %initial alpha
Y = mat6(:,18); %FWHM cond
mdl = LinearModel.fit(X,Y);
plot(mdl)
legend({'Data' 'Fit' 'Confidence Interval' '' 'Identity Line'},'location','best')
legend boxoff
set(gca,'XTick',20:20:110,'YTick',20:20:180)
ylim([20 180])
DrawIdentityLine(gca);
xlim([20 110])
title('Discrimination x Generalization', 'FontSize',14)
set(gca,'FontSize',12)
xlabel('Discrimination (alpha [deg])')
ylabel('FWHM Conditioning')
axis square
box off
%% new way (tobe worked on)
%%
x  = (rand(30,1));
y  = (rand(30,1));
figure(1);
hold off
plot(x,y,'o')
hold on;
%
i = [1:30];
for n = 1:1000
    ii     = randsample(i,30,1);
    B(:,n) = [x(ii) ones(30,1)]\y(ii);
end
%
X2 = [linspace(0,1,100)' ones(100,1)];
Y = X2*B;
Y = Y';
YY = prctile(Y,[.5 99.5])';
hold on;
plot(X2(:,1),YY(:,1))
plot(X2(:,1),YY(:,2),'r')
hold off;


%% graph alpha x discrimination (binning)
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat');
sp = [2 4];
linewidth = 3;
subplot(sp(1),sp(2),[3 4 7 8])
text(.5,5.5,'C','FontSize',20);
% bin perception, plot kappa
% 600ms
mat = mat6;
col = [1 3];
p = prctile(mean(mat(:,col),2),linspace(0,100,4));
outcome_m = [];outcome_s=[];
for n = 1:length(p)-1   
    if n ==1
        i = (mean(mat(:,col),2) >= p(n)) & (mean(mat(:,col),2) <= p(n+1));
    else
        i = (mean(mat(:,col),2) >= p(n)) & (mean(mat(:,col),2) <= p(n+1));
    end
    sum(i)
    ii(:,n) = i;
    mean(mean(mat(i,[1 3]),2))   
    binalpha(n) = mean(mean(mat(i,[1 3]),2));
    outcome_m(n) = mean(mean(mat(i,[18]),2));
    outcome_s(n) = std(mean(mat(i,[18]),2))./sqrt(sum(i)); 
end
bar(1:3,binalpha,'FaceColor',[.7 .7 .7],'EdgeColor','none')%bar(1:3,[2.59 .995 .147],'FaceColor',[.3 .3 .3])
bar(1:3,outcome_m,'FaceColor','w');hold on;
errorbar(1:3,outcome_m,outcome_s,'k.','LineWidth',linewidth)
set(gca,'XTickLabel',{'good','moderate','weak'})
xlabel('discrimination performance')
ylabel('fear specificity (kappa)')
axis square
box off
xlim([0 4])
hold off

%
for n=1:length(ii)
    group(n) = find(ii(n,:));
end
figure
boxplot(mat6(:,18))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 2 - Task Demands x Subject Idiosyncrasy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
p           = Project;
mask = p.getMask('ET_discr');
subjects = intersect([Project.subjects_1500],find(mask));
mask = p.getMask('PMF');
subjects = intersect(subjects,find(sum(mask,2)==4));
fix = Fixmat(subjects,1:5);
v = [];
c = 0;
for ph = 1:5
    c = c+1;
    v{c} = {'phase' ph};
end
fix.getmaps(v{:});
fix.plot
fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 5]);
fix.plot
v = [];
c = 0;
for sub = [6 15 31]
    for ph = 1:5
        c = c+1;
        v{c} = {'subject' sub 'phase' ph 'deltacsp' fix.realcond};
    end
end
fix.getmaps(v{:});
fix.plot

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 3 - beneficial locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% discrimination
clear all
p           = Project;
mask = p.getMask('ET_discr');
subjects = intersect([Project.subjects_1500],find(mask));
mask = p.getMask('PMF');
subjects = intersect(subjects,find(sum(mask,2)==4));
g = Group(subjects);[mat tags] = g.parameterMat;
clear g
fix = Fixmat(subjects, 5);
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmaps = fix.vectorize_maps;

%normal correlation
for n = 1:length(subjmaps)
    [r(n),pval(n)] = corr(subjmaps(n,:)',mean(mat(:,[2 4]),2));
end
imagesc(reshape(r,[50 50]),[-1 1])
for n = find(pval<0.05);[y x] = ind2sub([50 50],n);text(x,y,'x','FontSize',20);end
fix.maps = reshape(r,[50 50]);
fix.plot
% correct it by mean fixmap
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = Scale(mean(fix.maps,3));
rc       = r.*fix.maps(:)';
fix.maps = reshape(rc,[50 50]);
fix.plot
% get the numbers
[mini,i] = min(rc)
r(i)
pval(i)
[y,x] = ind2sub([50 50], i)
% -> minimal pixel for rc doesn't get significant in the original
% correlation..

%% fear specificity
clear all
p           = Project;
mask = p.getMask('ET_feargen');
subjects = intersect([Project.subjects_1500],find(mask));
mask = p.getMask('RATE');
subjects = intersect(subjects,find(mask));
g = Group(subjects);[mat tags] = g.parameterMat; mises = g.loadmises; mat = [mat mises];
clear g
fix = Fixmat(subjects, 1);
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmaps = fix.vectorize_maps;

%normal correlation
for n = 1:length(subjmaps)
    [r(n),pval(n)] = corr(subjmaps(n,:)',mat(:,13));
end
imagesc(reshape(r,[50 50]),[-1 1])
for n = find(pval<0.05);[y x] = ind2sub([50 50],n);text(x,y,'x','FontSize',20);end
fix.maps = reshape(r,[50 50]);
fix.plot
% correct it by mean fixmap
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = Scale(mean(fix.maps,3));
rc       = r.*fix.maps(:)';
fix.maps = reshape(rc,[50 50]);
fix.plot
% get the numbers
[maxi,i] = max(r)
r(i)
pval(i)
[y,x] = ind2sub([50 50], i)
% -> minimal pixel for rc doesn't get significant in the original
% correlation..

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig X - Hierarchical Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
fix             = Fixmat(subjects,4);

% get the data
g = Group(subjects);
[mat tags] = g.parameterMat;
mises = g.loadmises;
mat = [mat mises];
valid_pmf = ismember(subjects,find(sum(p.getMask('PMF'),2)==4));
valid_rate = ismember(subjects,find(p.getMask('RATE')));
mat(~valid_pmf,1:11) = NaN;
mat(~valid_rate,12:end) = NaN;
clear g
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
[branch_id,order] = fix.dendrogram(3,mat(:,13));
%% 
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\dendro_branch_id_subjects.mat')
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\mat_N27_for_dendro.mat')
g1 = Group(sub1);
g2 = Group(sub2);
mat1 = g1.loadmises;
mat2 = g2.loadmises;
%% DENDROGRAM BEHAVIORAL CLUSTER OUTCOMES
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\mat_N27_for_dendro.mat','mat','mat1','mat2','sub1','sub2','subjects','tags','branch_id')
fonts = 12;
clf
subplot(4,2,1:2);
errorbar([nanmean(mat1(:,2)) nanmean(mat2(:,2))],[nanstd(mat1(:,2))./sqrt(sum(~isnan(mat1(:,2)))) nanstd(mat2(:,2))./sqrt(sum(~isnan(mat2(:,2))))],'k','LineWidth',2)

% for n = 1:8;  fwhm1(n) = vM2FWHM(mat1(n,2));end
% for n = 1:17; fwhm2(n) = vM2FWHM(mat2(n,2));end
% subplot(4,2,1:2);
% errorbar([nanmean(fwhm1) nanmean(fwhm2)],[nanstd(fwhm1)./sqrt(sum(~isnan(mat1(:,2)))) nanstd(fwhm2)./sqrt(sum(~isnan(mat2(:,2))))],'k','LineWidth',2)
set(gca,'XTick',1:2,'XTickLabel',{'Cluster 1' 'Cluster 2'},'YTick',0:5:15,'FontSize',fonts)
ylabel('Fear Specificity (Kappa)','FontSize',fonts)
box off
line([2.1 2.1],[nanmean(mat1(:,2)) nanmean(mat2(:,2))],'LineWidth',2)
text(2.15,mean([nanmean(mat1(:,2)) nanmean(mat2(:,2))]),'*','FontSize',30);
%% SCR for dendrogram
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_cut_fordendro.mat','scr_data_cut','scr_data','subjects','mask')
scr1 = scr_data(:,1:8,ismember(subjects,sub1));% -repmat(scr_data_cut(:,9,ismember(subjects,sub1)),[1 8 1]);
scr2 = scr_data(:,1:8,ismember(subjects,sub2));% -repmat(scr_data_cut(:,9,ismember(subjects,sub2)),[1 8 1]);

subplot(4,2,3);
SetFearGenColors;
plot(mean(scr1,3),'LineWidth',2)
hold on;
plot(squeeze(mean(scr_data(:,9,ismember(subjects,sub1)),3)),'k:','LineWidth',2)
ylim([0 1])
axis square;box off;
ylabel('SCR [muS]','FontSize',fonts)
set(gca,'XTick',[0 250 550],'XTickLabel',{'0','2.5','5.5'},'YTick',[0 .5 1],'FontSize',fonts)
line([250 250],[0 1],'LineStyle','--','Color',[.2 .2 .2],'LineWidth',1.5)
line([550 550],[0 1],'LineStyle','--','Color',[.2 .2 .2],'LineWidth',1.5)
xlim([150 650])

subplot(4,2,4)
SetFearGenColors;
plot(mean(scr2,3),'LineWidth',2)
hold on;
plot(squeeze(mean(scr_data(:,9,ismember(subjects,sub2)),3)),'k:','LineWidth',2)
ylim([0 1])
axis square;box off;
set(gca,'XTick',[0 250 550],'XTickLabel',{'0','2.5','5.5'},'YTick',[0 .5 1],'FontSize',fonts)
line([250 250],[0 1],'LineStyle','--','Color',[.2 .2 .2],'LineWidth',1.5)
line([550 550],[0 1],'LineStyle','--','Color',[.2 .2 .2],'LineWidth',1.5)
xlim([150 650])

scr1 = squeeze(mean(scr_data_cut(:,1:8,ismember(subjects,sub1))))';
scr2 = squeeze(mean(scr_data_cut(:,1:8,ismember(subjects,sub2))))';

subplot(4,2,5);
b = bar(-135:45:180,mean(scr1));SetFearGenBarColors(b);
hold on;
e = errorbar(-135:45:180,mean(scr1),std(scr1)./sqrt(length(scr1)),'k.','LineWidth',1.5);
x = linspace(-135,180,1000);
plot(x,Tuning.VonMises(x,0.1395,2.2234,1.4489,0.1292),'k','LineWidth',1.5)
axis square;box off;
line(xlim,repmat(mean(mean(scr_data_cut(:,9,ismember(subjects,sub1)))),[2 1]),'Color','k','LineWidth',1.5,'LineStyle',':');
ylim([0 .6])
xlim([-185 230])
ylabel('SCR tuning [muS]','FontSize',fonts)
set(gca,'XTick',[0 180],'XTickLabel',{'CS+' 'CS-'},'YTick',[0 .3 .6],'FontSize',fonts)


subplot(4,2,6);
b = bar(-135:45:180,mean(scr2));SetFearGenBarColors(b);
hold on;
e = errorbar(-135:45:180,mean(scr2),std(scr2)./sqrt(length(scr2)),'k.','LineWidth',1.5);
x = linspace(-135,180,1000);
plot(x,Tuning.VonMises(x,0.1798,2.5709,-6.3190,.2803),'k','LineWidth',1.5)
axis square;box off;
line(xlim,repmat(mean(mean(scr_data_cut(:,9,ismember(subjects,sub2)))),[2 1]),'Color','k','LineWidth',1.5,'LineStyle',':');
ylim([0 .6])
xlim([-185 230])
set(gca,'XTick',[0 180],'XTickLabel',{'CS+' 'CS-'},'YTick',[0 .3 .6],'FontSize',fonts)

%% alpha
p = Project;
alphmat1 = g1.parameterMat;
alphmat2 = g2.parameterMat;
invalid_pmf1 = ~ismember(g1.ids,find(sum(p.getMask('PMF'),2)==4));
alphmat1(invalid_pmf1,1:11) = nan;
invalid_pmf2 = ~ismember(g2.ids,find(sum(p.getMask('PMF'),2)==4));
alphmat2(invalid_pmf2,1:11) = nan;
alpha1 = mean(alphmat1(~invalid_pmf1,[1 3]),2);
alpha2 = mean(alphmat2(~invalid_pmf2,[1 3]),2);
subplot(4,2,7:8)
errorbar([mean(alpha1) mean(alpha2)],[std(alpha1)./sqrt(length(alpha1)) std(alpha2)./sqrt(length(alpha2))],'k','LineWidth',2)
set(gca,'XTick',1:2,'XTickLabel',{'Cluster 1' 'Cluster 2'},'FontSize',fonts)
ylabel('initial alpha [deg]','FontSize',fonts)
box off
[h,p,ci,stats] = ttest2(alpha1,alpha2,'tail','right')
line([2.1 2.1],[mean(alpha1) mean(alpha2)],'Color','k','LineWidth',2)
text(2.15,mean([mean(alpha1) mean(alpha2)]),'*','Color','k','FontSize',30);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig X - RSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/Users/onat/Documents/Code/Matlab/fancycarp/')
clear all;
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_600);
fix             = Fixmat(subjects,[1 2 3 4 5]);

%% plot and save single subject fixations maps RSA
tsub = length(unique(fix.subject));
out  = [];
cormat_t = nan(8,8,tsub);
cormat_b = nan(8,8,tsub);
pval_t = nan(8,8,tsub);
pval_b = nan(8,8,tsub);
subc            = 0;
% fix.kernel_fwhm = 25;
% fix.kernel_fwhm = 37*0.8;

for subject = unique(fix.subject);
    subc = subc + 1
    %creaete the query cell
    v = [];
    c = 0;
    for ph = [2 4]
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject' subject};
        end
    end
    % plot and save fixation maps
    fix.getmaps(v{:});
    maps = fix.maps;
%     fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:16),3),[1 1 16]);%correct for baseline
    fix.maps   = maps(:,:,9:end) - repmat(mean(maps(:,:,9:end),3),[1 1 8]);%take out the average    
    %         fix.plot
    %         SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/AllFixations/baselinecorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
    %         saveas(gcf,sprintf('C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/baselineANDMeancorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
    [cormat_t(:,:,subc) pval_t(:,:,subc)] = fix.corr;
    fix.maps   = maps(:,:,1:8) - repmat(mean(maps(:,:,1:8),3),[1 1 8]);%take out the average
    [cormat_b(:,:,subc) pval_b(:,:,subc)] = fix.corr;
end
% end
%% RSA plot
subplot(1,2,1);
imagesc(mean(cormat_b,3),[-.15 .15]);
axis square;colorbar;
set(gca,'fontsize',15)
axis off
subplot(1,2,2);
imagesc(mean(cormat_t,3),[-.15 .15])
axis square;colorbar
set(gca,'fontsize',15)
axis off

%% RSA by single fixations
nfix = 1:4;
tsub = length(unique(fix.subject));
tfix = length(fix);
cormat_t = nan(8,8,tsub,tfix);
cormat_b = nan(8,8,tsub,tfix);

for nf = nfix(:)';
    subc   = 0;
    for subject = unique(fix.subject);
        subc = subc + 1;
        %creaete the query cell
        v = [];
        c = 0;
        for ph = [2 4]
            for cond = -135:45:180
                c    = c+1;
                v{c} = {'phase', ph, 'deltacsp' cond 'fix' nf 'subject' subject};
            end
        end
        % plot and save fixation maps
        fix.getmaps(v{:});
        maps                 = fix.maps;
        fix.maps             = maps(:,:,1:8)   - repmat(mean(maps(:,:,1:8),3),[1 1 8]);%correct for baseline
        cormat_b(:,:,subc,nf)   = fix.corr;
        fix.maps   = maps(:,:,9:end) - repmat(mean(maps(:,:,9:end),3),[1 1 8]);%take out the average
        cormat_t(:,:,subc,nf)   = fix.corr;
    end
end
 %%
for nf = 1:4;
    subplot(2,4,nf)
    imagesc([nanmean(cormat_b(:,:,:,nf),3) ],[-.15 .15]);
    title(sprintf('Fixation No: %d',nf))
    if nf == 1
        ylabel('Baseline')
    end
    axis square
    axis off
    colorbar
    subplot(2,4,nf+4)
    imagesc([nanmean(cormat_t(:,:,:,nf),3)],[-.15 .15]);
    axis square
    colorbar
    axis off
    if nf == 1
        ylabel('Test')
    end
end
