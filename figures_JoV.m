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
%three exemplary subjects / discr
v = [];
c = 0;
for ph = [1 4]
    for sub = [6 15 31]
        c = c+1;
        v{c} = {'subject' sub 'phase' ph};
    end
end
fix.getmaps(v{:});
fix.plot
%%
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
%% histogram of alpha and kappa for both durations
clf
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate2pmf_N54.mat')
ind = subs15;
sp = [2 4];
crit = 22.5; % for later noshifter vs shifter improvement

subplot(sp(1),sp(2),1)%FWHM
[y,xi] = ksdensity(mat(ind,18),'bandwidth',10);
plot(xi,y,'b','LineWidth',2);hold on;
[y,xi] = ksdensity(mat(ind,19),'bandwidth',10);
plot(xi,y,'k','LineWidth',2)
axis square;box off
hold off
% xlim([min(abs(mat(ind,16))) max(abs(mat(:,16)))])\
xlim([0 180])
% set(gca,'YTick',yticks,'XTick',xticks,'XTickLabel',{'wide' 'sharp'})
ylabel('density')
title('FWHM')
text(10,7.8e-3,'A','FontSize',20);
subplot(sp(1),sp(2),2)%MU
[y,xi] = ksdensity(abs(mat(ind,15)),'bandwidth',7);
plot(xi,y,'b','LineWidth',2);hold on;
[y,xi] = ksdensity(abs(mat(ind,16)),'bandwidth',7);
plot(xi,y,'k','LineWidth',2)
legend('Conditioning','Generalization')
legend boxoff
axis square;box off
xlim([min(abs(mat(ind,16))) max(abs(mat(:,16)))])
% set(gca,'YTick',yticks,'XTick',xticks,'XTickLabel',{'wide' 'sharp'})
ylabel('density')
title('MU')
line([crit crit],ylim,'Color','k')
hold off
%
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat')
% load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\max4.mat')
ind = logical(abs(mat6(:,16))<= crit);
subplot(sp(1),sp(2),5)%alpha results
text(.5,78,'B','FontSize',20);
ylabel('threshold alpha (deg)')
b = bar([1 2 4 5],mean(mat6(:,1:4)));
SetFearGenBarColors(b,[1 0 0;.8 0 0;0 1 1;0 .8 .8]')
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
hold off;
%% graph alpha x discrimination (binning)
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat');
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
    mean(mean(mat(i,[1 3]),2))   
    binalpha(4,n) = mean(mean(mat(i,[1 3]),2));
    outcome_m(n) = mean(mean(mat(i,[12]),2));
    outcome_s(n) = std(mean(mat(i,[12]),2))./sqrt(sum(i)); 
end
bar(1:3,outcome_m,'FaceColor','r');hold on;
bar(1:3,[2.59 .995 .147],'FaceColor',[.3 .3 .3])
errorbar(1:3,outcome_m,outcome_s,'k.','LineWidth',linewidth)
set(gca,'XTickLabel',{'good','moderate','weak'})
xlabel('discrimination performance')
ylabel('fear specificity (kappa)')
axis square
box off
xlim([0 4])
hold off
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 2 - beneficial locations
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


%% SCR for dendrogram
clf
scr1 = scr_data2(:,1:8,1:8);% -repmat(scr_data2(:,9,1:8),[1 8 1]);
scr2 = scr_data2(:,1:8,9:end);% -repmat(scr_data2(:,9,9:end),[1 8 1]);

subplot(2,2,1);
SetFearGenColors;
plot(mean(scr1,3),'LineWidth',2)
xlim([150 700])
ylim([0 1])
axis square;box off;
ylabel('SCR [muS]')
set(gca,'XTick',[250 550],'XTickLabel',{'2.5','5.5'},'YTick',[0 .5 1],'FontSize',12)

subplot(2,2,2)
SetFearGenColors;
plot(mean(scr2,3),'LineWidth',2)
xlim([150 700])
ylim([0 1])
axis square;box off;
set(gca,'XTick',[250 550],'XTickLabel',{'2.5','5.5'},'YTick',[0 .5 1],'FontSize',12)


scr1 = scr_bars(1:length(sub1),1:8,3);%-repmat(scr_bars(1:length(sub1),9,3),[1 8]);
scr2 = scr_bars(9:end,1:8,3);%-repmat(scr_bars(9:end,9,3),[1 8]);

subplot(2,2,3);
b = bar(1:8,mean(scr1));SetFearGenBarColors(b);
hold on;
e = errorbar(1:8,mean(scr1),std(scr1)./sqrt(length(scr1)),'k.','LineWidth',1.5);
axis square;box off;
line(xlim,repmat(mean(scr_bars(1:8,9,3)),[2 1]),'Color','k','LineWidth',1.3,'LineStyle',':');
ylim([0 .6])
ylabel('SCR tuning [muS]')
set(gca,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'},'YTick',[0 .3 .6],'FontSize',12)


subplot(2,2,4);
b = bar(1:8,mean(scr2));SetFearGenBarColors(b);
hold on;
e = errorbar(1:8,mean(scr2),std(scr2)./sqrt(length(scr2)),'k.','LineWidth',1.5);
axis square;box off;
line(xlim,repmat(mean(scr_bars(9:end,9,3)),[2 1]),'Color','k','LineWidth',1.3,'LineStyle',':');
ylim([0 .6])
set(gca,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'},'YTick',[0 .3 .6],'FontSize',12)


%% SCR for dendrogram (not nulltrial corrected)
scr1 = scr_bars(1:length(sub1),1:8,3);
scr2 = scr_bars(9:end,1:8,3);

clf
subplot(3,2,1);
b = bar(1:8,mean(scr1));SetFearGenBarColors(b);
hold on;
e = errorbar(1:8,mean(scr1),std(scr1)./sqrt(length(scr1)),'k.','LineWidth',1.5);
axis square;box off;
set(gca,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'})

subplot(3,2,2);
b = bar(1:8,mean(scr2));SetFearGenBarColors(b);
hold on;
e = errorbar(1:8,mean(scr2),std(scr2)./sqrt(length(scr2)),'k.','LineWidth',1.5);
axis square;box off;
set(gca,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'})
EqualizeSubPlotYlim(gcf)


scr1 = scr_data2(:,1:8,1:8);
scr2 = scr_data2(:,1:8,9:end);

subplot(3,2,3);
SetFearGenColors;
plot(mean(scr1,3),'LineWidth',2)
xlim([0 800])
axis square;box off;

subplot(3,2,4)
SetFearGenColors;
plot(mean(scr2,3),'LineWidth',2)
xlim([0 800])
axis square;box off;

subplot(3,2,5)
b = bar(1:8,mean(mean(scr1,1),3));
hold on;
SetFearGenBarColors(b);
e = errorbar(1:8,mean(mean(scr1,1),3),std(mean(scr1,1))./sqrt(16),'k.','LineWidth',1.5);
axis square;box off;
subplot(3,2,6)
b = bar(1:8,mean(mean(scr2,1),3));
hold on;
SetFearGenBarColors(b);
e = errorbar(1:8,mean(mean(scr2,1),3),squeeze(std(mean(scr2,1)))./sqrt(16),'k.','LineWidth',1.5);


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
