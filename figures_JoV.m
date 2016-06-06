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
%histogram of alpha and kappa
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat');
subplot(2,2,1)

subplot(2,2,2)
hist(mean(mat15(:,[1 3]),2),9);
axis square;box off
xlabel('alpha before')
subplot(2,2,3)
hist(mat15(:,12),9);
axis square;box off
set(gca,'YTick',0:5:10)
xlabel('kappa cond')
subplot(2,2,4)
hist(mat15(:,13),9);
axis square;box off
set(gca,'YTick',0:5:10)
xlabel('kappa test')
EqualizeSubPlotYlim(gcf)

%% graph alpha x discrimination (binning)
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat');
linewidth = 3;
x_pos = [1.05:1:3.05;.95:1:2.95;1.05:1:3.05;.95:1:2.95];
% bin perception, plot kappa
clf
subplot(1,2,1)
% 1500ms
col = [1 3];
mat = mat15;
p = prctile(mean(mat(:,col),2),linspace(0,100,4));

outcome_m = [];outcome_s=[];
for n = 1:length(p)-1
    i = (mean(mat(:,col),2) > p(n)) & (mean(mat(:,col),2) <= p(n+1));
    sum(i);
    binalpha(1,n) = mean(mean(mat(i,[1 3]),2));
    outcome_m(n) = mean(mean(mat(i,[12]),2));
    outcome_s(n) = std(mean(mat(i,[12]),2))./sqrt(sum(i));
end

errorbar(x_pos(1,:),outcome_m,outcome_s,':','LineWidth',linewidth)
hold on;
% 600ms
mat = mat6;
p = prctile(mean(mat(:,col),2),linspace(0,100,4));

outcome_m = [];outcome_s=[];
for n = 1:length(p)-1    
    i = (mean(mat(:,col),2) > p(n)) & (mean(mat(:,col),2) <= p(n+1));
    sum(i);
    binalpha(2,n) = mean(mean(mat(i,[1 3]),2));
    outcome_m(n) = mean(mean(mat(i,[12]),2));
    outcome_s(n) = std(mean(mat(i,[12]),2))./sqrt(sum(i)); 
end
errorbar(x_pos(2,:),outcome_m,outcome_s,'LineWidth',linewidth)
legend('1500 ms','600 ms')
legend boxoff
%
subplot(1,2,2)
title('Test Phase')
% 1500ms
col = [1 3];
mat = mat15;
p = prctile(mean(mat(:,col),2),linspace(0,100,4));

outcome_m = [];outcome_s=[];
for n = 1:length(p)-1
    i = (mean(mat(:,col),2) > p(n)) & (mean(mat(:,col),2) <= p(n+1));
    sum(i);
    binalpha(3,n) = mean(mean(mat(i,[1 3]),2));
    outcome_m(n) = mean(mean(mat(i,[13]),2));
    outcome_s(n) = std(mean(mat(i,[13]),2))./sqrt(sum(i));
end
errorbar(x_pos(3,:),outcome_m,outcome_s,':','LineWidth',linewidth)
hold on;
% 600ms
mat = mat6;
p = prctile(mean(mat(:,col),2),linspace(0,100,4));

outcome_m = [];outcome_s=[];
for n = 1:length(p)-1    
    i = (mean(mat(:,col),2) > p(n)) & (mean(mat(:,col),2) <= p(n+1));
    sum(i);
    binalpha(4,n) = mean(mean(mat(i,[1 3]),2));
    outcome_m(n) = mean(mean(mat(i,[13]),2));
    outcome_s(n) = std(mean(mat(i,[13]),2))./sqrt(sum(i)); 
end
errorbar(x_pos(4,:),outcome_m,outcome_s,'LineWidth',linewidth)
legend('1500 ms','600 ms')
legend boxoff
for n=1:2
subplot(1,2,n)
set(gca,'XTick',1:3,'XTicklabel',{'good','medium','bad'},'YTick',2:2:10);
ylim([0 10])
ylabel('FearGen specificity (kappa)')
xlabel('Discrimination performance')
box off
axis square
end
subplot(1,2,1);title('Conditioning')
subplot(1,2,2);title('Test phase')



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
