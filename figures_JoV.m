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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% behavioral results
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% big figure for behavioral results
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\ratings_600.mat','ratings');
sp = [2 5];
for n = 2:4
    subpl(n) =  subplot(sp(1),sp(2),n);
    b = bar(-135:45:180,mean(ratings(:,:,n)));
    hold on;
    e = errorbar(-135:45:180,mean(ratings(:,:,n)),std(ratings(:,:,n))./sqrt(size(ratings,1)),'k.');
    set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',[0 5 10],'FontSize',12)
    SetFearGenBarColors(b)
    set(e,'LineWidth',2,'Color','k')
    ylim([0 10])
    xlim([-180 225])
    axis square
    box off
end
subplot(sp(1),sp(2),2);ylabel('Rating of p(shock)','FontSize',12)
%add Groupfit line
params = [5.1801 0.7288 deg2rad(0.8714) 1.8312; 4.4396 1.1171 deg2rad(-0.0598) 1.9506]; %600ms cond;test
% params = [5.4773,1.3072,deg2rad(-8.0356),1.7061;4.6438,2.2964,deg2rad(-1.4551),1.9557]; %1500 cond;test
x = -150:0.1:195;
subplot(sp(1),sp(2),2);
plot(x,mean(mean(ratings(:,:,2))),'k-','LineWidth',2)
subplot(sp(1),sp(2),3);
plot(x,VonMises(deg2rad(x),params(1,1),params(1,2),params(1,3),params(1,4)),'k-','LineWidth',2)
line([0 180],[9 9],'Color','k','LineWidth',2)
text(45,9,'***','FontSize',20)
subplot(sp(1),sp(2),4);
plot(x,VonMises(deg2rad(x),params(2,1),params(2,2),params(2,3),params(2,4)),'k-','LineWidth',2)
line([0 180],[9 9],'Color','k','LineWidth',2)
text(45,9,'***','FontSize',20)
subplot(sp(1),sp(2),2);
title('Free Viewing','FontSize',14)
subplot(sp(1),sp(2),3);
title('Conditioning','FontSize',14)
subplot(sp(1),sp(2),4)
title('Generalization','FontSize',14)

%% scr
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_600_zscores.mat')
ylims = [-1 1];%600 [-1 1]  %1500: [-1 1.8]
Yticks = [-1 0 1];%600 [-.8 0 .8], %1500: [-1 0 1]

% plot
for n = [7 8 9]
    ind = n-7;
    subplot(sp(1),sp(2),n);
    b = bar(-135:45:180,mean(DD([1:8]+8*ind,:),2));
    hold on;
    ylim(ylims)
    xlim([-180 225])
    e = errorbar(-135:45:180,mean(DD([1:8]+8*ind,:),2),std(DD([1:8]+8*ind,:),0,2)./sqrt(size(DD,2)),'k.');
    set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',Yticks,'FontSize',12)
    SetFearGenBarColors(b)
    set(e,'LineWidth',2,'Color','k')
    axis square
    box off
end
subplot(2,5,8);ylim([0 2.2]);set(gca,'YTick',0:1:2.0)
%cond phase CS+ CS-
subplot(sp(1),sp(2),8);
line([0 180],repmat(1.8,[1 2]),'Color','k','LineWidth',2)
text(90,1.8,'***','FontSize',20)
subplot(sp(1),sp(2),9);
line([0 180],repmat(.7,[1 2]),'Color','k','LineWidth',2)
text(90,.7,'*','FontSize',20)


subplot(sp(1),sp(2),7);ylabel('SCR (z-Score)','FontSize',12)
subplot(sp(1),sp(2),9);
hold on;
params = fit_results.params;
x = -150:0.1:195;
plot(x,VonMises(deg2rad(x),params(1),params(2),deg2rad(params(3)),params(4)),'k-','LineWidth',2)


%% pmf at 1 and 5
sp = [2 5];
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\pdiff_data.mat')
pdiffgroup = mean(pdiff,3);
sdgroup = std(pdiff,0,3)./sqrt(size(pdiff,3));
minn = 5;
maxn = 50;
dotsize  = squeeze(Scale(nansum(OutOfNum,3))*(maxn-minn)+minn);

xhd = linspace(0,135,1000);
xc = 0;
for x = 0:11.25:135;
    xc = xc+1;
    subplot(sp(1),sp(2),[1 6])
    errorbar(x,pdiffgroup(xc,2),sdgroup(xc,2),'.','Color','c','LineWidth',2,'MarkerEdgeColor','c','MarkerSize',dotsize(xc,2,1)) %before, CS-
    hold on
    % plot(xhd,PAL_Weibull(out.params1(2,:),xhd),'LineWidth',2,'Color','c')
    errorbar(x,pdiffgroup(xc,1),sdgroup(xc,1),'.','Color','r','LineWidth',2,'MarkerEdgeColor','r','MarkerSize',dotsize(xc,1,1)) %before, CS+
    % plot(xhd,PAL_Weibull(out.params1(1,:),xhd),'LineWidth',2,'Color','r')
    box off
    xlim([-10 145])
    set(gca,'XTick',[0 45 90 135],'YTick',0:.2:1,'FontSize',12)
    ylabel('p(different)','FontSize',12)
    xlabel('delta X [deg]')
    title('Discrimination pre','FontSize',14)
    ylim([0 1]);
    %     plot(params1(1,1),PAL_Weibull(params1(1,:),params1(1,1)),'k.','MarkerSize',20)
    
    subplot(sp(1),sp(2),[5 10])
    errorbar(x,pdiffgroup(xc,4),sdgroup(xc,4),'.','Color','c','LineWidth',2,'MarkerEdgeColor','c','MarkerSize',dotsize(xc,2,2)) %after, CS-
    hold on
    % plot(xhd,PAL_Weibull(out.params1(4,:),xhd),'LineWidth',2,'Color','c')
    errorbar(x,pdiffgroup(xc,3),sdgroup(xc,3),'.','Color','r','LineWidth',2,'MarkerEdgeColor','r','MarkerSize',dotsize(xc,1,2)) %after, CS+
    % plot(xhd,PAL_Weibull(out.params1(3,:),xhd),'LineWidth',2,'Color','r')
    box off
    xlim([-10 145])
    set(gca,'XTick',[0 45 90 135],'YTick',0:.2:1,'FontSize',12)
    ylabel('p(different)','FontSize',12)
    xlabel('delta X [deg]')
    title('Discrimination post','FontSize',14)
    ylim([0 1])
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histogram of alpha and kappa
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 12;

load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate2pmf_N54.mat')
ind = subs6;
sp = [2 4];
crit = 25; % for later noshifter vs shifter improvement
binsize = 8;

clf
subplot(sp(1),sp(2),1)%FWHM
hold off
b = linspace(0,180,binsize);
histf(mat(ind,18),b,'FaceColor','b','EdgeColor','none','FaceAlpha',.6);
hold on;
histf(mat(ind,19),b,'FaceColor','k','EdgeColor','none','FaceAlpha',.6);
axis square;box off
hold off
legend('Conditioning','Generalization','orientation','vertical','location','best')
ylabel('counts','FontSize',fontsize)
title('FWHM','FontSize',fontsize+2)
text(10,7.8e-3,'A','FontSize',fontsize+8);
xlim([0 200])
xlabel('deg','FontSize',fontsize)
set(gca,'YTick',[0 5 10],'FontSize',fontsize)

subplot(sp(1),sp(2),2)%MU
hold off
b = linspace(0,40,binsize-1);
histf(abs(mat(ind,15)),b,'FaceColor','b','EdgeColor','none','FaceAlpha',.6);
hold on;
histf(abs(mat(ind,16)),b,'FaceColor','k','EdgeColor','none','FaceAlpha',.6);
axis square;box off
hold off
ylabel('counts','FontSize',fontsize)
title('abs(Mu)','FontSize',fontsize+2)
hold off
xlim([-5 45])
xlabel('deg','FontSize',fontsize)
set(gca,'YTick',[0 5 10],'FontSize',fontsize)

%% second row - thresholds and improvement
fontsize = 12;
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat')

sp = [2 4];
crit = 22.5; % for later noshifter vs shifter improvement

ind = logical(abs(mat6(:,16))<= crit);

subplot(sp(1),sp(2),5)%alpha results
ylabel('threshold alpha (deg)','FontSize',fontsize)
b = bar([1 2 4 5],mean(mat6(:,1:4)));
title('discrimination thresholds','FontSize',fontsize+2)
SetFearGenBarColors(b,[1 0 0;.6 0 0;0 1 1;0 .6 .6]')
line([1 2],[72 72],'Color','k','LineWidth',2)
text(1,75,'***','FontSize',fontsize+4)
hold on;
e = errorbar([1 2 4 5],mean(mat6(:,1:4)),std(mat6(:,1:4))./sqrt(length(mat6)),'k.');
set(e,'LineWidth',2)
set(gca,'xtick',[1.5 4.5],'xticklabel',{'CS+' 'CS-'},'FontSize',fontsize)
set(e,'LineWidth',2)
axis square;box off;
xlim([0 6])
text(0,70,'C','FontSize',fontsize+8);

subplot(sp(1),sp(2),6)%csp improvement shifters vs no shifters
b = bar(1:2,[mean(mat6(ind,11)) mean(mat6(~ind,11))]);
set(b,'BarWidth',.5,'FaceColor','k')
title('discrimination improvement','FontSize',fontsize+2)
hold on;box off
e = errorbar(1:2,[mean(mat6(ind,11)) mean(mat6(~ind,11))],[std(mat6(ind,11))./sqrt(sum(ind)) std(mat6(~ind,11))./sqrt(sum(~ind))],'k.','LineWidth',1.5);
e = errorbar(1,mean(mat6(ind,11)),std(mat6(ind,11))./sqrt(sum(ind)),0,'w.','LineWidth',1.5);
axis square
xlim([0 3])
set(gca,'XTick',[1 2],'XTicklabel',{'centered' 'shifted'},'YTick',-10:10:20,'FontSize',fontsize)
ylim([-10 20])
ylabel('corrected improvement (deg)','FontSize',fontsize)
text(.5,78,'C','FontSize',20);
hold off;

subplot(sp(1),sp(2),2)%MU
l=line([crit crit],[0 10]);set(l,'Color','r','LineWidth',1.5);

%% correlation
% subplot(sp(1),sp(2),[3 4 7 8])
% text(20,180,'C','FontSize',20);
% hold on;
% plot(mat6(:,17),mat6(:,18),'bo','MarkerFaceColor','b','MarkerSize',10)
% l = lsline;set(l,'LineWidth',2)
% box off
% xlabel('initial alpha')
% ylabel('FWHM Cond')
% % regression
% X = mat6(:,17); %initial alpha
% Y = mat6(:,18); %FWHM cond
% mdl = LinearModel.fit(X,Y);
% plot(mdl)
% legend({'Data' 'Fit' 'Confidence Interval' '' 'Identity Line'},'location','best')
% legend boxoff
% set(gca,'XTick',20:20:110,'YTick',20:20:180)
% ylim([20 180])
% DrawIdentityLine(gca);
% xlim([20 110])
% title('Discrimination x Generalization', 'FontSize',14)
% set(gca,'FontSize',12)
% xlabel('Discrimination (alpha [deg])')
% ylabel('FWHM Conditioning')
% axis square
% box off
%% new way (bootstrap regression)
clear all
sp = [2 4];
fontsize = 12;
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat')
x  = mat6(:,17);
y  = mat6(:,18);
subplot(sp(1),sp(2),[3 4 7 8])
plot(x,y,'o','MarkerFaceColor','b')
hold on;
%
i = 1:length(x);
for n = 1:10000
    ii     = randsample(i,length(x),1);
    B(:,n) = [x(ii) ones(length(x),1)]\y(ii);
end
%
X2 = [linspace(min(x),max(x),100)' ones(100,1)];
plot(X2(:,1),X2*mean(B,2),'r','LineWidth',2)
Y = X2*B;
Y = Y';
YY = prctile(Y,[2.5 97.5])';
% plot(X2(:,1),YY(:,1),'r')
% plot(X2(:,1),YY(:,2),'r')
plot(linspace(min(x),max(x),100),linspace(min(x),max(x),100),'b','LineWidth',2)
xlim([min(x)-5 max(x)+5])
ylim([20 210])
%shaded area
x1 = X2(:,1)';
y1 = YY(:,1)';
x2 = X2(:,1)';
y2 = YY(:,2)';
fill([x1 fliplr(x2)],[y1 fliplr(y2)],'r','FaceAlpha',.5,'EdgeColor','none')
text(30,190,'B','FontSize',fontsize+8);
ylabel('FWHM Cond [deg]','FontSize',fontsize)
xlabel('Initial Alpha [deg]','FontSize',fontsize)
title('Discrimination x Fear Generalization','FontSize',fontsize+2)
box off

set(gca,'XTick',30:20:100,'YTick',20:40:200,'FontSize',fontsize)
hold off

%% bootstrap sigmoid
clear all
close all
nboot = 100000;
sp = [2 4];
fontsize = 12;
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\rate_and_pmf_N48.mat')
x0  = mat6(:,17);
x = x0 - min(x0);
y0  = mat6(:,18);
subplot(1,2,1);
plot(x0,y0,'o','MarkerFaceColor','b')
hold on;
y = y0./180;
subplot(1,2,2);plot(x,y,'o','MarkerFaceColor','k')
ndata = 1:length(x);
% run fit
clear params LL fit;
for n = 1:nboot
    fprintf('Bootstrap n = %g. \n',n)
    ii(:,n)     = randsample(ndata,length(x),1);
    %[params(:,n) LL(n)] = FitLogsCumNorm(x(ii(:,n)),y(ii(:,n)),1);
    [o] = FitLogistic_manually(x(ii(:,n)),y(ii(:,n)));
    fit(:,n)    = o.fitup;
    LL(:,n)     = o.Likelihood;
    params(:,n) = o.Est(1:4);
end
params1 = params;
params(1,:) = params(1,:) + min(x0);
%%
% prepare figure
X2 = linspace(min(x0),max(x0),100);
for n = 1:nboot
    Y(:,n) = PAL_CumulativeNormal(params(:,n),X2)';
    Y0(:,n) = (PAL_CumulativeNormal(params(:,n),X2)*180)';
end
Y0 = Y0';
YY = prctile(Y0,[2.5 97.5]);
plot(X2,prctile(Y0,50),'r','LineWidth',2);
%% plot
subplot(sp(1),sp(2),[3 4 7 8])
plot(x0,y0,'o','MarkerFaceColor','b')
hold on;
plot(X2,prctile(Y0,50),'r','LineWidth',2);
%shaded area
x1 = X2;
y1 = YY(1,:);
x2 = X2;
y2 = YY(2,:);
fill([x1 fliplr(x2)],[y1 fliplr(y2)],'r','FaceAlpha',.5,'EdgeColor','none')
xlim([min(x0)-10 max(x0)+10])
ylim([min(y0)-10 max(y0)+10])
axis square
box off
set(gca,'YTick',20:40:200,'XTick',20:20:100)
text(13,185,'B','FontSize',20);
%% illustration for new threshold .63 vs .5
guess = .2;
lapse = .1;

clf;
PlotTransparentLine(linspace(0,120,1000)',PAL_Weibull([45 3 .2 .1],linspace(0,120,1000))',.4,'r','LineWidth',3);
hold on;
% noise = randn(10,1)'*.03.*randsample([-1 1],10,true);
e = errorbar(linspace(0,100,10),PAL_Weibull([45 3 .2 .1],linspace(0,120,10))+noise,noise+.02,'ko','MarkerFaceColor','k','LineWidth',2);
ylim([0 1]);
xlim([-5 105]);
axis square
set(gca,'YTick',.2:.2:1,'XTick',0:20:100,'FontSize',14)
line([45 45],ylim)
line(xlim,repmat(PAL_Weibull([45 3 guess lapse],45),[2 1]))
line(repmat(PAL_Weibull([45 3 guess lapse],(1-guess-lapse)*.5+guess,'inverse'),[2 1]),ylim)
line(xlim,repmat((1-guess-lapse)*.5+guess,[2 1]))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% don't mention the peakshift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph alpha x discrimination (binning of M and SD)
clear all;p = Project;
subjects = intersect(p.subjects_600,find(sum(p.getMask('PMF'),2)==4));
subjects = intersect(subjects,find(p.getMask('RATE')));
g = Group(subjects);
mat =  g.parameterMat;
fontsize =  14;
binsize = 8;
sp = [1 4];
figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scatterplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
subplot(sp(1),sp(2),1)
plot(mat(:,17),mat(:,18),'bo','MarkerFaceColor','b')
l = lsline;
set(l,'LineWidth',2,'Color','b')
xlim([min(mat(:,17))-10 max(mat(:,17))+10])
ylim([min(mat(:,18))-10 max(mat(:,18))+10])
set(gca,'XTick',0:50:100,'YTick',0:50:150,'fontsize',fontsize)
axis square
box off
xlabel('initial threshold \alpha','fontsize',fontsize)
ylabel('FWHM Conditioning','fontsize',fontsize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bin perception, plot FWHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 600ms

bincol = [17];%which col for bins;
outc = 18; %which col for outcome (cond = 18, test=19)
p = prctile(mean(mat(:,bincol),2),linspace(0,100,4));
outcome_m = [];outcome_s=[];binalpha = [];binsd = [];
for n = 1:length(p)-1
    if n ==1
        i = (mean(mat(:,bincol),2) >= p(n)) & (mean(mat(:,bincol),2) <= p(n+1));
    else
        i = (mean(mat(:,bincol),2) >= p(n)) & (mean(mat(:,bincol),2) <= p(n+1));
    end
    sum(i)
    ii(:,n) = i;
    binalpha(n) = mean(mean(mat(i,[1 3]),2));
    binsd(n) = std(mean(mat(i,[1 3]),2));
    outcome_m(n) = mean(mean(mat(i,outc),2));
    outcome_s(n) = std(mean(mat(i,outc),2))./sqrt(sum(i));
end
subplot(sp(1),sp(2),2)
text(-.5,200,'B','FontSize',20);
bar(1:3,outcome_m,.6,'FaceColor',[.3 .3 .3],'EdgeColor','none');hold on;
plot(1:3,binalpha,'r*','Markersize',10)
% bar(1:3,binalpha,'FaceColor',[.7 .7 .7],'EdgeColor','none')%bar(1:3,[2.59 .995 .147],'FaceColor',[.3 .3 .3])
errorbar(1:3,outcome_m,outcome_s,'k.','LineWidth',2)
set(gca,'XTickLabel',{'good','moderate','weak'},'fontsize',fontsize-2)
set(gca,'YTick',0:100:200,'fontsize',fontsize)
xlabel('discrimination performance','fontsize',fontsize)
ylabel('M FWHM Conditioning','fontsize',fontsize)
axis square
box off
xlim([0 4])
hold off
subplot(sp(1),sp(2),3)
bar(1:3,outcome_s,.6,'FaceColor',[.5 .5 .5],'EdgeColor','none');hold on;
set(gca,'XTickLabel',{'good','moderate','weak'},'fontsize',fontsize-2)
set(gca,'YTick',[0:10:20],'fontsize',fontsize)
xlabel('discrimination performance','fontsize',fontsize)
ylabel('SEM FWHM Conditioning','fontsize',fontsize)
ylim([0 21])
axis square
box off
xlim([0 4])
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% threshold bars, before after
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
sp = [1 4];
p = Project;
g = Group(intersect(p.subjects_600,find(sum(p.getMask('PMF'),2)==4)));
mat = g.parameterMat;
subplot(sp(1),sp(2),4)%alpha results
ylabel('threshold alpha (deg)','FontSize',fontsize)
b = bar([1 2 4 5],mean(mat(:,[1 2 3 4])));
SetFearGenBarColors(b,[1 0 0;.6 0 0;0 1 1;0 .6 .6]')
line([1 2],[mean(mat(:,1))+10 mean(mat(:,1))+10],'Color','k','LineWidth',2)
text(1,mean(mat(:,1))+13,'***','FontSize',fontsize+4)
hold on;
e = errorbar([1 2 4 5],mean(mat(:,[1 2 3 4])),std(mat(:,[1 2 3 4]))./sqrt(length(g.ids)),'k.');
set(e,'LineWidth',2)
ylim([35 80])
set(gca,'xtick',[1.5 4.5],'xticklabel',{'CS+' 'CS-'},'FontSize',fontsize,'YTick',[20 40 60 80])
ylabel('perceptual threshold \alpha [deg]','Fontsize',fontsize)
set(e,'LineWidth',2)
axis square;box off;
xlim([0 6])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 2 - Task Demands x Subject Idiosyncrasy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
p           = Project;
mask = p.getMask('ET_discr');
subjects = intersect([Project.subjects_1500],find(mask));
mask = p.getMask('ET_feargen');
subjects = intersect(subjects,find(mask));
fix = Fixmat(subjects,1:5);
v = [];
c = 0;
for sub = [6 13 10]
    for ph = 1:5
        c = c+1;
        v{c} = {'subject' sub 'phase' ph 'deltacsp' [0 180 18000]};
    end
    c=c+1;
    v{c} = {'subject' sub 'deltacsp' [ 0 180 18000]};
end

for ph = 1:5
    c = c+1;
    v{c} = {'phase' ph 'deltacsp' [0 180 18000]};
end
c = c+1;
v{c} = {'deltacsp' [0 180 18000]};
fix.getmaps(v{:});
fix.plot('linear',[4 6])


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
fix = Fixmat(subjects, 1);
fix.getsubmaps;
% fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmaps = fix.vectorize_maps;

%normal correlation
for n = 1:length(subjmaps)
    [r(n),pval(n)] = corr(subjmaps(n,:)',mean(mat(:,[1 3]),2),'type','Spearman');
end
ss = [sqrt(length(r)) sqrt(length(r))];
clf
imagesc(reshape(r,ss),[-1 1])
for n = find(pval<0.001);[y x] = ind2sub(ss,n);text(x,y,'x','FontSize',8);end
fix.maps = reshape(r,ss);
fix.plot
% correct it by mean fixmap
fix.getsubmaps;
% fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = Scale(mean(fix.maps,3));
rc       = r.*fix.maps(:)';
fix.maps = reshape(rc,ss);
fix.plot
% get the numbers
[mini,i] = min(rc)
r(i)
pval(i)
[y,x] = ind2sub([50 50], i)
% -> minimal pixel for rc doesn't get significant in the original
% correlation..

%% weight fixmaps, then correlate
for n = 1:25;subjmaps(:,n) = fix.maps(:,n).*scaling;end
fix.getsubmaps;
fix.maps = fix.vectorize_maps;
scaling = mean(fix.maps,2);
for n = 1:25;subjmaps(:,n) = fix.maps(:,n).*scaling;end

for n = 1:length(subjmaps)
    [r(n),pval(n)] = corr(subjmaps(n,:)',mean(mat(:,[1 3]),2),'type','Spearman');
end
rb(pval>.01)=NaN;
fix.maps = reshape(rb,[500 500]);
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
% fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmaps = fix.vectorize_maps;

for n = 1:25; param(n,1) = vM2FWHM(mat(n,13));end
%normal correlation
for n = 1:length(subjmaps)
    [r(n),pval(n)] = corr(subjmaps(n,:)',param,'type','Spearman');
end
clf
ss = [sqrt(length(r)) sqrt(length(r))];
imagesc(reshape(r,ss),[-.5 .5])
for n = find(pval<0.05);[y x] = ind2sub(ss,n);text(x,y,'x','FontSize',8);end
figure
fix.maps = reshape(r,ss);
fix.plot
% correct it by mean fixmap
fix.getsubmaps;
% fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = Scale(mean(fix.maps,3));
rc       = r.*fix.maps(:)';
fix.maps = reshape(rc,ss);
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
subjects        = intersect(find(mask),Project.subjects_1500);
fix             = Fixmat(subjects,[1 2 3 4 5]);

%% plot and save single subject fixations maps RSA
tsub = length(unique(fix.subject));
cormat = nan(16,16,tsub);
pval = nan(16,16,tsub);
subc            = 0;
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
    fix.maps(:,:,1:8)     = maps(:,:,1:8) - repmat(mean(maps(:,:,1:8),3),[1 1 8]);%take out the average
    fix.maps(:,:,9:end)   = maps(:,:,9:end) - repmat(mean(maps(:,:,9:end),3),[1 1 8]);%take out the average
    [cormat(:,:,subc) pval(:,:,subc)] = fix.corr;
end

%% RSA plot but fisherz transformed
cormatz = fisherz_inverse(median(fisherz(cormat),3));
cormatz = CancelDiagonals(cormatz,0);
figure
imagesc(cormatz)
axis square;colorbar;
set(gca,'fontsize',15)
axis off
%% RSA by single fixations
nfix = 1:4;
tsub = length(unique(fix.subject));
tfix = length(fix);
cormat_t = nan(8,8,tsub,tfix);
cormat_b = nan(8,8,tsub,tfix);
subc   = 0;

for subject = unique(fix.subject);
    for nf = nfix(:)';
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

%% SVM of CS+ vs CS-
r4 = load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\findingparams\SVM_simulations_ph4_NEV18_fwhm5.mat','AVE','HP','result','R');
r2 = load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\svm_analysis\findingparams\SVM_simulations_ph2_NEV18_fwhm5.mat','AVE','HP','result','R');
fix = Fixmat(6,3);
fix.maps = cat(3,r2.AVE,r4.AVE);
fix.plot
cbh = colorbar;
set(cbh,'YTick',-.04:.02:.04,'FontSize',12)
%% SVM classified as CSP
close all
x = -135:45:180;
xhd = linspace(-165,210,100);
figure
subplot(1,2,1)
m = mean(mean(r2,2),3);
sd = std(mean(r2,2),0,3);
sem = sd./sqrt(size(r2,3));
b = bar(x,m);
SetFearGenBarColors(b);
hold on;
errorbar(x,m,sem,'k.','LineWidth',2);
plot(xhd,t2.fit_results.fitfun(xhd,t2.fit_results.params(1:4)),'k','LineWidth',2)
ylabel('Classified as CS+','FontSize',12)
title('Free Viewing','fontsize',14)

subplot(1,2,2);
m = mean(mean(r4,2),3);
sd = std(mean(r4,2),0,3);
sem = sd./sqrt(size(r4,3));
b = bar(x,m);
SetFearGenBarColors(b);
hold on;
errorbar(x,m,sem,'k.','LineWidth',2);
plot(xhd,t4.fit_results.fitfun(xhd,t4.fit_results.params(1:4)),'k','LineWidth',2)
title('Generalization','fontsize',14)

EqualizeSubPlotYlim(gcf);
for n = 1:2;
    subplot(1,2,n);
    set(gca,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'},'YTick',.3:.1:.7,'FontSize',12)
    ylim([.3 .7])
    box off
    axis square
    xlim([-180 225])
    line([-180 225],[.5 .5],'LineStyle',':','LineWidth',2,'Color','k')%.5 8.5 for 1:8 xscale
end
SaveFigure('classified_as_csp.eps')