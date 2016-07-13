%general response to 8 faces, for 1500 and 600ms
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scrmat.mat')
p = Project;
mask = p.getMask('SCR');
any(~ismember(ss,find(sum(mask,2)==3)))%there's none, so we didn't need to exclude any of the 54 people.

ind15 = ismember(ss,Project.subjects_1500);%n = 25
ind6 = ismember(ss,Project.subjects_600);% n = 29
scr = m;
scr15 = mean(scr(ind15,:));
[h,pval,ci,stats] = ttest(scr(ind15,3),scr(ind15,4));
scr6 = mean(scr(ind6,:));
subplot(1,2,1)
b = bar(scr15);
SetFearGenBarColors(b);
subplot(1,2,2)
b = bar(scr6);
SetFearGenBarColors(b);
EqualizeSubPlotYlim(gcf)

%% plot scr to 8 faces BCT
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_600BCT.mat','scr_bars','subjects')
ylims = [0 1.1];
Yticks = [0 0.4 0.8 1.2];
% or
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_1500BCT.mat','subjects','scr_bars')
Yticks = [0 0.3 0.6];
ylims = [0 0.7];
%% plot
clf
for n = [1 3]
subplot(1,3,n);
b = bar(1:8,mean(scr_bars(:,1:8,n)));
hold on;
ylim(ylims)
xlim([0 9])
e = errorbar(1:8,mean(scr_bars(:,1:8,n)),std(scr_bars(:,1:8,n))./sqrt(size(scr_bars,1)),'k.');
set(gca,'XTick',1:8,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',Yticks)
SetFearGenBarColors(b)
set(e,'LineWidth',2,'Color','k')
axis square
box off
end
subplot(1,3,1);ylabel('SCR [muS]')
%cond special treatment
subplot(1,3,2);
b(1) = bar(4,mean(scr_bars(:,4,2)));
hold on;
e(1) = errorbar(4,mean(scr_bars(:,4,2)),std(scr_bars(:,4,2))./sqrt(size(scr_bars,1)),'k.');
ylim(ylims)
xlim([0 9])
set(gca,'XTick',1:8,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',Yticks);
b(2) =  bar(8,mean(scr_bars(:,8,2)));
e(2) = errorbar(8,mean(scr_bars(:,8,2)),std(scr_bars(:,8,2))./sqrt(size(scr_bars,1)));
set(b(1),'FaceColor',[1 0 0],'EdgeColor','w');
set(b(2),'FaceColor','c','EdgeColor','w');
set(e,'LineWidth',2,'Color','k');
axis square
box off
try
    for n = 1:size(scr_bars,3)
        subplot(1,3,n)
        line(xlim,repmat(mean(scr_bars(:,9,n)),[2 1]),'Color','k','LineWidth',1.3,'LineStyle',':')
    end
end

t = supertitle('SCR 600 ms');set(t,'FontSize',14)


%% is that connected to improvement in pmf or ratings?
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scrmat.mat')
p = Project;
%% pmf
mask = p.getMask('PMF');mask = find(sum(mask,2)==4);
subjects = intersect(mask,[Project.subjects_600 Project.subjects_1500]);
subjects = intersect(subjects,ss);
sub15 = ismember(subjects,Project.subjects_1500);
sub6  = ismember(subjects,Project.subjects_600);
g = Group(subjects);
scr_response = m(ismember(ss,subjects),4)-m(ismember(ss,subjects),8);
[mat tags] = g.parameterMat;
%% rate
mask = p.getMask('RATE');mask = find(mask);
subjects = intersect(mask,[Project.subjects_600 Project.subjects_1500]);
subjects = intersect(subjects,ss);
sub15 = ismember(subjects,Project.subjects_1500);
sub6  = ismember(subjects,Project.subjects_600);
g = Group(subjects);
g.getSI(8);
scr_response = m(ismember(ss,subjects),4)-m(ismember(ss,subjects),8);
[mat tags] = g.parameterMat;
clear g;

%% does SCR correlate with other behavioral parameters?
param = mat(:,14);
[r,pval] = corr(param,scr_response)
% look at durations separately
clf
plot(scr_response,param,'k.');hold on;l=lsline;set(l,'LineWidth',2)
plot(scr_response(sub15),param(sub15),'b.','MarkerSize',20);hold on;l=lsline;set(l,'LineWidth',2)
plot(scr_response(sub6),param(sub6),'r.','MarkerSize',20);l=lsline;set(l,'LineWidth',2)
a(1)=ylabel('CSP impr');
a(2)=xlabel('SCR (CSP-CSN)');
set(a,'FontSize',15);
%% plot single subject scr tunings.
figure;
c=0;
for n = 1:length(scr_response)
    c=c+1;
subplot(7,8,c)
b = bar(m(n,:));
SetFearGenBarColors(b);
end
EqualizeSubPlotYlim(gcf)

%% all durations
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\pmf_data_N51.mat')
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scrmat.mat')
valid = intersect(subjects,ss);
%optionally restrict duration
valid = intersect(valid, Project.subjects_600);
scr_response = m(:,4)-m(:,8);
%or
%scr_response = mean(m,2);
[r,pval] = corr(mat(ismember(subjects,valid),11),scr_response(ismember(ss,valid)))
plot(mat(ismember(subjects,valid),11),scr_response(ismember(ss,valid)),'bo')
%only for people with higher SCR for CSP than CSN
ss1 = ss(find(m(:,4)>m(:,8)));
valid = intersect(subjects,ss1)
[r,pval] = corr(mat(ismember(subjects,valid),11),scr_response(ismember(ss,valid)));
plot(mat(ismember(subjects,valid),11),scr_response(ismember(ss,valid)),'bo')
%general SCR level?
scr = mean(m,2);

%% interaction with fear tunings?
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scrmat.mat')
p = Project;
mask = p.getMask('RATE');
mask = intersect(find(mask),[Project.subjects_1500 Project.subjects_600]);% N = 54
valid = intersect(mask,ss);
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\ratings_mises.mat')
scr_response = mean(m,2);
[r,pval] = corrcoef([mat(ismember(mask,valid),12) mat(ismember(mask,valid),13) mat(ismember(mask,valid),14) ...
    mat(ismember(mask,valid),15) mat(ismember(mask,valid),16)...
    scr_response(ismember(ss,valid))])


%% conditioning manipulation check
p = Project;
mask = find(sum(p.getMask('SCR'),2)==3);
subjects = intersect(mask,Project.subjects_1500);

%% SCR graphs, not bars
scr_data = NaN(800,9,length(subjects));
sc=0;
for sub = subjects(:)'
    fprintf('Working on subject %d .. \n',sub)
    sc=sc+1;
    s = Subject(sub);
    s.scr.cut(s.scr.findphase('test$'));
    s.scr.run_ledalab;
    scr_data(:,:,sc) = s.scr.ledalab.mean(1:800,1:9);
end
%% collect SCR bars
phasenames = {'base$' 'cond$' 'test$'};
phaseconds = [9 3 9];
scr_bars = NaN(length(subjects),9,3);
for ph = 1:3
    sc = 0;
    for sub = subjects(:)'
        fprintf('Working on phase %d, subject %d .. \n',ph,sub)
        sc=sc+1;
        s = Subject(sub);
        s.scr.cut(s.scr.findphase(phasenames{ph}));
        s.scr.run_ledalab;        
        s.scr.plot_tuning_ledalab(1:phaseconds(ph));
        if ismember(ph,[1 3])
            scr_bars(sc,:,ph) = s.scr.fear_tuning;
        else
            scr_bars(sc,[4 8 9],ph) = s.scr.fear_tuning;
        end
        close all
        clear s
    end
end
%% SCR graphs, not bars
scr_data = NaN(800,9,length(subjects),3);
phasenames = {'base$' 'cond$' 'test$'};
phaseconds = [9 3 9];
sc = 0;
for ph = 1:3
    for sub = subjects(:)'
        fprintf('Working on subject %d .. \n',sub)
        sc=sc+1;
        s = Subject(sub);
        s.scr.cut(s.scr.findphase(phasenames{ph}));
        s.scr.run_ledalab;
        if ismember(ph,[1 3])
            scr_data(:,:,sc,ph) = s.scr.ledalab.mean(1:800,1:9);
        else
            scr_data(:,[4 8 9],sc,ph) = s.scr.ledalab.mean(1:800,1:9);
        end
    end
end
%% SCR graphs, not bars
scr_data0 = NaN(800,9,length(subjects));
sc=0;
for sub = subjects(:)'
    fprintf('Working on subject %d .. \n',sub)
    sc=sc+1;
    s = Subject(sub);
    s.scr.cut(s.scr.findphase('base$'));
    s.scr.run_ledalab;
    scr_data0(:,:,sc) = s.scr.ledalab.mean(1:800,1:9);
    
end
scr_data2 = NaN(800,9,length(subjects));
sc=0;
for sub = subjects(:)'
    fprintf('Working on subject %d .. \n',sub)
    sc=sc+1;
    s = Subject(sub);
    s.scr.cut(s.scr.findphase('test$'));
    s.scr.run_ledalab;
    scr_data2(:,:,sc) = s.scr.ledalab.mean(1:800,1:9);
end
scr_data1 = NaN(800,9,length(subjects));
sc=0;
for sub = subjects(:)'
    fprintf('Working on subject %d .. \n',sub)
    sc=sc+1;
    s = Subject(sub);
    s.scr.cut(s.scr.findphase('cond$'));
    s.scr.run_ledalab;
    scr_data1(:,[4 8 9],sc) = s.scr.ledalab.mean(1:800,1:3);    
end

scr_data = cat(4,scr_data0,scr_data1,scr_data2);

%check window
figure;
SetFearGenColors;
plot(mean(scr_data(:,1:8,:,3),3),'LineWidth',2);hold on;
plot(mean(scr_data(:,9,:,3),3),'k','LineWidth',2);
%cut part that's interesting for mean bars
scr_data_cut = scr_data(250:550,:,:,:);%take window of 200:500 for 600, and 250:550 for 1500;
data = squeeze(mean(scr_data_cut(:,1:8,:,:)));
dataz = nanzscore(data);
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_600.mat','scr_data','subjects','scr_data_cut','data','dataz');
%plot single tunings
for n=1:length(subjects);subplot(5,6,n);b=bar(dataz(:,n,3));SetFearGenBarColors(b);end
%plot group
for n = 1:3;subplot(1,3,n);b=bar(mean(dataz(:,:,n),2));SetFearGenBarColors(b);hold on;errorbar(mean(dataz(:,:,n),2),std(data(:,:,n),0,2)./sqrt(size(dataz,2)),'k.','LineWidth',2');axis square;box off;end
EqualizeSubPlotYlim(gcf);
%% plot single fits
nsub = 29;
x = -135:45:180;
for n = 1:nsub
    clear data;
    data.x = x;
    data.y = scr_bars(1:8,n,3)';
    data.ids = n;
    tuning = Tuning(data);
    tuning.SingleSubjectFit(8);
    fit{n} = tuning.singlesubject{1};
    signif(n) = fit{n}.pval > -log10(0.05);
end
%plot the bars plus curves where significant
for n = 1:nsub
subplot(floor(sqrt(nsub)),ceil(sqrt(nsub)),n);
b = bar(x,scr_bars(1:8,n,3));
SetFearGenBarColors(b);
set(gca,'XTick',[0 180],'XTickLabel',{'CS+','CS-'})
hold on;
x_HD = linspace(min(x),max(x),1000);
if signif(n)
    plot(x_HD,fit{n}.fitfun(x_HD,fit{n}.Est),'k-','linewidth',2);
end
end
%groupfit the tuned people
data.x = repmat(x,[sum(signif) 1]);
data.y = scr_bars(logical(signif),1:8,3);
data.ids = subjects(logical(signif));

% % EqualizeSubPlotYlim(gcf)
% for n = 1:length(g.ids);valid(n) = g.tunings.scr.singlesubject{n}.pval > -log10(0.05);end
% for n = find(valid)
% subplot(5,6,n)
% bla = ylim;
% text(135,bla(2)*0.8,'X','FontSize',20)
% end
%% make Groupfit and plot it
g.tunings.scr.GroupFit(8)
% ylims = [0 0.55];%1500
ylims = [0 0.8]; %600
% Yticks = 0:0.1:0.5;%1500
Yticks = [0 0.2 0.4 0.6 0.8];%600

clf
n = 1;
x_HD = linspace(min(g.tunings.scr.x(n,:)),max(g.tunings.scr.x(n,:)),1000);
b = bar(-135:45:180,mean(g.tunings.scr.y));
axis square
xlim([-180 225])
ylim(ylims)
hold on;
SetFearGenBarColors(b);
plot(x_HD,g.tunings.scr.groupfit.fitfun(linspace(-135,180,1000),g.tunings.scr.groupfit.Est),'k','LineWidth',2)
errorbar(-135:45:180,mean(g.tunings.scr.y),std(g.tunings.scr.y)./sqrt(length(g.ids)),'k.','LineWidth',2)
set(gca,'XTicklabel',-135:45:180,'YTick',Yticks)
title('Groupfit Von Mises SCR 600 ms');
%compare it to mean(single subject fit) Curve
plot(x_HD,g.tunings.scr.groupfit.fitfun(linspace(-135,180,1000),mean(g.tunings.scr.params)),'b','LineWidth',2)

%% 
p = Project;
mask = find(sum(p.getMask('SCR'),2)==3);
subs_scr = intersect(mask,Project.subjects_1500);

mask = find(sum(p.getMask('PMF'),2)==4);
subs_pmf = intersect(mask,Project.subjects_1500);

subs = intersect(subs_pmf,subs_scr);
g = Group(subs); mat = g.parameterMat;

crit = (scr_bars(:,4,3)-scr_bars(:,8,3))./(scr_bars(:,4,3)+scr_bars(:,8,3));
crit = params(:,1);
[r,pval]=corr(mat(:,11),crit(ismember(subs_scr,subs)))

%% just putting this here
                i                    = (self.ledalab.x(:,1) >= 2.5)&(self.ledalab.x(:,1) <= 5.5);
                self.fear_tuning     = mean(self.ledalab.mean(i,:));%take out the average in that window
                self.fear_tuning = self.fear_tuning(:,conds);