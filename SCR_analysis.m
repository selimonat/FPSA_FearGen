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
b=bar(scr15);
SetFearGenBarColors(b);
subplot(1,2,2)
b=bar(scr6);
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
subjects = intersect(mask,Project.subjects_600);

phasenames = {'base$' 'cond$' 'test$'};
phaseconds = [9 3 9];
bdnf_problems = [2 13 18 31 32 39 50];
%%

for ph = 1:3
    sc = 0;
    for sub = subjects(~ismember(subjects,bdnf_problems))
        fprintf('Working on phase %d, subject %d .. \n',ph,sub)
        sc=sc+1;
        s = Subject(sub);
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

clear all
p = Project;
subjects =

phasenames = {'base$' 'cond$' 'test$'};
phaseconds = [8 2 8];
bdnf_problems = [2 13 18 31 32 39 50];


%% 


