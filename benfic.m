clear all;
p = Project;
subjects = intersect(p.subjects_1500,find(sum(p.getMask('PMF'),2)==4));
subjects = intersect(subjects,find(p.getMask('RATE')));
subjects = intersect(find(p.getMask('ET_feargen').*p.getMask('ET_discr')),subjects);
%subjects = intersect(subjects,find(sum(p.getMask('SCR'),2)==3));

fix  = Fixmat(subjects,[1 3 4 5]);

g = Group(subjects);
mat = g.parameterMat;


% fix.kernel_fwhm = 72;
M = [];
for ns = unique(fix.subject)
    query = {'subject' ns 'deltacsp' [0 180 18000]};
    fix.getmaps(query);
    M     = cat(3,M,fix.maps);
end
%

% load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_beneficial_locations.mat')
a    = mat(:,17);
b    = mat(:,19);
% c    = [DD(12,:)-DD(16,:)]';
d    = mat(:,11);
e    = mat(:,20);
% f    = [DD(20,:)-DD(24,:)]';
% % exclude sub 46 bc no usable eye data
% a(8) = [];
% b(8) = [];
%
titletags = {'alpha pre' 'FWHM test' 'SCR Cond CS+>CS-' 'CS+ corr impr' 'SI' 'SCR Test CS+>CS-'};

fix.maps = M;
type = 'Spearman';
[r_a, p_a]     = corr(fix.vectorize_maps',a,'Type',type);
[r_b, p_b]     = corr(fix.vectorize_maps',b,'Type',type);
% [r_c, p_c]     = corr(fix.vectorize_maps',c,'Type',type);
[r_d, p_d]     = corr(fix.vectorize_maps',d,'Type',type);

[r_e, p_e]     = corr(fix.vectorize_maps',e,'Type',type);
% [r_f, p_f]     = corr(fix.vectorize_maps',f,'Type',type);
% crit  = .10;
% r_a(p_a>crit) = nan;
% r_b(p_b>crit) = nan;
% r_c(p_c>crit) = nan;
% r_d(p_d>crit) = nan;
% r_e(p_e>crit) = nan;
% r_f(p_f>crit) = nan;


fix.maps = reshape([r_a,r_b,r_d,r_e],[size(M,1) size(M,2) 4]);
figure;fix.plot;%title(mat2str(tags{n}),'interpreter','none');
for n=1:6
    subplot(2,3,n)
    t=title(titletags{n});
    set(t,'FontSize',14)
end
% drawnow
%
% subplot(2,3,1);title('initial alpha');
% subplot(2,3,2);title('FWHM test');
% subplot(2,3,3);title('SCR Cond');
% subplot(2,3,4);title('alpha corr.improvement');
% subplot(2,3,5);title('FWHM SI');
% subplot(2,3,6);title('SCR Test');

%% same with t-test on median split
crit = mat(:,18);
good = crit<=median(crit);
t = nan(250000,1);
% FM = fix.vectorize_maps;
for n = 1:length(fix.vectorize_maps)
    [h,p,ci,stats] = ttest2(FM(n,good),FM(n,~good));
    t(n) = stats.tstat;
end
%% only phase for which parameter is specific
clear all;
p = Project;
subjects = intersect(p.subjects_1500,find(sum(p.getMask('PMF'),2)==4));
subjects = intersect(subjects,find(p.getMask('RATE')));
subjects = intersect(find(p.getMask('ET_feargen').*p.getMask('ET_discr')),subjects);
subjects = intersect(subjects,find(sum(p.getMask('SCR'),2)==3));

fix  = Fixmat(subjects,[1 3 4 5]);
g    = Group(subjects);
mat  = g.parameterMat;

param(:,1) = mat(:,17);
param(:,2) = mat(:,19);
param(:,3) = mat(:,11);
param(:,4) = mat(:,20);



phaseind = [1 4 1 4];
for np = 1:4
    fprintf('Computing correlations %g/%g. \n',np,size(param,2))
    M = [];
    for ns = unique(fix.subject)
        query = {'subject' ns 'deltacsp' [0 180 18000] 'phase' phaseind(np)};
        fix.getmaps(query);
        M     = cat(3,M,fix.maps);
    end
    fix.maps = M;
     [r(:,np) pval(:,np)] = corr(fix.vectorize_maps',param(:,np),'Type','Spearman');
end
fix.maps = reshape(r,[500 500 size(param,2)]);
fix.plot
%

%% link found beneficial maps to behavior
subjmaps = reshape(M,[500*500 22]);
disc = r_a;
disc(disc>=0) = 0;
disc(isnan(disc)) = 0;
disc = -disc;
subjload = disc'*subjmaps;
[rho,pval] = corr(mat(:,17),subjload')


%% bootstrap with rank correlation

clear all;
p        = Project;
subjects = intersect(p.subjects_1500,find(sum(p.getMask('PMF'),2)==4));
subjects = intersect(subjects,find(p.getMask('RATE')));
subjects = intersect(find(p.getMask('ET_feargen').*p.getMask('ET_discr')),subjects);

fix  = Fixmat(subjects,[1 3 4 5]);
g    = Group(subjects);
mat  = g.parameterMat;

% fix.kernel_fwhm = 72;
M = [];
for ns = unique(fix.subject)
    query = {'subject' ns 'deltacsp' [0 180 18000]};
    fix.getmaps(query);
    M     = cat(3,M,fix.maps);
end
%

% load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_beneficial_locations.mat')
behav    = mat(:,[17 19 11 20]);
fix.maps = M;
superfun = @(x,y) corr(x,y,'type','pearson');

R = [];
for nparam = 1:size(behav,2)    
    fprintf('Param: %03d\n',nparam);
    tic
    R(:,:,nparam) = bootstrp(1000,superfun,fix.vectorize_maps',behav(:,nparam));
    toc
end

%% 
clear all;
p = Project;
subjects15 = intersect(p.subjects_1500,find(sum(p.getMask('PMF'),2)==4));
subjects15 = intersect(subjects15,find(p.getMask('RATE')));
subjects15 = intersect(find(p.getMask('ET_feargen').*p.getMask('ET_discr')),subjects15);
subjects6 = intersect(p.subjects_600,find(sum(p.getMask('PMF'),2)==4));
subjects6 = intersect(subjects6,find(p.getMask('RATE')));
subjects6 = intersect(find(p.getMask('ET_feargen').*p.getMask('ET_discr')),subjects6);

fix  = Fixmat([subjects15' subjects6'],[1 3 4 5]);

g15 = Group(subjects15);
mat15 = g15.parameterMat;
g6 = Group(subjects6);
mat6 = g6.parameterMat;

% fix.kernel_fwhm = 72;
M1 = [];
for ns = subjects15(:)'
    query = {'subject' ns 'deltacsp' [0 180 18000]};
    fix.getmaps(query);
    M1     = cat(3,M1,fix.maps);
end
M2 = [];
for ns = subjects15(:)'
    query = {'subject' ns 'deltacsp' [0 180 18000]};
    fix.getmaps(query);
    M2     = cat(3,M2,fix.maps);
end
M3 = [];
for ns = subjects6(:)'
    query = {'subject' ns 'deltacsp' [0 180 18000]};
    fix.getmaps(query);
    M3     = cat(3,M3,fix.maps);
end
M4 = [];
for ns = subjects6(:)'
    query = {'subject' ns 'deltacsp' [0 180 18000]};
    fix.getmaps(query);
    M4     = cat(3,M4,fix.maps);
end
%
%

% load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_beneficial_locations.mat')
a    = mat15(:,17);
b    = mat15(:,18);
% c    = [DD(12,:)-DD(16,:)]';
d    = mat6(:,17);
e    = mat6(:,18);
% f    = [DD(20,:)-DD(24,:)]';
% % exclude sub 46 bc no usable eye data
% a(8) = [];
% b(8) = [];
%

M1 = reshape(M1,[500*500,22]);
M2 = reshape(M2,[500*500,22]);
M3 = reshape(M3,[500*500,24]);
M4 = reshape(M4,[500*500,24]);

type = 'Pearson';
[r_a, p_a]     = corr(M1',a,'Type',type);
[r_b, p_b]     = corr(M2',b,'Type',type);
% [r_c, p_c]     = corr(fix.vectorize_maps',c,'Type',type);
[r_d, p_d]     = corr(M3',d,'Type',type);

[r_e, p_e]     = corr(M4',e,'Type',type);
% [r_f, p_f]     = corr(fix.vectorize_maps',f,'Type',type);
% crit  = .10;
% r_a(p_a>crit) = nan;
% r_b(p_b>crit) = nan;
% r_c(p_c>crit) = nan;
% r_d(p_d>crit) = nan;
% r_e(p_e>crit) = nan;
% r_f(p_f>crit) = nan;

fix.maps = reshape([r_a,r_b,r_d,r_e],[500 500 4]);
figure;fix.plot;%title(mat2str(tags{n}),'interpreter','none');
for n=1:6
    subplot(2,3,n)
    t=title(titletags{n});
    set(t,'FontSize',14)
end
% drawnow
%
% subplot(2,3,1);title('initial alpha');
% subplot(2,3,2);title('FWHM test');
% subplot(2,3,3);title('SCR Cond');
% subplot(2,3,4);title('alpha corr.improvement');
% subplot(2,3,5);title('FWHM SI');
% subplot(2,3,6);title('SCR Test');

