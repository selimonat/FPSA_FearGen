%% get the data first;
addpath('/Users/onat/Documents/Code/Matlab/fancycarp/')
clear all;
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
fix             = Fixmat(subjects,[1 2 3 4 5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%MEGA SUBJECT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test how different fwhm parameters change the fixmaps and corrmats
v = [];
c = 0;
for ph = [2 4]
    for cond = -135:45:180
        c = c+1;
        v{c} = {'phase', ph, 'deltacsp' cond};
    end
end
% plot and save some fixation maps
for k = 37;%%Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
    fix.kernel_fwhm = k;
    fix.getmaps(v{:});
    
%     fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:16),3),[1 1 16]);%correct for baseline    
%     fix.maps    = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
    fix.maps    = fix.maps(:,:,1:8) - repmat(mean(fix.maps(:,:,1:8),3),[1 1 8]);%take out the average
    fix.plot;
    supertitle(mat2str(k),1);
%     SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/MegaSubject/AllFixations/unitize/1500/nocorrection/Maps_Scale_%04.4g.png',k));    
    out= fix.corr;
    figure(1000);clf
    imagesc(out);colorbar
%     SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/MegaSubject/AllFixations/unitize/1500/nocorrection/CorrMat_Scale_%04.4g.png',k));
    title(mat2str(k));        
end
%% For a given fwhm, compute fixmaps and cormat for different fixation indices.
fix.kernel_fwhm = 25;
for nfix = [1 2 3 4 5 6];%there are about 200 7th fixations (not enough)
    %creaete the query cell
    v = [];
    c = 0;
    for ph = [2 4]
        for cond = -135:45:180
            c = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'fix' nfix(:)};
        end
    end
    % plot and save fixation maps       
    fix.getmaps(v{:});
    fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
%     fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
    fix.plot;
    supertitle(mat2str(k),1);
%     SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/MegaSubject/SingleFixations/1500/Maps_nFix_%02d.png',nfix(1),nfix(2)));
    out = fix.corr;
    figure(10001);
    imagesc(out(9:end,9:16));colorbar    
%     SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/MegaSubject/SingleFixations/1500/Corrmat_nFix_%02d.png',nfix(1),nfix(2)));
    title([mat2str(fix.kernel_fwhm) '-' mat2str(nfix)] );
    figure(10002)
    bla(nfix,:) = out(12,9:end);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%SINGLE SUBJECT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All the above is done by considering all subjects as one single megasubject i.e. the analysis
% is not done individually.
% compute single subject correlation matrices
%% plot and save single subject fixations maps RSA
tsub = length(unique(fix.subject));
out  = [];
% for k = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
fix.kernel_fwhm = 25;
fix.kernel_fwhm = 37*0.8;
cormat_t = nan(8,8,tsub);
cormat_b = nan(8,8,tsub);
pval_t = nan(8,8,tsub);
pval_b = nan(8,8,tsub);
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
figure
subplot(1,2,1);
imagesc(median(cormat_b,3),[-.1 .1]);
axis square;colorbar;
set(gca,'fontsize',15)
axis off
subplot(1,2,2);
imagesc(median(cormat_t,3),[-.1 .1])
axis square;colorbar
set(gca,'fontsize',15)
axis off

%% RSA plot but fisherz transformed
for n = 1:27
    cormat_tz(:,:,n) = reshape(fisherz(cormat_t(:,:,n)),[8 8 1]);
    cormat_bz(:,:,n) = reshape(fisherz(cormat_b(:,:,n)),[8 8 1]);
end
base = reshape(ifisherz(median(cormat_bz,3)),[8 8 1]);
base(logical(eye(8,8))) = 1;
test = reshape(ifisherz(median(cormat_tz,3)),[8 8 1]);
test(logical(eye(8,8))) = 1;
figure
subplot(1,2,1);
imagesc(base,[-.1 .1]);
axis square;colorbar;
set(gca,'fontsize',15)
axis off
subplot(1,2,2);
imagesc(test,[-.1 .1]);
axis square;colorbar
set(gca,'fontsize',15)
axis off

%% try the same as MDS
tsub = length(unique(fix.subject));
dismat_t = zeros(8,8,tsub);
subc            = 0;
for subject = unique(fix.subject);
    subc = subc + 1
    %creaete the query cell
    v = [];
    c = 0;
    for ph = [4]
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject' subject};
        end
    end
    % plot and save fixation maps
    fix.getmaps(v{:});
    fix.maps   = fix.maps(:,:,1:8) - repmat(mean(fix.maps(:,:,1:8),3),[1 1 8]);%baseline, take out the average
    fixmap = fix.vectorize_maps;
    
    dismat_t(:,:,subc) = squareform(pdist(fixmap','correlation'));
    
end
%%
sc = 0;
for sub = 1:27%[1:22 25:27]
    sc=sc+1
    sub
    mds(:,:,sc)   = mdscale(dismat_t(:,:,sub),2,'criterion','metricstress');
end
% mds(:,:) = mdscale(mean(dismat_t,3),2,'criterion','metricstress');
%%
% plot plot the MDS results
colors = GetFearGenColors;
for n=1:length(mds)
    for f = 1:8
        set(gcf,'DefaultAxesColorOrder',colors);
        plot(mds(f,1,n),mds(f,2,n),'wo','MarkerFaceColor',colors(f,:),'MarkerSize',20)
        hold on
    end
end
hold on
axis square
%%


%% compute single subject maps and similarity matrices, and finally average them across subjects.
tsub = length(unique(fix.subject));
cormat_t = nan(8,8,tsub,6);
cormat_b = nan(8,8,tsub,6);
for k = 25%Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
    fix.kernel_fwhm = k;
    for nfix = [1 2 3 4 5 6];%there are about 200 7th fixations (not enough)
        subc   = 0;
        for subject = unique(fix.subject);
            subc = subc + 1;
            %creaete the query cell
            v = [];
            c = 0;
            for ph = [2 4]
                for cond = -135:45:180
                    c    = c+1;
                    v{c} = {'phase', ph, 'deltacsp' cond 'fix' nfix 'subject' subject};
                end
            end
%             figure(10);clf;
            % plot and save fixation maps
            fix.getmaps(v{:});
            maps                 = fix.maps;
            fix.maps             = maps(:,:,1:8)   - repmat(mean(maps(:,:,1:8),3),[1 1 8]);%correct for baseline
            cormat_b(:,:,subc,nfix)   = fix.corr;
            fix.maps   = maps(:,:,9:end) - repmat(mean(maps(:,:,9:end),3),[1 1 8]);%take out the average            
            cormat_t(:,:,subc,nfix)   = fix.corr;            
        end
    end
end
%%
for nf = 1:4;
    subplot(2,4,nf)
    imagesc([ nanmean(cormat_b(:,:,:,nf),3) ],[-.2 .1]);
    title(sprintf('Fixation No: %d',nf))
    if nf == 1
        ylabel('Baseline')
    end
    axis square
    axis off
    colorbar
    subplot(2,4,nf+4)
    imagesc([ nanmean(cormat_t(:,:,:,nf),3)],[-.2 .1]);
    axis square
    colorbar
    axis off
    if nf == 1
        ylabel('Test')
    end
end
SaveFigure('/Users/onat/Desktop/SimilarityOFConditionsFixByFix','-r150');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare Single subject across task similarity
sub_c = 0;
tsub  = length(unique(fix.subject));
cormat = nan(5,5,tsub);
for subs = unique(fix.subject);    
    sub_c = sub_c+1;
    v    = [];
    c    = 0;
    for ph = [1 2 3 4 5]
        c    = c+1;
        v{c} = {'phase', ph , 'deltacsp' fix.realcond 'subject' subs};
    end
    fix.getmaps(v{:});    
    fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 5]);
    cormat(:,:,sub_c) = fix.corr;
end
mat = reshape(ifisherz(mean(reshape(fisherz(cormat),[5 5 27]),3)),[5 5]);
mat = CancelDiagonals(mat,1);
imagesc(mat)
colorbar
%% COMPARE FIXATION PATTERNS ACROSS THE 5 PHASES (MEGASUBJECT)
sub_c  = 0;
v      = [];
c      = 0;
for ph = [1 2 3 4 5]
    c    = c+1;
    v{c} = {'phase', ph , 'deltacsp' [0 180 18000] };
end
fix.getmaps(v{:});
cormat_group= fix.corr;
%%
v    = [];
c    = 0;
for ph = [1 2 3 4 5]
    for subs = unique(fix.subject);
        c    = c+1;
        v{c} = {'phase', ph , 'deltacsp' fix.realcond 'subject' subs};
    end
end
fix.getmaps(v{:});
maps = fix.maps;

%     fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 5]);
cormat= fix.corr;
imagesc(cormat);
%%
for nx = 1:5
    ix = [1:27]+27*(nx-1);
    for ny = 1:5
        iy = [1:27]+27*(ny-1);
        bla(ny,nx) = mean(diag(cormat(iy,ix)))        
    end
end


%% Plot simple fixation counts for CS+ and CS? and for different phases.
p               = Project;
mask            = p.getMask('ET');
subjects        = find(mask(:,1));
subjects        = intersect(subjects,Project.subjects_1500);
M = [];S = [];
for np = 1:5
    fix             = Fixmat(subjects,np);
    [c i ]          = fix.histogram;
    close all;
    M(1,np)           = mean(squeeze(c(:,ismember(i,[0]),:)));
    M(2,np)           = mean(squeeze(c(:,ismember(i,[180 18000]),:)));
    S(1,np)           =  std(squeeze(c(:,ismember(i,[0]),:))./sqrt(28));    
    S(2,np)           =  std(squeeze(c(:,ismember(i,[180 18000]),:))./sqrt(28));
end
%
errorbar(M(1,:),S(1,:),'bo-','linewidth',3);
hold on
errorbar(M(2,:),S(2,:),'ro-','linewidth',3);
hold off;
box off
set(gca,'xtick',1:5,'xticklabel',{'Disc1' 'B' 'C' 'T' 'Disc2'});ylabel('Fix Count');
grid on;
xlim([.5 5.5]);
legend({'CS+' 'CS?'});legend boxoff
SaveFigure(sprintf('/Users/onat/Desktop/FixCount1500.png'));
%% Plot simple fixation counts for 1st vs 2nd face of a trial in the discrimination task
clear all
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = find(mask);
subjects        = intersect(subjects,Project.subjects_1500);
M = [];S = [];
%cut out the odd (2nd facec of trial) Trials
first  = [1:2:400]';
second = [2:2:400]';
for np=[1 5]
    cc = 0;
    for fs = {first second};
        fix             = Fixmat(subjects,np);
        cc=cc+1;
        v = {'trialid' fs{1}};
        fix.UpdateSelection(v{:})
        fix.ApplySelection
        [c i]          = fix.histogram;
        close all;
        M(cc,np)           = mean(nanmean(c,2));%mean across cond and subjects
        S(cc,np)           = mean(nanstd(c,0,2))/length(unique(fix.subject));
    end
end
%
figure
errorbar([1 2],M(1,[1 5]),S(1,[1 5]),'bo-','linewidth',3);
hold on;
errorbar([1.2 2.2],M(2,[1 5]),S(2,[1 5]),'ro-','linewidth',3)
set(gca,'xtick',[1.1,2.1],'xticklabel',{'Disc1' 'Disc2'});ylabel('Fix Count');
legend('1st Face','2nd Face')
%% plot corresponding fixmaps for first vs second face in discrimination
fix = Fixmat(g.ids,[1 5]);
c=0;
v=[];
for fs = {first second};
        c=c+1;
        v{c} = {'trialid' fs{1}};
end
fix.getmaps(v{:});
fix.plot
%% Plot the discrimination thresholds as a scatterhist
p        = Project;
mask     = find(p.getMask('ET_feargen').*prod(p.getMask('PMF'),2));
subjects = intersect(Project.subjects_1500,mask);
g1500    = Group(subjects);
subjects = intersect(Project.subjects_600,mask);
g600    = Group(subjects);
%
clf;figure(1)
dummy = g1500;
scatterhist(dummy.pmf.csp_before_alpha,dummy.pmf.csp_after_alpha,'nbins',15);
axis square;
xlim([0 100]);ylim([0 100])
DrawIdentityLine(gca);
xlabel('before');ylabel('after');title('1500ms');
SaveFigure(sprintf('/Users/onat/Desktop/1500.png'),'-r250');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare fixation maps for GOOD vs. BAD discriminators during Discr-1
p        = Project;
mask     = find(p.getMask('ET_Discr').*prod(p.getMask('PMF'),2));
subjects = intersect(Project.subjects_1500,mask);%subjects
fix      = Fixmat(subjects,1);
g        = Group(subjects);%get the data for these subjects
param    = g.pmf.csn_before_alpha + g.pmf.csp_before_alpha;%the parameter of interest
i        = param < median(param);%median
%
c = 0;v = [];
for subs = {subjects(i) subjects(~i)}%good and then bad
    c       = c+1;
    v{c}    = {'phase', 1 , 'deltacsp' [0 18000] 'subject' subs{1} 'fix' 3};
end
%
fix.kernel_fwhm = 35;
fix.getmaps(v{:});
% fix.plot('log');
%
fix.maps = -diff(fix.maps,1,3);%bad - good.
fix.plot('linear');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fixation Counts and behavior
[c index ]      = fix.histogram;
%simple correlations
[mat labels] = g.parameterMat
%SI
[r p] = corr(nanmean(c,2),mat(:,end))
plot(nanmean(c,2),mat(:,end),'bo','MarkerFaceColor','b','LineWidth',2);ls=lsline;set(ls,'LineWidth',2)
%threshold
[r p] = corr(nanmean(c,2),mean(mat(:,1:4),2))
plot(nanmean(c,2),mean(mat(:,1:4),2),'bo','MarkerFaceColor','b','LineWidth',2);ls=lsline;set(ls,'LineWidth',2)
%%
close all;
subs            = {subjects(i) subjects(~i)}; %good and then bad
m               = [mean(c(ismember(index.sub,subs{1}),index.cond == 0)) mean(c(ismember(index.sub,subs{2}),index.cond == 0))];
s               = [std(c(ismember(index.sub,subs{1}),index.cond == 0))./sqrt(length(subs{1})) std(c(ismember(index.sub,subs{2}),index.cond == 0))./sqrt(length(subs{2}))];
errorbar(m,s);
set(gca,'xtick',[1 2],'xticklabel',{'Good' 'Bad'});ylabel('fixation count');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's do the same as before, but first colloect individual subject fixation maps and then median split the fixation maps
p        = Project;
mask     = find(p.getMask('ET_feargen').*prod(p.getMask('PMF'),2));
subjects = intersect(Project.subjects_1500,mask);%subjects
fix      = Fixmat(subjects,1);
g        = Group(subjects);%get the data for these subjects
param    = (g.pmf.csp_before_alpha + g.pmf.csn_before_alpha)./2;%the parameter of interest
i        = param < median(param);%median
%
c = 0;v = [];
for subs = subjects(:)'%good and then bad
    c       = c+1;
    v{c}    = {'phase', 1 , 'deltacsp' [0 18000] 'subject' subs};
end
fix.kernel_fwhm = 35;
fix.getmaps(v{:});
fix.maps        = cat(3,mean(fix.maps(:,:,i),3),mean(fix.maps(:,:,~i),3));

fix.plot('log');
%
fix.maps = -diff(fix.maps,1,3);%bad - good.
fix.plot('linear');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute correlation between fixation counts and discrimination threshold.
p        = Project;
mask     = find(p.getMask('ET_feargen').*prod(p.getMask('PMF'),2));
subjects = intersect([Project.subjects_1500],mask);%subjects
fix      = Fixmat(subjects,1);
%
g        = Group(subjects);%get the data for these subjects
param    = -(g.pmf.csp_before_alpha(:) + g.pmf.csn_before_alpha(:));%the parameter of interest
%
c = 0;v = [];
for subs = subjects(:)'%good and then bad
    c       = c+1;
    v{c}    = {'phase', 1 , 'deltacsp' [0 18000] 'subject' subs};
end
fix.unitize     = 1;
fix.kernel_fwhm = 35;
fix.getmaps(v{:});
M = fix.vectorize_maps';
for i = 1:size(M,2)
    if sum(abs(M(:,i))) ~= 0
        r(i) = corr2(M(:,i),param);
    else
        r(i) = 0;
    end
end
fix.maps = reshape(r,size(fix.maps,1),size(fix.maps,2));
fix.plot
%% Compute correlation between fixation counts and generalization width.
p        = Project;
mask     = find(p.getMask('ET_feargen').*prod(p.getMask('RATE'),2));
subjects = intersect([Project.subjects_1500],mask);%subjects
subjects = setdiff(subjects,9)
fix      = Fixmat(subjects,4);
%
g        = Group(subjects);%get the data for these subjects
g.getSI(3);

%%
[M i]       = g.parameterMat;
param       = M(:,end);%sharpening index
%
c = 0;v = [];
for subs = subjects(:)'%good and then bad
    c       = c+1;
    v{c}    = {'phase', 4  'subject' subs };
end
fix.unitize     = 1;
fix.kernel_fwhm = 35;
fix.getmaps(v{:});
M = fix.vectorize_maps';
for i = 1:size(M,2)
    if sum(abs(M(:,i))) ~= 0
        r(i) = corr2(M(:,i),param);
    else
        r(i) = 0;
    end
end
fix.maps = reshape(r,size(fix.maps,1),size(fix.maps,2));
fix.plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get fixation counts, as well as mean duration for subject (like Kanan et al, 2015)
load('C:/Users/onat/Google Drive/EthnoMaster/data/midlevel/singletrialfixmaps/labels1500.mat')
fixnum_m = NaN(length(unique(labels.sub)),5);
fixnum_sd = NaN(length(unique(labels.sub)),5);
fixdur_m = NaN(length(unique(labels.sub)),5);
fixdur_sd = NaN(length(unique(labels.sub)),5);
cc=0;
for ns=unique(labels.sub)
    cc=cc+1;
    for phase = 1:5
        fix            = Fixmat(ns,phase);
        [c i]          = fix.histogram;
        close all;
        fixnum_m(cc,phase)        = mean(c);%mean fix per phase
        fixnum_sd(cc,phase)       = std(c);%std fix per phase
        fixdur_m(cc,phase)        = mean(fix.stop-fix.start);%mean fix per phase
        fixdur_sd(cc,phase)       = std(double(fix.stop-fix.start));%std fix per phase
    end
end
%
errorbar(M(1,:),S(1,:),'bo-','linewidth',3);
hold on
errorbar(M(2,:),S(2,:),'ro-','linewidth',3);
hold off;
box off
set(gca,'xtick',1:5,'xticklabel',{'Disc1' 'B' 'C' 'T' 'Disc2'});ylabel('Fix Count');
grid on;
xlim([.5 5.5]);
legend({'CS+' 'CS?'});legend boxoff
SaveFigure(sprintf('/Users/onat/Desktop/FixCount1500.png'));
%% PLOT FIXMAPS for different phases broken down to CS+ and CS?
p        = Project;
mask     = find(p.getMask('ET_feargen'));
subjects = intersect(Project.subjects_1500,mask);%subjects
fix      = Fixmat(subjects,[2 3 4]);
fix.unitize = 1;
%%
v     = [];
c     = 0;
v     = {'phase', 2 , 'deltacsp' [0 180]}
fix.kernel_fwhm = 35;
fix.getmaps(v)
blank = fix.maps;
v     = {{'phase', 4 , 'deltacsp' 0} {'phase', 4 , 'deltacsp' 180}};
fix.getmaps(v{:})
fix.maps = fix.maps - repmat(blank,[1 1 2]);;
% fix.maps = -diff(fix.maps,1,3)
fix.plot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correlation of subject fixation patterns between different phases

subc=0;
for subs = subjects(:)'
    subc=subc+1;
    c = 0;v = [];
    for ph = 1:5
        c       = c+1;
        v{c}    = {'phase', ph , 'subject' subs};
    end
    fix.getmaps(v{:});
    corrmat(:,:,subc) = fix.corr;
end

%% Which parts of the fixation map correlate with good alpha
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
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmaps = fix.vectorize_maps;

%normal correlation
for n = 1:length(subjmaps)
    [r(n),pval(n)] = corr(subjmaps(n,:)',mean(mat(:,[1 3]),2));
end
imagesc(reshape(r,[50 50]),[-1.5 1.5])
for n = find(pval<0.05);[y x] = ind2sub([50 50],n);text(x,y,'x','FontSize',20);end
fix.maps = reshape(r,[50 50]);
fix.plot

%bootstrap
for n = 1:length(subjmaps)
    data = [subjmaps(n,:)',mean(mat(:,[1 3]),2)];
    if ~any(sum(data)==0)
        n
        [dummy dummy2] = CorrelationMatrix(data, 100,0.05,'bootci');
        cmat(n) = dummy(1,2);
        pmat(n) = dummy2(1,2);
    else
        cmat(n) = NaN;
        pmat(n) = 0;
    end
end
%mask: multiply with average fixmap
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = Scale(mean(fix.maps,3));
fix.maps = reshape(r,[50 50]).*fix.maps;
fix.plot


% search pixel that is significant
[y,x] = ind2sub([50 50],intersect(find(r<0),find(pval<0.01)));
imagesc(reshape(r,[50 50]).*fix.maps)
text(x,y,'x','FontSize',20,'Color','w')






%binary mask
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = mean(fix.maps,3);
bla = fix.maps;
fix.maps = bla > prctile(bla(:),5);
fix.maps = reshape(cmat,[50 50]).*fix.maps;
fix.plot

%% Which parts of the fixation map correlate with an increased generalization.
clear all
p           = Project;
mask        = find(p.getMask('ET_feargen').*p.getMask('RATE'));
subjects    = intersect([Project.subjects_1500],mask);%subjects
g           = Group(subjects);%get the data for these subjects
[mat tags]  = g.parameterMat;
mises = g.loadmises; mat = [mat mises];
clear g;
fix         = Fixmat(subjects,4);
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
subjmaps = fix.vectorize_maps;

%normal correlation
param       = mat(:,13);%kappa test
for n  = 1:25; param(n,1) = vM2FWHM(mat(n,13));end

for n = 1:length(subjmaps)
    [r(n),pval(n)] = corr(subjmaps(n,:)',param);
end
imagesc(reshape(r,[50 50]),[-1 1])
fix.maps = reshape(r,[50 50]);fix.plot
for n = find(pval<0.05);[y x] = ind2sub([50 50],n);text(x,y,'x','FontSize',20);end
fix.maps = reshape(r,[50 50]);
fix.plot

%bootci
for n = 1:length(subjmaps)
    data = [subjmaps(n,:)',mat(:,14)];
    if ~any(sum(data)==0)
        n
        [dummy dummy2] = CorrelationMatrix(data, 100,0.05,'bootstrp');
        dummy(isnan(dummy)) = 0;
        cmat(n) = mean(dummy);
%         pmat(n) = dummy2;
    else
        cmat(n) = NaN;
        pmat(n) = 0;
    end
end

%mask: multiply with average fixmap
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = Scale(mean(fix.maps,3));
fix.maps = reshape(r,[50 50]).*fix.maps;
fix.plot
%binary mask
fix.getsubmaps;
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
fix.maps = mean(fix.maps,3);
bla = fix.maps;
fix.maps = bla > prctile(bla(:),5);
fix.maps = reshape(cmat,[50 50]).*fix.maps;
fix.plot

%% median split good vs bad generalizers
i           = param > median(param);%median

c = 0;
v = [];
for subs = {subjects(i) subjects(~i)}%good and then bad
    c       = c+1;
    v{c}    = {'phase', 4 , 'deltacsp' [0] 'subject' subs{1}};
end
%
fix.kernel_fwhm = 35;
fix.getmaps(v{:});
fix.plot;
% fix.maps = -diff(fix.maps,1,3);fix.plot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Gaussian Tuned locations
p           = Project;
mask        = find(p.getMask('ET_feargen'));
subjects    = intersect(Project.subjects_1500,mask);%subjects
fix         = Fixmat(subjects,[2 4]);
%%
c           = 0;
v           = [];
for nph = 4
    for cond = -135:45:180;%good and then bad
        c    = c+1;
        v{c} = {'phase', nph , 'deltacsp' cond};
    end
end
fix.maptype = 'bin';
fix.binsize = 15;
fix.unitize = 0;
fix.getmaps(v{:});
M           = fix.vectorize_maps;

x = -135:45:180;
for npix = 1:size(M,1)
    y = M(npix,:);
    if sum(y) > .005
        try
            dummy   = FitGauss(x,y,3);
            Amplitude(npix) = dummy.Est(1);
            pval(npix)      = dummy.pval;
        catch
            Amplitude(npix) = NaN;
            pval(npix)      = NaN;
        end
    else
        Amplitude(npix) = NaN;
        pval(npix)      = NaN;
    end
end
fix.maps = reshape(pval,33,33);
fix.plot;





    









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%good subject? bad subject? %regarding Discrimination Task Fixations
p               = Project;
mask            = p.getMask('ET');
subjects        = find(mask(:,1));
subjects        = intersect(subjects,Project.subjects_1500);
fix             = Fixmat(subjects,[1 5]);

subjects = intersect(Project.subjects_1500,find((ETmask(:,1)==1)));
g = Group(subjects);
fix             = Fixmat(subjects,[2 3 4]);
%GOOOD vs BAD
subject_alpha=[];
for i = 1:length(g.ids)
    subject_alpha = [subject_alpha; mean(g.pmf.params1(:,1,i),1)];
end
med = median(subject_alpha);
bad = g.ids(find(subject_alpha>=med));
good  = g.ids(find(subject_alpha<med));

%%
tsub = length(unique(fix.subject));
k = 20;
    fix.kernel_fwhm = k;
    covmat          = nan(16,16,tsub);
    cormat          = nan(16,16,tsub);
    subc            = 0;
    %create the query cell
    v = [];
    c = 0;
    for subjects = {good bad}
        for ph = [2 4]
            for cond = -135:45:180
                c    = c+1;
                v{c} = {'phase', ph, 'deltacsp' cond 'subject' subjects{1}};
            end
        end
    end
    % plot and save fixation maps
    fix.getmaps(v{:});
    %%
    %for one group only
    %fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
    %fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
   
    %%
    %%Good - Bad:
    %first Baseline correction, but within Group
    fix.maps(:,:,1:16)             = fix.maps(:,:,1:16)    - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for group_baseline
    fix.maps(:,:,17:32)            = fix.maps(:,:,17:32)   - repmat(mean(fix.maps(:,:,17:24),3),[1 1 16]);%correct for group_baseline
    %Mean Testphase correction (within group)
    fix.maps(:,:,9:16) = fix.maps(:,:,9:16) - repmat(mean(fix.maps(:,:,9:16),3),[1 1 8]);
    fix.maps(:,:,25:32) = fix.maps(:,:,25:32) - repmat(mean(fix.maps(:,:,25:32),3),[1 1 8]);
    %take the diff (Good - Bad)
    fix.maps(:,:,1:16) = fix.maps(:,:,1:16) - fix.maps(:,:,17:32);
    fix.maps(:,:,17:32) = []; %need 2x8cond for fix.plot
    fix.plot
    t=supertitle('Good - Bad Base&Mean corrected (600ms)',1);
    set(t,'FontSize',14)
  %% 
%SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/AllFixations/baselinecorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
saveas(gcf,'C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/good_vs_bad/within_group_corrected/1500_Good-Bad_k20_baseANDmeancorr.png');


%%
%%%improvement
%what is subject's improvement? For CSP? or mean(CSP/CSN)?
%%CSP
load(sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep));
subjects = intersect(Project.subjects_600,find((ETmask(:,1)==1)));
g = Group(subjects);
fix             = Fixmat(subjects,[2 3 4]);

alpha_csp_impr = squeeze(g.pmf.params1(1,1,:) - g.pmf.params1(3,1,:));
csp_improvers = g.ids(find(alpha_csp_impr>median(alpha_csp_impr)));
csp_nonimpr   = g.ids(find(alpha_csp_impr<=median(alpha_csp_impr)));

k = 20;
fix.kernel_fwhm = k;
%create the query cell
v = [];
c = 0;
for subjects = {csp_improvers csp_nonimpr}
    for ph = [2 4]
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject' subjects{1}};
        end
    end
end
% plot and save fixation maps
fix.getmaps(v{:});
uncorr=fix.maps;%save it to have it for inbetween saving steps
%%
%first Baseline correction, but within Group
fix.maps(:,:,1:16)             = fix.maps(:,:,1:16)    - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for group_baseline
fix.maps(:,:,17:32)            = fix.maps(:,:,17:32)   - repmat(mean(fix.maps(:,:,17:24),3),[1 1 16]);%correct for group_baseline
%Mean Testphase correction (within group)
fix.maps(:,:,9:16) = fix.maps(:,:,9:16) - repmat(mean(fix.maps(:,:,9:16),3),[1 1 8]);
fix.maps(:,:,25:32) = fix.maps(:,:,25:32) - repmat(mean(fix.maps(:,:,25:32),3),[1 1 8]);
%take the diff (Improved - Non-Improved)
fix.maps(:,:,1:16) = fix.maps(:,:,1:16) - fix.maps(:,:,17:32);
fix.maps(:,:,17:32) = []; %need 2x8cond for fix.plot
fix.plot
t=supertitle('Improved - Non-Improved Base&Mean corrected (600ms)',1);
set(t,'FontSize',14)

saveas(gcf,'C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/improvement_CSP/within_group_corrected/improved_k20_uncorrected.png');

%%
%Improvement concerning both CSP and CSN (their mean)
load(sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep));
subjects = intersect(Project.subjects_1500,find((ETmask(:,1)==1)));
g = Group(subjects);
fix             = Fixmat(subjects,[2 3 4]);

alpha_both_impr = squeeze(mean(g.pmf.params1(1:2,1,:)) - mean(g.pmf.params1(3:4,1,:)));
csp_improvers = g.ids(find(alpha_both_impr>median(alpha_both_impr)));
csp_nonimpr   = g.ids(find(alpha_both_impr<=median(alpha_both_impr)));
%%

k = 20;
fix.kernel_fwhm = k;
%create the query cell
v = [];
c = 0;
for subjects = {csp_improvers csp_nonimpr}
    for ph = [2 4]
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject' subjects{1}};
        end
    end
end
% plot and save fixation maps
fix.getmaps(v{:});
uncorr=fix.maps;%save it to have it for inbetween saving steps

%first Baseline correction, but within Group
fix.maps(:,:,1:16)             = fix.maps(:,:,1:16)    - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for group_baseline
fix.maps(:,:,17:32)            = fix.maps(:,:,17:32)   - repmat(mean(fix.maps(:,:,17:24),3),[1 1 16]);%correct for group_baseline
%Mean Testphase correction (within group)
fix.maps(:,:,9:16) = fix.maps(:,:,9:16) - repmat(mean(fix.maps(:,:,9:16),3),[1 1 8]);
fix.maps(:,:,25:32) = fix.maps(:,:,25:32) - repmat(mean(fix.maps(:,:,25:32),3),[1 1 8]);
%take the diff (Improved - Non-Improved)
fix.maps(:,:,1:16) = fix.maps(:,:,1:16) - fix.maps(:,:,17:32);
fix.maps(:,:,17:32) = []; %need 2x8cond for fix.plot
fix.plot
t=supertitle('Improved - Non-Improved Base&Mean corrected (1500ms)',1);
set(t,'FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% show fixations for different phases
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_600);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
fix             = Fixmat(subjects,[1 2 3 4 5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Discrimination Task
%% get the data for discrimination task
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_600);
fix             = Fixmat(subjects,[1 5]);
%% get the data for discrimination task
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
fix             = Fixmat(subjects,[1 5]);
%% get the query to compare CS+ vs .CS?
fix.kernel_fwhm = 20;
%create the query cell
v = [];
c = 0;
for ph = [1 5]%before and after
    for cond = [0,18000]%CS+ and CS?
        c    = c+1;
        v{c} = {'phase', ph, 'deltacsp' cond };
    end
end
% plot and save fixation maps
fix.getmaps(v{:});
dummy = fix.maps;
%fix.plot;
m = mean(dummy(:,1:2),2);
%r= corrcoef(dummy-repmat(m,[1 4])
%% Median split according to discrimination threshold
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),subjects);
fix             = Fixmat(subjects,[1 5]);
g               = Group(subjects');
subject_alpha=[];
for i = 1:length(g.ids)
    subject_alpha = [subject_alpha; mean(g.pmf.params1(1:2,1,i),1) mean(g.pmf.params1(3:4,1,i),1)];
end

med_before   = median(subject_alpha(:,1));
bad_before   = g.ids(find(subject_alpha(:,1)>=med_before));
good_before  = g.ids(find(subject_alpha(:,1)<med_before));
med_after    = median(subject_alpha(:,2));
bad_after    = g.ids(find(subject_alpha(:,2)>=med_after));
good_after   = g.ids(find(subject_alpha(:,2)<med_after));


mean_goodbefore=mean(subject_alpha(ismember(g.ids,good_before),1))
mean_badbefore=mean(subject_alpha(ismember(g.ids,bad_before),1))
mean_goodafter=mean(subject_alpha(ismember(g.ids,good_before),2))
mean_badafter=mean(subject_alpha(ismember(g.ids,bad_before),2))
%% one map per subject and phase - 135 instances
%create the query cell
cond = [-16875:1125:18000,-135:45:180];
v = [];
c = 0;
for ph = 1:5
    for sub = unique(fix.subject)
        c    = c+1;
        v{c} = {'phase' ph 'subject' sub  'deltacsp' cond };
    end
end
% plot and save fixation maps
fix.getmaps(v{:});
fix.maps   = imresize(fix.maps,.1,'method','bilinear');
datamatrix = fix.vectorize_maps;
%% where do good vs bad discriminators look? (independent from condition)
%GOOD Discrimination BEFORE
%create the query cell
v = [];
c = 0;
for subjects = {good_before bad_before}
    c    = c+1;
    v{c} = {'phase',1, 'deltacsp' [0,18000] 'subject' subjects{1}};
end
% plot and save fixation maps
fix.getmaps(v{:});
%% where do good vs bad discriminators look? (independent from condition)
%GOOD Discrimination AFTER
%create the query cell
v = [];
c = 0;
for subjects = {good_after bad_after}
    c    = c+1;
    v{c} = {'phase',5, 'deltacsp' [0,18000] 'subject' subjects{1}};
end
% plot and save fixation maps
fix.getmaps(v{:});

%% 
%CS+ vs CS_ after
%create the query cell
v = [];
c = 0;
for subjects = {good_after bad_after}
    for phase = [1 5]
        for cond = [0 18000]
            c    = c+1;
            v{c} = {'phase',phase, 'deltacsp' cond 'subject' subjects{1}};
        end
    end
end
% plot and save fixation maps
fix.getmaps(v{:});
uncorr=fix.maps;
bcorr = cat(3,[uncorr(:,:,3:4)-repmat(mean(uncorr(:,:,1:2),3),[1 1 2])],[uncorr(:,:,7:8)-repmat(mean(uncorr(:,:,5:6),3),[1 1 2])]);%correct for individual_baseline (which is own CSPandCSN mean in phase1);
%% high improvers vs. low improvers.
%mean improvement in both CSP and CSN
alpha_impr = squeeze(mean(g.pmf.params1(1:2,1,:)) - mean(g.pmf.params1(3:4,1,:)));
improvers = g.ids(find(alpha_impr>median(alpha_impr)));
nonimpr   = g.ids(find(alpha_impr<=median(alpha_impr)));

mean_improvers=mean(alpha_impr(ismember(g.ids,improvers)));
mean_nonimpr=mean(alpha_impr(ismember(g.ids,nonimpr)));
%% 
%create the query cell
v = [];
c = 0;
for subjects = {improvers nonimpr}
    for phase = [1 5]
        c    = c+1;
        v{c} = {'phase', phase, 'deltacsp' [0, 18000] 'subject' subjects{1}};
    end
end

% plot and save fixation maps
fix.getmaps(v{:});
fix.maps = fix.maps(:,:,[2,4]) - repmat(mean(fix.maps(:,:,[1,3]),3),[1 1 2]);

%%
%improvement in CSP only
alpha_csp_impr = squeeze(g.pmf.params1(1,1,:) - g.pmf.params1(3,1,:));
csp_improvers = g.ids(find(alpha_csp_impr>median(alpha_csp_impr)));
csp_nonimpr   = g.ids(find(alpha_csp_impr<=median(alpha_csp_impr)));
%% high guessing rate vs low guessing rate
for i = 1:length(g.ids)
    guess(i) = mean(g.subject{i}.pmf.params1(:,3));
end
med = median(guess);
high_guess = g.ids(find(guess>med))';
low_guess  = g.ids(find(guess<=med))';
%create the query cell
v = [];
c = 0;
for subjects = {high_guess low_guess}
        c    = c+1;
        v{c} = {'subject' subjects{1}};
        v{c}
end
fix.getmaps(v{:});
fix.plot

%%
%show 5 phases next to each other
%% get the data for discrimination task
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
fixd            = Fixmat(subjects,[1 5]);
fixfg           = Fixmat(subjects,[2 3 4]);

fixd.getmaps({'phase',1,'deltacsp',[0,18000]},{'phase',5,'deltacsp',[0,18000]});
fixfg.getmaps({'phase',2,'deltacsp',[0,180]},{'phase',3,'deltacsp',[0,180]},{'phase',4,'deltacsp',[0,180]});
discr   = fixd.maps;
feargen = fixfg.maps;
raw     = cat(3,feargen,discr);
blank   = mean(raw,3);
corr    = raw - repmat(blank, [ 1 1 5]);
%%
clear all;
p         = Project;
rate_mask = find(p.getMask('RATE'));
pmf_mask  = find(sum(p.getMask('PMF'),2) == 4);
subs      = intersect(intersect(rate_mask,pmf_mask),Project.subjects_1500);
g         = Group(subs);
g.getSI(3);
[mat labels] = g.parameterMat;
%%
pmf          = [mean(mean(mat(:,[1 2]))) mean(mean(mat(:,[5 6])))];
ratings      = [mean(mean(mat(:,[end-2:end-1])))];
x = linspace(0,180,10000);
plot(x,1-PAL_Weibull(pmf,x),'r')
hold on;plot(x,make_gaussian1d(x,1,ratings,0),'g--');

%%
[mat labels] =g.parameterMat;
[cmat mask] = CorrelationMatrix(mat,10000,0.01,'bootci');
set(gca,'yticklabel',labels,'ytick',1:length(labels));
%%


%%
p=Project;
pmf_mask  = find(sum(p.getMask('PMF'),2) == 4);
subs      = intersect(pmf_mask,Project.subjects_1500);
subs      = intersect(subs,find(p.getMask('ET_discr')));
g         = Group(subs);
g.getSI(3);
[mat labels] =g.parameterMat;

fix = Fixmat(subs,[1 5]);


%% fixation maps for CSP improvers
med = median(mat(:,9));
good_imp = g.ids(find(mat(:, 9) > med));
bad_imp  = g.ids(find(mat(:,9) <= med));

% create the query cell
v = [];
c = 0;
for subjects = {good_imp bad_imp} 
        c    = c+1;
        v{c} = {'phase', [5], 'deltacsp' [0] 'subject' subjects{1}};
end
% plot and save fixation maps
fix.getmaps(v{:});
fix.plot
fix.maps = -diff(fix.maps,1,3);
fix.plot
% create the query cell
v = [];
c = 0;
for sub = subs' 
    for phase = [1]
        c    = c+1;
        v{c} = {'phase', phase, 'deltacsp' [0 18000] 'subject' sub  };
    end
end

fix.getmaps(v{:});
% fix.maps = fix.maps(:,:,26:50);
m  = fix.vectorize_maps;

for i = 1:size(m,1)
   r(i)=corr2(m(i,:)',median(mat(:,[1 3]),2));
end

fix.maps = reshape(r,[500 500]);
fix.plot;

%% fixations per pixel, correlation with delta alpha
med = median(mean(mat(:,[1 3]),2));
good_imp = g.ids(find(mean(mat(:,[1 3]),2) < med));
bad_imp  = g.ids(find(mean(mat(:,[1 3]),2) >= med));

v = [];
c = 0;
for sub = {good_imp bad_imp}
    for phase = [1]
        c    = c+1;
        v{c} = {'phase', phase, 'deltacsp' [0 18000] 'subject' sub{1}};
    end
end

fix.getmaps(v{:});
fix.maps=-diff(fix.maps,1,3);
fix.plot
%% subject fixation maps
v = [];
c = 0;
for sub = g.ids
        c    = c+1;
        v{c} = {'subject' sub};
end
fix.getmaps(v{:});
%% phase fixation maps
v = [];
c = 0;
for ph = 1:5
        c    = c+1;
        v{c} = {'phase' ph 'subject' 6};
end
fix.getmaps(v{:});
%% Can generalization pattern predict SI?
% w = model.SVs'*model.sv_coef;
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\singletrialfixmaps\1500\alpha_before_N25\discr_hp.mat')
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\singletrialfixmaps\1500\alpha_before_N25\labels.mat')
p = Project;
mask = p.getMask('RATE');
subjects = intersect(find(mask),unique(labels.sub));

g          = Group(subjects);
g.getSI(3);
[mat tags] = g.parameterMat;
alpha_bef  = mean(mat(:,[1 3]),2);
fix        = Fixmat(g.ids,4);
fix.kernel_fwhm = 37;
%get subjmaps
v = [];
c = 0;
for sub = unique(fix.subject)
        c    = c+1;
        v{c} = {'subject' sub};
end
fix.getmaps(v{:});
fix.maps   = imresize(fix.maps,.1,'method','bilinear');
subjmap    = fix.vectorize_maps;
hpload     = discr_hp(:)'*subjmap;

[r,p] = corr(hpload',g.SI)
figure;
plot(hpload,g.sigma_test,'b.','MarkerSize',30)
l = lsline; set(l,'LineWidth',2);
xlabel('subjmap phase 4 x discr hyperplane')
ylabel('SI')

%% confusion of 4-5 and 4-1
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\svm_analysis\versiona733ddc\phases_insubject_rand0\result.mat')
a = squeeze(mean(result,3));%bootstraps
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,5]);%scale by rowsums
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\singletrialfixmaps\N27\labels.mat')
sub_c = unique(labels.sub);
p = Project;
mask = p.getMask('RATE');
g = Group(intersect(find(mask),sub_c));
g.getSI(3);
[mat tags] = g.parameterMat;
conf45 = squeeze(scaled(4,5,ismember(sub_c,g.ids)));
conf41 = squeeze(scaled(4,1,ismember(sub_c,g.ids)));
%correlate confusion with sigma_test and SI
[r,p]=corr(conf45,mat(:,13))
[r,p]=corr(conf45,mat(:,14))
[r,p]=corr(conf41,mat(:,13))
[r,p]=corr(conf41,mat(:,14))
%% 
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('ET_discr');
subjects = intersect(subjects,find(mask));
fix = Fixmat(subjects,1:5);

%% correlation matrix for subjects (even / odd trials)
N = length(unique(fix.subject));
v = [];
c = 0;
for sub=unique(fix.subject)
    c    = c+1;
    v{c} = {'subject' sub 'trialid' 0:2:400};%even trials
end
for sub=unique(fix.subject)
    c    = c+1;
    v{c} = {'subject' sub 'trialid' 1:2:400};%odd trials
end
fix.getmaps(v{:});
fix.maps   = imresize(fix.maps,.1,'method','bilinear');
corrmat    = fix.corr;
corrmat    = corrmat(1:N,(N+1):end);
%%
% same for phases seperately
N = length(unique(fix.subject));
phases = 1:5;
ed = nan(N,N,length(phases));
vecmaps=[];

for ph = phases
    v = [];
    c = 0;
    for sub=unique(fix.subject)
        c    = c+1;
        v{c} = {'subject' sub 'trialid' 0:2:400 'phase' ph};%even trials
    end
    for sub=unique(fix.subject)
        c    = c+1;
        v{c} = {'subject' sub 'trialid' 1:2:400 'phase' ph};%odd trials
    end
    fix.getmaps(v{:});
    fix.maps   = imresize(fix.maps,.1,'method','bilinear');
    vecmaps = [vecmaps fix.vectorize_maps];
end

%% subjects x phases
N = length(unique(fix.subject));
phases = 1:5;
vecmaps=[];
for ph = phases
    v = [];
    c = 0;
    for sub=unique(fix.subject)
        c    = c+1;
        v{c} = {'subject' sub 'phase' ph};%even trials
    end
    fix.getmaps(v{:});
    fix.maps   = imresize(fix.maps,.1,'method','bilinear');
    vecmaps = [vecmaps fix.vectorize_maps];
end
corrmat = corr(vecmaps);
Z = linkage(vecmaps','average','correlation');
[H,T,outperm]=dendrogram(Z);

tree = linkage(corrmat,'average','correlation')
D = pdist(corrmat,'correlation');
leafOrder = optimalleaforder(tree,D);

subplot(2,1,1);[H,T,outperm] = dendrogram(tree,0);
subplot(2,1,2);[H1,T1,outperm1] = dendrogram(tree,0,'Reorder',leafOrder);

%colorcoding
rgb = circshift(hsv(27),[3 0]);
subcol=repmat(rgb,5,1);
phcol=[];for n=1:5;phcol = [phcol; repmat([(6-n)/5 (6-n)/5 (6-n)/5],27,1)];end
comb = phcol.*subcol;
fgcol = [repmat([1 0 0],27,1); repmat([0 1 0],27,1); repmat([0 0.6 0],27,1);repmat([0 0.3 0],27,1);repmat([1 0 0],27,1)];
for n=1:135;hold on;...
        plot(n,4,'o','MarkerFaceColor',subcol(n,:),'Color',subcol(n,:));...
        plot(n,3,'o','MarkerFaceColor',phcol(n,:),'Color',phcol(n,:));...
        plot(n,2,'o','MarkerFaceColor',comb(n,:),'Color',comb(n,:));...
end
hold on;
for n=1:135;plot(n,1,'o','MarkerFaceColor',subcol(outperm1(n),:),'Color',subcol(outperm1(n),:));hold on;end
for n=1:135;plot(n,0,'o','MarkerFaceColor',phcol(outperm1(n),:),'Color',phcol(outperm1(n),:));hold on;end
for n=1:135;plot(n,-1,'o','MarkerFaceColor',comb(outperm1(n),:),'Color',comb(outperm1(n),:));hold on;end
ylim([-1 5])

%% ET_discr and PMF - subjects
%
p = Project;
mask = p.getMask('ET_discr');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('PMF');
subjects = intersect(find(sum(mask,2)==4),subjects);
fix = Fixmat(subjects,1);
g = Group(subjects);
[mat tags ]= g.parameterMat;

[~,ind] = sort(mean(mat(:,[1 3]),2));
%% ET_feargen and SI - subjects
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('RATE');
subjects = intersect(find(mask),subjects);
fix = Fixmat(subjects,4);
g = Group(subjects);
g.getSI(3);
[mat tags ]= g.parameterMat;
[~,ind] = sort(mat(:,14));
%% compute dendrograms, plot fixationmaps acc. to leafOrders and sorted by ind
N = length(unique(fix.subject));
v = [];
c = 0;
for sub=unique(fix.subject)
    c    = c+1;
    v{c} = {'subject' sub};
end
fix.getmaps(v{:});
fix.maps   = imresize(fix.maps,.1,'method','bilinear');
vecmaps = fix.vectorize_maps;
%corrmat = corr(vecmaps);
tree = linkage(vecmaps','average','correlation');
D = pdist(vecmaps','correlation');
leafOrder = optimalleaforder(tree,D);

% subplot(2,1,1);[H,T,outperm] = dendrogram(tree,0);
[H1,T1,outperm1] = dendrogram(tree,0,'Reorder',leafOrder);o
%% sort by mean SCR-response
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('SCR');
subjects = intersect(find(sum(mask,2)==3),subjects);
fix = Fixmat(subjects,4);
fix.getsubmaps;

%% find y-axis cut threshold for given number k of clusters in dendrogram
k = 20;
clusters = cluster(tree,'maxclust',k);
t = sort(tree(:,3));
th = t(size(tree,1)+2-k);
[H,T,order] = dendrogram(tree,0,'reorder',order,'colorthreshold', th);
%% kmeans
%take all ET subjects, collect single fixations.
%Run k-means for a range of k
IDX = [];
data = fix.vectorize_maps';
for k=1:20%length(unique(fix.subject))
    fprintf('Running k = %d...\n',k)
    IDX(:,k)=kmeans(data,k,'Distance','correlation');
    [S,H] = silhouette(data, IDX(:,k));
    close all
    sil(k)=mean(S); %The mean silhoette value for two groups
end
figure;
plot(1:k, sil,'ok-','MarkerFaceColor','k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dendrogram feargen phase 4
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('SCR');
subjects = intersect(find(sum(mask,2)==3),subjects);
mask = p.getMask('RATE');
subjects = intersect(find(mask),subjects);
mask = p.getMask('PMF');
subjects = intersect(find(mask),subjects);

g = Group(subjects);
g.getSI(3);
[mat tags ] = g.parameterMat;
fix = Fixmat(subjects,4);
fix.getsubmaps
fix.maps   = imresize(fix.maps,.1,'method','bilinear');
% dendrogram with colored clusters
figure;cccc
[H,T,order,tree] = fix.dendrogram;


cluster(1).subs = [10 17 2 4 6 7 19 22 16 21 8 18 12 13];
cluster(2).subs = [15 11 1 9 3 5 20 14];
cluster(1).leafs = [1 2 3 4 5 6 8 9 13 14 17 18 20]; 
cluster(2).leafs = [7 10 11 12 15 16 19];
%% dendrogram discrimination phase 1
p = Project;
mask = p.getMask('ET_discr');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('SCR');
subjects = intersect(find(sum(mask,2)==3),subjects);
mask = p.getMask('RATE');
subjects = intersect(find(mask),subjects);
mask = p.getMask('PMF');
subjects = intersect(find(mask),subjects);

g = Group(subjects);
g.getSI(3);
[mat tags ] = g.parameterMat;
fix = Fixmat(subjects,1);
fix.getsubmaps
fix.maps   = imresize(fix.maps,.1,'method','bilinear');
% dendrogram with colored clusters
[H,T,order,tree] = fix.dendrogram;

cluster(1).subs = [6 10 15 3 21 11 4 1 12 16 18];
cluster(2).subs = [7 14 5 2 23 17 8 19 9 22 13 20];
cluster(1).leafs = [4 7 9 12 14 16 17 18 20 21]; 
cluster(2).leafs = [1 2 3 5 6 8 10 11 13 15 19];
%% subjects that we have all the data for
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('ET_discr');
subjects = intersect(find(mask),subjects);
mask = p.getMask('SCR');
subjects = intersect(find(sum(mask,2)==3),subjects);
mask = p.getMask('RATE');
subjects = intersect(find(mask),subjects);
mask = p.getMask('PMF');
subjects = intersect(find(sum(mask,2)==4),subjects);

g = Group(subjects);
g.getSI(3);
[mat tags ] = g.parameterMat;
fix = Fixmat(subjects,1:5);
fix.getsubmaps
fix.maps   = imresize(fix.maps,.1,'method','bilinear');
% dendrogram with colored clusters
[H,T,order,tree] = fix.dendrogram;

cluster(1).subs = [13 4 2 8 1 10 19 9 14];
cluster(2).subs = [12 3 21 15 5 16 16 18 20 7 17 11];
cluster(1).leafs = [4 6 11 12 15 16 17 19]; 
cluster(2).leafs = [1 2 3 5 7 8 9 10 13 14 18];
%OR
cluster(1).subs = [4 8 13 2 1 10];
cluster(2).subs = [12 3 5 19 9 14 16];
cluster(3).subs = [6 15 21 17 7 20 11 18];
cluster(1).leafs = [7 10 12 15 18]; 
cluster(2).leafs = [3 5 11 14 16 17];
cluster(3).leafs = [1 2 4 6 8 9 13];
%% plot the cluster parameters.
k = size(cluster,2);
rgbs = [1 0 0;0 0 1];%rgbs =hsv(k);
for cl = 1:k;
    for n = cluster(cl).leafs
        set(H(n),'Color',rgbs(cl,:))
    end
end
set(H(end),'Color','k')
set(H,'LineWidth',2)
axis off

fix.plotband(subjects(order));

% grouped Fixmaps for clusters from dendrogram
figure;
c=0;
v=[];
for cl = 1:k
    c=c+1;
    fix.getmaps({'subject' subjects(cluster(cl).subs) 'deltacsp' fix.realcond})
    h = subplot(1,k,c);imagesc(fix.maps);
    axis square
    axis off
    subplotChangeSize(h,.01,.01);
    caxis([0 8e-5])
    colorbar off
end
% SI mean barplots
figure;
for n =1:k
    subplot(1,k,n);
%     b=bar(mean(g.SI(cluster(n).subs)));
    b=bar(mean(g.sigma_test(cluster(n).subs)));
    hold on;
%     e = errorbar(mean(g.SI(cluster(n).subs)),std(g.SI(cluster(n).subs))./sqrt(length(cluster(n).subs)),'.-');
    e = errorbar(mean(g.sigma_test(cluster(n).subs)),std(g.sigma_test(cluster(n).subs))./sqrt(length(cluster(n).subs)),'.-');
    set(b,'FaceColor',rgbs(n,:));
    set(e,'LineWidth',1.5,'Color','k')
    axis square
    ylim([0 60])
    set(gca,'XTickLabel',{''},'FontSize',14)
end

% alpha mean bar plot
alpha = mean(mat(:,[1 3]),2);
figure;
for n =1:k
    subplot(1,k,n);
    b=bar(mean(alpha(cluster(n).subs)));
    hold on;
    e = errorbar(mean(alpha(cluster(n).subs)),std(alpha(cluster(n).subs))./sqrt(length(cluster(n).subs)),'.-');
    set(b,'FaceColor',rgbs(n,:));
    set(e,'LineWidth',1.5,'Color','k')
    axis square
    ylim([0 70])
    set(gca,'XTickLabel',{''},'FontSize',14)
end

%SCR
g.ModelSCR('test$',3)
for n = 1:length(g.ids)-1
    ind = g.subject{n}.scr.findphase('test$');
    g.subject{n}.scr.cut(ind);
    g.subject{n}.scr.run_ledalab;
    scr(:,:,n) = g.subject{n}.scr.ledalab.mean(1:800,1:8);
end
for n=[cluster(1).subs cluster(2).subs];subplot(5,5,n);plot(repmat(g.subject{n}.scr.ledalab.x(1:800,1),[1 8]),scr(:,:,n));end
%exemplary SCRs
subs = [2 3];%4 is for cluster 1, 3 for cluster 2
figure;
for n = 1:k
subplot(1,k,n)
plot(repmat(g.subject{subs(n)}.scr.ledalab.x(1:800,1),[1 8]),mean(scr(:,:,subs(n)),3),'LineWidth',2.5)
axis square
set(gca,'FontSize',14)
xlim([-1.1 7.1])
end
EqualizeSubPlotYlim(1)
%Group SCR Curve
figure;
for n = 1:k
subplot(1,k,n)
t0 = find(g.subject{1}.scr.ledalab.x(:,1)==0);
const     = squeeze(mean(mean(mean(scr(1:t0,:,cluster(n).subs),1),2),3));
plot(repmat(g.subject{1}.scr.ledalab.x(1:800,1),[1 8]),mean(scr(:,:,cluster(n).subs),3),'LineWidth',2.5)
axis square
set(gca,'FontSize',14)
xlim([-1.1 7.1])
end
EqualizeSubPlotYlim(1)
% SCR mean bar plots
% av_resp = squeeze(mean(mean(scr,1),2));
%OR
% const     = squeeze(mean(mean(scr(1:t0,:,:),1),2));
% av_resp = squeeze(mean(mean(scr,1),2))-const;
% %OR
av_resp = squeeze(mean(scr,1));
av_resp = squeeze(av_resp(4,:)./av_resp(8,:));

figure;
for n =1:k
    subplot(1,k,n);
    b=bar(mean(av_resp(cluster(n).subs)));
    hold on;
    e = errorbar(mean(av_resp(cluster(n).subs)),std(av_resp(cluster(n).subs))./sqrt(length(cluster(n).subs)),'.-');
    set(b,'FaceColor',rgbs(n,:));
    set(e,'LineWidth',1.5,'Color','k')
    axis square
%     ylim([0 1.5])
    set(gca,'XTicklabel',{''},'YTick',[0 .5 1 1.5],'FontSize',14)
end
%% random useful code snippets
valid = prod([g.tunings.rate{3}.pval;g.tunings.rate{4}.pval] > -log10(0.05));
mat((ismember(unique(fix.subject),subjects)==0),[1:11]) = nan;
mat((ismember(unique(fix.subject),subjects)==0),[12:14])= nan;
%% latest dendrogram, k=3 clusters
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_600);
mask = p.getMask('PMF');
subjects = intersect(find(sum(mask,2)==4),subjects);
fix = Fixmat(subjects,1);
fix.getsubmaps;
fix.maps        = imresize(fix.maps,0.1,'method','bilinear');

g = Group(subjects);[mat tags]=g.parameterMat;mises = g.loadmises;clear g;mat = [mat mises];

%now we need to prepare the correct data for these guys
mask = p.getMask('RATE');
invalid_r = ~ismember(subjects,find(mask));
mask = p.getMask('PMF');
invalid_a = ~ismember(subjects,find(sum(mask,2)==4)); 
mat(invalid_a,1:11) = NaN;
mat(invalid_r,12:end)= NaN;

for nc = 11
    [branch_id,ids] = fix.dendrogram(2,mean(mat(:,nc),2));
    supertitle(tags(nc),1,'interpreter','none');
end
%[h p stats] = ttest2(mat(branch_id==1,14),mat(branch_id==2,14)) % gives p = 0.045
%% plot the relevant parameters for cluster 1 and 2
%bars alpha
subplot(4,2,1:2);
for n=1:2
bar(n,nanmean(mean(mat(branch_id==n,[1 3]),2)));
hold on;
errorbar(n,nanmean(mean(mat(branch_id==n,[1 3]),2)),nanstd(mean(mat(branch_id==n,[1 3]),2))./sqrt(sum(branch_id==n)),'k.','LineWidth',2);
end
xlim([0.5 2.5])
ylabel('initial alpha')
ylim([40 70])
set(gca,'XTick',[],'XTicklabel',[]);box off;
%bars kappa
subplot(4,2,3:4);
for n=1:2
bar(n,nanmean(mat(branch_id==n,13)));
hold on;
errorbar(n,nanmean(mat(branch_id==n,13)),nanstd(mat(branch_id==n,13))./sqrt(sum(branch_id==n)),'k.','LineWidth',2);
end
xlim([0.5 2.5]);box off;
ylabel('kappa test')
set(gca,'XTick',[])
%SCR Graphs
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_fordendroN27.mat')
pl(1)=subplot(4,2,5);
SetFearGenColors; plot(nanmean(scr_datac(:,1:8,branch_id==1),3),'LineWidth',2);axis square
xlim([0 800]);ylim([-.2 1]);box off;set(gca,'XTick',[]);
ylabel('SCR [\muS]')
pl(2)=subplot(4,2,6);
SetFearGenColors; plot(nanmean(scr_datac(:,1:8,branch_id==2),3),'LineWidth',2);axis square
xlim([0 800]);ylim([-.2 1]);box off;set(gca,'XTick',[]);

scr_datac = scr_data(:,1:8,:)-repmat(scr_data(:,9,:),[1 8 1]);
%SCR tunings
pl(3)=subplot(4,2,7);
b=bar(1:8,squeeze(nanmean(scr_bars(branch_id==1,1:8,3),1)));SetFearGenBarColors(b);
hold on;
errorbar(1:8,squeeze(nanmean(scr_bars(branch_id==1,1:8,3),1)),...
    squeeze(nanmean(scr_bars(branch_id==1,1:8,3),1))./sqrt(sum(branch_id==1)),'k.','LineWidth',1.5)
ylim([0 1]);
ylabel('SCR [\muS]')
axis square; box off
set(gca,'XTick',[4 8],'XTicklabel',{'CS+' 'CS-'})
pl(4)=subplot(4,2,8);
b=bar(1:8,squeeze(nanmean(scr_bars(branch_id==2,1:8,3),1)));SetFearGenBarColors(b);
hold on;
errorbar(1:8,squeeze(nanmean(scr_bars(branch_id==2,1:8,3),1)),...
    squeeze(nanmean(scr_bars(branch_id==2,1:8,3),1))./sqrt(sum(branch_id==2)),'k.','LineWidth',1.5)
axis square;box off
ylim([0 1]);
set(gca,'XTick',[4 8],'XTicklabel',{'CS+' 'CS-'})





%%
% load('C:\Users\user\Dropbox\feargen_hiwi\dump\m.mat')
% load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scrmat.mat')
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_fordendro.mat','subjects','branch_id','scr_data','scr_datac','scr_bars_dendro')
sub1 = subjects(branch_id==1);
sub2 = subjects(branch_id==2);
ph = 2;
scr_crit = (scr_bars_dendro(:,4,ph)-scr_bars_dendro(:,8,ph))./(scr_bars_dendro(:,4,ph)+scr_bars_dendro(:,8,ph));
%or
scr_crit = nanmean(scr_bars_dendro(:,1:8,ph),2)-scr_bars_dendro(:,9,ph);
%or
scr_crit = mean(scr_bars_dendro(:,1:8,ph),2);

[branch_id,order] = fix.dendrogram(3,scr_crit);
[h,p,ci,stats] = ttest2(scr_crit(branch_id==1),scr_crit(branch_id==2),'tail','left');
% or as bars...
clf
for n=1:2
    subplot(1,2,n);bar(1,mean(scr_crit(branch_id==n)));hold on;errorbar(1,mean(scr_crit(branch_id==n)),std(scr_crit(branch_id==n))./sqrt(sum(branch_id==n)),'k.');
    axis square
end
EqualizeSubPlotYlim(gcf)
%% dendrogram of phases 1:5
p = Project;
mask = p.getMask('ET_feargen');
subjects = intersect(find(mask),Project.subjects_1500);
mask = p.getMask('ET_discr');
subjects = intersect(find(mask),subjects);
fix = Fixmat(subjects,1:5);
g = Group(subjects);[mat tags] = g.parameterMat;mises = g.loadmises; mat(:,12:16) = mises;
mask = p.getMask('RATE');
invalid_r = ~ismember(subjects,find(mask));
mask = p.getMask('PMF');
invalid_a = ~ismember(subjects,find(sum(mask,2)==4)); 
mat(invalid_a,1:11) = NaN;
mat(invalid_r,12:14)= NaN;


fix.getsubmaps; 
fix.maps        = imresize(fix.maps,0.1,'method','bilinear');

load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\scr_1500BCT.mat')
scr_crit = (scr_bars(:,4,3)-scr_bars(:,8,3))./(scr_bars(:,4,3)+scr_bars(:,8,3));
[branch_id,ids] = fix.dendrogram(2,scr_crit);


%plot SCR tunings for both clusters
subplot(1,2,1); b=bar(nanmean(scr_bars(branch_id==1,1:8,3)));axis square;SetFearGenBarColors(b);hold on;
errorbar(nanmean(scr_bars(branch_id==1,1:8,3)),nanstd(scr_bars(branch_id==1,1:8,3))./sqrt(length(branch_id)),'k.');
subplot(1,2,2); b=bar(nanmean(scr_bars(branch_id==2,1:8,3)));axis square;SetFearGenBarColors(b);hold on;
errorbar(nanmean(scr_bars(branch_id==2,1:8,3)),nanstd(scr_bars(branch_id==1,1:8,3))./sqrt(length(branch_id)),'k.');
EqualizeSubPlotYlim(gcf)
%correct for mean arousal
figure
scr_bars_mc = scr_bars(:,1:8,3)-repmat(nanmean(scr_bars(:,1:8,3),2),[1 8]);
subplot(1,2,1); b=bar(nanmean(scr_bars_mc(branch_id==1,:)));axis square;SetFearGenBarColors(b);hold on;
errorbar(nanmean(scr_bars_mc(branch_id==1,:)),nanstd(scr_bars_mc(branch_id==1,:))./sqrt(length(branch_id)),'k.');
subplot(1,2,2); b=bar(nanmean(scr_bars_mc(branch_id==2,:)));axis square;SetFearGenBarColors(b);hold on;
errorbar(nanmean(scr_bars_mc(branch_id==2,:)),nanstd(scr_bars_mc(branch_id==1,:))./sqrt(length(branch_id)),'k.');
EqualizeSubPlotYlim(gcf)
% SCR after cond (only 4 and 8)
figure;
subplot(1,2,1); b=bar([4 8],nanmean(scr_bars(branch_id==1,[4 8],2)));axis square;xlim([0 9]);
subplot(1,2,2); b=bar([4 8],nanmean(scr_bars(branch_id==2,[4 8],2)));axis square;xlim([0 9]);
EqualizeSubPlotYlim(gcf)
% SCR curves
scr_data = NaN(800,11,length(subjects));
sc=0;
for sub = subjects(:)'
        fprintf('Working on subject %d .. \n',sub)
        sc=sc+1;
        try
        s = Subject(sub);
        s.scr.cut(s.scr.findphase('test$'));
        s.scr.run_ledalab;
        scr_data(:,:,sc) = s.scr.ledalab.mean;
        end
end

figure;
subplot(1,2,1);SetFearGenColors; plot(nanmean(scr_data(:,1:8,branch_id==1),3),'LineWidth',2);axis square
subplot(1,2,2);SetFearGenColors; plot(nanmean(scr_data(:,1:8,branch_id==2),3),'LineWidth',2);axis square
EqualizeSubPlotYlim(gcf)


%% find beneficial locations that correlate with discrimination/alpha
clear all
p = Project;
mask = p.getMask('ET_discr');
subjects = intersect(find(mask),[Project.subjects_1500]);
mask = p.getMask('PMF');
subjects = intersect(find(sum(mask,2)==4),subjects);
g = Group(subjects);
[mat tags] = g.parameterMat;
clear g

fix = Fixmat(subjects,1);
fix.getsubmaps;
fix.maps        = imresize(fix.maps,0.1,'method','bilinear');
submaps = fix.vectorize_maps;
r = NaN(length(submaps),1);
pval = NaN(length(submaps),1);
for n = 1:length(submaps)
    [r(n),pval(n)] = corr(submaps(n,:)',mean(mat(:,[1 3]),2));
end
fix.maps = reshape(r,[50 50]);
fix.plot

%% freezing behavior?
clear all
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
fix = Fixmat(subjects,[2 4]);
%
[a,b] = fix.histogram;
close all;
for pc = [1 2]%phase count
histmat(:,:,pc) = zscore(a(:,1:8,pc)')';
subplot(1,2,pc)
b = bar(mean(histmat(:,:,pc)));
hold on;
SetFearGenBarColors(b);
e = errorbar(mean(histmat(:,:,pc)),std(histmat(:,:,pc))./sqrt(length(histmat)),'k.');
axis square
set(gca,'XTicklabel',{'' '' '' 'CS+' '' '' '' 'CS-'})
end
subplot(1,2,1)
ylabel('mean zscore(nfix)')
%%
figure;
subplot(1,2,1);imagesc(a(:,1:8));
subplot(1,2,2);imagesc(histmat);

%sort histmat by SCR amplitude
scr_bars(2,:)=[];
histmat(end,:)=[];
a(end,:) = [];
scr_ampl = scr_bars(:,1) - scr_bars(:,2);
[~,i] = sort(scr_ampl);
subplot(1,2,1);imagesc(a(i,1:8));
subplot(1,2,2);imagesc(histmat(i,:));

%% SVM results (svm_analysis option 14 and 15)
%insubject:
w0 = w;
w = mean(w0,3);
num=60;
% fix.getsubmaps;
% fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
% meanmap = Scale(mean(fix.maps,3));

for n = 1:size(w,2);
    hp(:,:,n) = squeeze(eigen(:,1:num)*w(:,n));
end
% for n = 1:27;fix.maps(:,:,n) = reshape(hp(:,n),[50 50]);end
for n = 1:27;fix.maps(:,:,n) = reshape(hp(:,n),[50 50]).*meanmap;end
for n = 1:27;fix.maps(:,:,n) = reshape(hp(:,n),[50 50]);end
fix.plot

a = mean(result,4);
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,2]);
mean(scaled,3);


%% take difference CS+ - CS- as simple fixation map
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
fix = Fixmat(subjects,2:4);

subjects = unique(fix.subject);
pc = 0;
for ph = unique(fix.phase)
    pc = pc+1;
    sc=0;
    for sub = subjects(:)'
        sc=sc+1;
        c = 0;
        v = [];
        for cs = [0 180]
            c    = c+1;
            v{c} = {'deltacsp' cs 'subject' sub 'phase' ph};
        end
    fix.getmaps(v{:});
    maps(:,:,sc,pc) = fix.maps(:,:,1) - fix.maps(:,:,2);%-mean(fix.maps,3);
    end
end
fix.maps = squeeze(mean(maps,3));%over subjects
fix.plot
%weights from overall fixation probability within the phase...
c = 0;
v=[];
for ph = unique(fix.phase)
    c    = c+1;
    v{c} = {'phase' ph};
end
fix.getmaps(v{:});
meanmaps = Scale(fix.maps);
fix.maps = squeeze(mean(maps,3)).*meanmaps;
%% take difference CS+ NEIGHBORS- CS- as simple fixation map
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
fix = Fixmat(subjects,[2 4]);

pc = 0;
for ph = unique(fix.phase)
    pc = pc+1;
    sc=0;
    for sub = subjects(:)'
        sc=sc+1;
        c = 0;
        v = [];
        for cs = fix.realcond(:)' %[-45 45; 180 180]';neighbor vs CS-
            c    = c+1;
            v{c} = {'deltacsp' cs 'subject' sub 'phase' ph};
        end
    fix.getmaps(v{:});
    maps(:,:,sc,pc) = mean(fix.maps(:,:,[3 5]),3) - fix.maps(:,:,8);
    end
end
fix.maps = squeeze(mean(maps,3));%over subjects
fix.plot
%weights from overall fixation probability within the phase...
c = 0;
v=[];
for ph = unique(fix.phase)
    c    = c+1;
    v{c} = {'phase' ph};
end
fix.getmaps(v{:});
meanmaps = Scale(fix.maps);
fix.maps = squeeze(mean(maps,3)).*meanmaps;

%% is CS+ more similar between subjects than CS-?
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
fix = Fixmat(subjects,4);

c = 0;
v = [];
for cond = [0 180]
    for sub = subjects(:)'
        c    = c+1;
        v{c} = {'deltacsp' cond 'subject' sub};
    end
end
fix.getmaps(v{:});
fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
maps = fix.vectorize_maps;
for a = 1:27
    for b = 1:27
        if a < b
            rcsp(a,b) = corr(maps(:,a),maps(:,b));
        else
            rcsp(a,b) = NaN;
        end
    end
end
nanmean(nanmean(rcsp))
for a = 28:54
    for b = 28:54
        if a < b
            rcsn(a-27,b-27) = corr(maps(:,a),maps(:,b));
        else
            rcsn(a-27,b-27) = NaN;
        end
    end
end
nanmean(nanmean(rcsn))
%%
p = Project;
subjects = [Project.subjects_1500 Project.subjects_600];

c=0;
for n = subjects(:)'
    c=c+1;
    s = Subject(n);
    resp(c,1) = sum(s.paradigm{2}.out.response(s.paradigm{2}.presentation.oddball));
    resp(c,2) = sum(s.paradigm{3}.out.response(s.paradigm{3}.presentation.oddball));
    resp(c,3) = sum(s.paradigm{4}.out.response(s.paradigm{4}.presentation.oddball));
    clear s
end
%% entropy
%% try the same as MDS
tsub = length(unique(fix.subject));
subc            = 0;
for subject = unique(fix.subject);
    subc = subc + 1
    %creaete the query cell
    v = [];
    c = 0;
    for ph = [4]
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject' subject};
        end
    end
    fix.getmaps(v{:});
    entr(:,subc) = fix.entropy;
end
%% mean correlation of subjects, ,across phases
clear r; clear r_z; clear r_back;
subc            = 0;
for subject = unique(fix.subject);
    subc = subc + 1
    %creaete the query cell
    v = [];
    c = 0;
    for ph = 1:5
        c    = c+1;
        v{c} = {'phase', ph, 'subject' subject};
    end
    fix.getmaps(v{:});
    r = fix.corr;
    r = triu(r);
    r = CancelDiagonals(r,0);
    r_z = fisherz(r(r~=0));
    r_back(subc) = ifisherz(mean(r_z));
end
%% mean correlation of phases, across subjects
pc            = 0;
for ph = unique(fix.phase);
    pc = pc + 1
    %creaete the query cell
    v = [];
    c = 0;
    for sub = unique(fix.subject)
        c    = c+1;
        v{c} = {'phase', ph, 'subject' sub};
    end
    fix.getmaps(v{:});
    r = fix.corr;
    r = triu(r);
    r = CancelDiagonals(r,0);
    r_z = fisherz(r(r~=0));
    r_back(pc) = ifisherz(mean(r_z));
end

%% loadings of discrmap x cond
disc = rc;
disc(disc>=0) = 0;
disc(isnan(disc)) = 0;
disc = -disc;
%create the query cell for conditions
sc = 0;
for sub = unique(fix.subject)
    sc = sc+1;
    v = [];
    c = 0;
    for ph = 4
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject' sub};
        end
    end
    fix.getmaps(v{:});
    fix.maps = fix.maps;% - repmat(mean(fix.maps,3),[1 1 8]);
    fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
    condmaps = fix.vectorize_maps;
    discrload(:,sc) = disc * condmaps;
end
bar(mean(discrload,2))
%% only in good people (from dendrogram)
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\disc_posmap.mat')
load('C:\Users\user\Documents\Experiments\FearCloud_Eyelab\data\midlevel\dendro_branch_id_subjects.mat','sub2')
fix = Fixmat(sub2,4);
sc = 0;
for sub = unique(fix.subject)
    sc = sc+1;
    v = [];
    c = 0;
    for ph = 4
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject' sub};
        end
    end
    fix.getmaps(v{:});
    fix.maps = fix.maps;% - repmat(mean(fix.maps,3),[1 1 8]);
    fix.maps   = imresize(fix.maps,0.1,'method','bilinear');
    condmaps = fix.vectorize_maps;
    discrload(:,sc) = disc * condmaps;
end
bar(mean(discrload,2))

%% compare cs+ to cs- (only good people)
sc = 0;
for sub = unique(fix.subject)
    sc = sc+1;
    v = [];
    c = 0;
    for cond = [0 180]
        c    = c+1;
        v{c} = {'deltacsp', cond, 'subject', sub};
    end
    fix.getmaps(v{:});
    maps(:,:,sc) = fix.maps(:,:,1) - fix.maps(:,:,2);
end
fix.maps = mean(maps,3);
fix.plot
%% 5 x 27 Matrix - subject similarity within phases
clear all
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
fix = Fixmat(subjects,1:5);

pc = 0;
for ph = unique(fix.phase)
    pc = pc+1;
    v = [];
    c=0;
    for sub = subjects(:)'
        c    = c+1;
        v{c} = {'deltacsp' fix.realcond 'subject' sub 'phase' ph}; % like this, we have 27x27x5
    end
    fix.getmaps(v{:});
    rmat(:,:,pc) = fix.corr;
    rmat_f(:,:,pc) = reshape(fisherz(fix.corr),[27 27]);
end
ind = logical(CancelDiagonals(triu(ones(27,27)),0));
subplot(1,2,1)
imagesc(rmat(:,:,1))
axis square
subplot(1,2,2)
ind = logical(CancelDiagonals(triu(ones(27,27)),0));
imagesc(ind)
axis square
for n = 1:5
    dummy = rmat_f(:,:,n);
    corrs(:,n) = dummy(ind);
end
ifisherz(mean(mean(corrs)))%should give .53 
figure; bar(ifisherz(mean(corrs)))

%% phase similarity within subjects
clear all
p               = Project;
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
fix = Fixmat(subjects,1:5);

sc = 0;
c=0;
for sub = subjects(:)'
    sc = sc+1;
    pc=0;
    v=[];
    for ph = unique(fix.phase)
            pc    = pc+1;
            v{pc} = {'deltacsp' fix.realcond 'subject' sub 'phase' ph};
    end
 fix.getmaps(v{:});
%  fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 5]);
 rmat(:,:,sc) = fix.corr;
end

%% plot this mat
clf
cd('C:\Users\Lea\Documents\GitHub\globalfunctions\')
avemat = fisherz_inverse(mean(fisherz(rmat),3));
imagesc(avemat)
axis square
box off
cmap = zeros(256,3);
cmap(:,1) = linspace(0,1,256);
colormap(cmap)
caxis([min(avemat(:)) max(avemat(:))])
labels = {'Discrimination pre','Free viewing','Conditioning','Generalization','Discrimination post'};
set(gca,'XTick',1:5,'YTick',1:5,'YTickLabel',labels,'XTickLabel',{'D','F','C','G','D'},'FontSize',12)
cbh = colorbar;
set(cbh,'YTick',.8:.1:1,'FontSize',12)
%% compute means
rmat_f  = nan(5,5,27);
for n = 1:27; 
    rmat_f(:,:,n) = reshape(fisherz(rmat(:,:,n)),[5 5 1]);
end
% rmat_mean = CancelDiagonals(reshape(ifisherz(mean(rmat_f,3)),[5 5]),1); % average across the subjects, take the mean, inv fisher it...
% imagesc(rmat_mean)

ind = logical(CancelDiagonals(triu(ones(5,5)),0));
for n = 1:27;
    dummy = rmat(:,:,n);
    corrs(:,n) = dummy(ind);%this is then size (sum(1:5)-5 x 27)
end

ifisherz(mean(mean(fisherz(corrs)))) % gives r = .84 for un-cocktail-blank-corrected fixmaps, without diagonals
figure;bar(ifisherz(mean(fisherz(corrs))))
%% can this be related to feargen sharpness or sth?
gendiscr = squeeze(rmat(4,5,:));
load('C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\data\midlevel\misesmat.mat');
params = misesmat(unique(fix.subject),:);
fwhm = nan(27,1);
for n = 1:27; 
    if ~isnan(params(n,2))
    fwhm(n) = vM2FWHM(params(n,2));
    else
        fprintf('isnan, skipping  %d \n',n)
    end
end
[rho,pval] = corr(gendiscr(~isnan(fwhm)),fwhm(~isnan(fwhm)))