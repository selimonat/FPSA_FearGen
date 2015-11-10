%% get the data first;
addpath('/Users/onat/Documents/Code/Matlab/oop/')
clear all;
p               = Project;
mask            = p.getMask('ET');
subjects        = find(mask(:,1));
subjects        = intersect(subjects,Project.subjects_1500);
fix             = Fixmat(subjects,[2 3 4]);
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
for k = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
    fix.kernel_fwhm = k;
    fix.getmaps(v{:});
    
    fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:16),3),[1 1 16]);%correct for baseline    
%     fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
    fix.plot;
    supertitle(mat2str(k),1);
    SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/MegaSubject/AllFixations/unitize/1500/nocorrection/Maps_Scale_%04.4g.png',k));    
    out= fix.corr;
    figure(1000);clf
    imagesc(out);colorbar
    SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/MegaSubject/AllFixations/unitize/1500/nocorrection/CorrMat_Scale_%04.4g.png',k));
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
%% plot and save single subject fixations maps
tsub = length(unique(fix.subject));
out  = [];
% for k = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
fix.kernel_fwhm = 25;
cormat          = nan(16,16,tsub);
subc            = 0;
for subject = unique(fix.subject);
    subc = subc + 1;
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
    fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:16),3),[1 1 16]);%correct for baseline
    %         fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
    %         fix.plot
    %         SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/AllFixations/baselinecorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
    %         saveas(gcf,sprintf('C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/baselineANDMeancorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
    cormat(:,:,subc) = fix.corr;
end
% end
%% compute single subject maps and similarity matrices, and finally average them across subjects.
figure(10);set(gcf,'position',[440    52   932   746]);clf
tsub = length(unique(fix.subject));
for k = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
    fix.kernel_fwhm = k;
    for nfix = [1 2 3 4 5 6];%there are about 200 7th fixations (not enough)
        covmat = nan(16,16,tsub);
        cormat = nan(16,16,tsub);
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
            figure(10);clf;
            % plot and save fixation maps
            fix.getmaps(v{:});
            fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
%             fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average            
            covmat(:,:,subc)     = fix.cov;
            cormat(:,:,subc)     = CancelDiagonals(fix.corr,NaN);
        end
        missing = sum(isnan(mean(mean(covmat),2)));
        fprintf('%03d subjects with no fixations.\n',missing);
        clf;
        subplot(2,2,1);
        imagesc(nanmean(covmat(1:8,1:8,:),3));thincolorbar('vert');;title('covariance');ylabel('Baseline');
        subplot(2,2,3);
        imagesc(nanmean(covmat(9:end,9:end,:),3));thincolorbar('vert');;title('covariance');ylabel('Test');
        subplot(2,2,2);
        imagesc(nanmean(cormat(1:8,1:8,:),3));thincolorbar('vert');;title('correlation');ylabel('Baseline');
        subplot(2,2,4);
        imagesc(nanmean(cormat(9:end,9:end,:),3));thincolorbar('vert'); ;title('correlation');ylabel('Test');
        supertitle(sprintf('Scale:%4.4g, Fix: %02d, Missing: %02d, Total: %02d',k,nfix,missing,tsub),1)
        drawnow;        
        SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/SingleFixations/baselinecorrection/Corrmats_Scale_%04.4g_Fix_%02d.png',k,nfix));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPARE FIXATION PATTERNS ACROSS THE 5 PHASES
p               = Project;
mask            = p.getMask('ET');
subjects        = find(mask(:,1));
subjects        = intersect(subjects,Project.subjects_1500);
fix             = Fixmat(subjects,[1 2 3 4 5]);
% PLOT FIXMAPS for different phases
v = [];
c = 0;
for ph = [1 2 3 4 5]    
        c = c+1;
        v{c} = {'phase', ph};    
end
fix.getmaps(v{:})
fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 5]);
fix.plot;

%% compare MEGA subject -> unitizing vs. Unitized Single subjects' average
maps = [];
for subs = unique(fix.subject);    
    subs
    v    = [];
    c    = 0;
    for ph = [1 2 3 4 5]
        c    = c+1;
        v{c} = {'phase', ph , 'deltacsp' [0 180 18000] 'subject' subs};
    end
    fix.getmaps(v{:});
    maps = cat(4, maps, fix.maps - repmat(mean(fix.maps,3),[1 1 size(fix.maps,3)]));    
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
fix      = Fixmat(subjects,4);
%
g        = Group(subjects);%get the data for these subjects
g.getSI(3);
[M i]       = g.parameterMat;
%%
param       = M(:,end);%sharpening index
%
c = 0;v = [];
for subs = subjects(:)'%good and then bad
    c       = c+1;
    v{c}    = {'phase', 4  'subject' subs 'deltacsp' 0};
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
%% Which parts of the fixation map correlates with an increased generalization.
p           = Project;
mask        = find(p.getMask('ET_feargen').*p.getMask('RATE'));
subjects    = intersect(Project.subjects_1500,mask);%subjects
g           = Group(subjects);%get the data for these subjects
g.getSI(3);
fix         = Fixmat(subjects,[2 3 4]);
fix.unitize = 1;
[M i]       = g.parameterMat;
param       = M(:,end);%sharpening index
i           = param > median(param);%median
%
c = 0;
v = [];
for subs = {subjects(i) subjects(~i)}%good and then bad
    c       = c+1;
    v{c}    = {'phase', 4 , 'deltacsp' [0] 'subject' subs{1}};
end
%
fix.kernel_fwhm = 15;
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
%% 
% %create the query cell
% v = [];
% c = 0;
% for subjects = {good bad}    
%         for cond = [0,18000]
%             c    = c+1;
%             v{c} = {'phase', [1 5], 'deltacsp' cond 'subject' subjects{1}};
%         end
% end
% % plot and save fixation maps
% fix.getmaps(v{:});
% fix.plot
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
fix.plot
%% phase fixation maps
v = [];
c = 0;
for ph = 1:5
        c    = c+1;
        v{c} = {'phase' ph 'subject' 6};
en
fix.getmaps(v{:});
fix.plot

%% SI x hyperplane
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\singletrialfixmaps\1500\SI_N24\labels.mat')
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\singletrialfixmaps\1500\SI_N24\SI_hp.mat')
g = Group(unique(labels.sub));
g.getSI(3);
[mat tags] = g.parameterMat;
fix = Fixmat(g.ids,4);
%get subj fixmaps from above, then...
fix.maps        = imresize(fix.maps,.1,'method','bilinear');
subjmap = fix.vectorize_maps;
for sub = 1:size(subjmap,2)
    hpload(sub) = SI_hp(:)'*subjmap(:,sub);
end
[r,p]=corr(hpload',g.SI)

% discr x hyperplane
fix = Fixmat(unique(labels.sub),1);
%get subj fixmaps from above,then...
fix.maps        = imresize(fix.maps,.1,'method','bilinear');
subjmap = fix.vectorize_maps;
for sub = 1:size(subjmap,2)
    hpload(sub) = discr_hp'*subjmap(:,sub);
end
SIsubs = unique(labels.sub); %load the SI analysis labels here..
[r,p] = corr(hpload(ismember(g.ids,SIsubs))',g.SI(ismember(g.ids,SIsubs)))