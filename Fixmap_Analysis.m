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
%% Find tuned pixels
fix         = Fixmat(Project.subjects_1500,[2 3 4]);
fix.unitize = 0;
%% collect fixations
maps        = [];
for subject =unique(fix.subject)
    v = [];
    c = 0;
    for ph = [2 4]
        for cond = -135:45:180
            c    = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'subject', subject};
        end
    end
    fix.getmaps(v{:});
%     fix.maps = fix.maps - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);
    dummy    = fix.vectorize_maps;
    dummy(:,1:8) = [];
    dummy    = demean(dummy')';
    maps     = cat(3,maps,dummy);
end
data.x             = repmat(-135:45:180,[size(maps,1),1,size(maps,3)]);
data.y             = maps;
t                  = Tuning(data);
t.SingleSubjectFit(3);
%%
t   = [];
c   = 0;
for fwhm = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20)
    c = c+1;
    fix.kernel_fwhm      = fwhm;
    fix.maptype          = 'conv';
    fix.getmaps(v{:});
    fix.maps             = fix.maps - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
    %
    maps                 = fix.vectorize_maps;
    d.y                  = maps(:,9:16);
    d.y                  = d.y - repmat(mean(d.y,2),[1 8]);
    d.x                  = repmat(x,size(d.y,1),1);
    t{c}         = Tuning(d);
    t{c}.SingleSubjectFit(3);
end
%%
d.x  = repmat(-135:45:180,10,1);
for n = 1:10
d.y(n,:)  = make_gaussian_fmri_zeromean(-135:45:180,randsample(linspace(10,100,10),1),randsample(linspace(20,150,10),1));
end

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

