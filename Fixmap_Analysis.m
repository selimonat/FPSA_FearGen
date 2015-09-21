%% get the data first;
addpath('/Users/onat/Documents/Code/Matlab/oop/')
clear all;
load(sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep));
subjects = intersect(Project.subjects_1500,find(ETmask(:,1)==1));
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
    fix.maps             = fix.maps - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
    fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
    fix.plot;
    supertitle(mat2str(k),1);
    SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/600/baselineANDMeancorrection/Maps_Scale_%04.4g.png',k));    
    out= fix.corr;
    imagesc(out(9:end,9:16));colorbar
    SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/600/baselineANDMeancorrection/CorrMat_Scale_%04.4g.png',k));
    title(mat2str(k));        
end
%% For a given fwhm, compute fixmaps and cormat for different fixation indices.
fix = Fixmat(Project.subjects_1500,[2 3 4]);
%%
fix.kernel_fwhm = 20;
for nfix = [1 2 3 4 5 6];%there are about 200 7th fixations (not enough)
    %creaete the query cell
    v = [];
    c = 0;
    for ph = [2 4]
        for cond = -135:45:180
            c = c+1;
            v{c} = {'phase', ph, 'deltacsp' cond 'fix' nfix};
        end
    end
    % plot and save fixation maps       
    fix.getmaps(v{:});
    fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
    fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
    fix.plot;
    supertitle(mat2str(k),1);
    SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/1500/SingleFixations/baselineANDMeancorrection/Maps_nFix_%02d.png',nfix));
    out= fix.corr;
    imagesc(out(9:end,9:16));colorbar
    SaveFigure(sprintf('/Users/onat/Desktop/fearcloud/1500/SingleFixations/baselineANDMeancorrection/Corrmat_nFix_%02d.png',nfix));
    title([mat2str(fix.kernel_fwhm) '-' mat2str(nfix)] );
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
clear all;
load(sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep));
subjects = intersect(Project.subjects_1500,find(ETmask(:,1)==1));
fix             = Fixmat(subjects,[2 3 4]);

%% plot and save single subject fixations maps UNCORRECTED
tsub = length(unique(fix.subject));
for k = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
    fix.kernel_fwhm = k;
    covmat          = nan(16,16,tsub);
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
        fix.plot
        %SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/AllFixations/uncorrected/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
        saveas(gcf,sprintf('C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/uncorrected/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
        close all
    end
end
%% plot and save single subject fixations maps, BASELINE CORRECTED
tsub = length(unique(fix.subject));
for k = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
    fix.kernel_fwhm = k;    
    covmat          = nan(16,16,tsub);
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
        fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
        fix.plot
        %SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/AllFixations/baselinecorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
        saveas(gcf,sprintf('C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/baselinecorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
    end
end

%% plot and save single subject fixations maps, BASELINE CORRECTED
tsub = length(unique(fix.subject));
for k = Fixmat.PixelPerDegree*logspace(log10(.1),log10(2.7),20);
    fix.kernel_fwhm = k;    
    covmat          = nan(16,16,tsub);
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
        fix.maps             = fix.maps            - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
        fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average
        fix.plot
        %SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/AllFixations/baselinecorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
        saveas(gcf,sprintf('C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/baselineANDMeancorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
    end
end
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
%% Find tuned pixels
fix = Fixmat(Project.subjects_1500,[2 3 4]);
%%
x = -135:45:180;
v = [];
c = 0;
for ph = [2 4]
    for cond = -135:45:180
        c    = c+1;
        v{c} = {'phase', ph, 'deltacsp' cond};
    end
end
fix.kernel_fwhm = 20;
fix.maptype     = 'bin';
fix.getmaps(v{:});
fix.maps             = fix.maps - repmat(mean(fix.maps(:,:,1:8),3),[1 1 16]);%correct for baseline
fix.maps(:,:,9:16)   = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);%take out the average                    
%%
maps = fix.vectorize_maps;
i    = find(sum(diff(maps,1,2),2) ~= 0)
d.y  = maps(i,9:16);
d.x  = repmat(-135:45:180,size(d.y,1),1);
t    = Tuning(d);
t.SingleSubjectFit(3);
%%
d.x  = repmat(-135:45:180,10,1);
for n = 1:10
d.y(n,:)  = make_gaussian_fmri_zeromean(-135:45:180,randsample(linspace(10,100,10),1),randsample(linspace(20,150,10),1));
end

%%
%%good subject? bad subject?
load(sprintf('%smidlevel%ssubjmasks%sETmask.mat',Project.path_project,filesep,filesep));
subjects = intersect(Project.subjects_600,find((ETmask(:,1)==1)));
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
    t=supertitle('Good - Bad Base&Mean corrected (1500ms)',1);
    set(t,'FontSize',14)
  %% 
%SaveFigure(sprintf('/Users/onat/Desktop/fixationmaps/singlesubject/1500/AllFixations/baselinecorrection/Maps_Scale_%04.4g_Subject_%03d.png',k,subject));
saveas(gcf,'C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/good_vs_bad/within_group_corrected/1500_Good-Bad_k20_baseANDmeancorr.png');


%%
%%%improvement
%what is subject's improvement? For CSP? or mean(CSP/CSN)?

saveas(gcf,'C:/Users/onat/Desktop/Lea/fearcloud_after_Msc/Fixmats/singlesubject/1500/AllFixations/improvement_CSP/within_group_corrected/improved_k20_uncorrected.png');


