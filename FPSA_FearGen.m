function [varargout]=FPSA_FearGen(varargin);
% [varargout]=FPSA_FearGen(varargin);
%
% Complete analysis and figure generation pipeline for the FPSA manuscript,
% first version in Biorxiv located at http://biorxiv.org/content/early/2017/04/15/125682
%
% Using this code, it is possible to generate all the results and figures
% presented in the manuscript. Relies on Fancycarp toolbox (1), for dealing
% with fixation data. All dependencies can be installed by calling this
% function with 'download_project' argument.
%
% Requirements:
% A *nix system (tar, git, unzip, find), in case you want to download data
% with matlab using FPSA_FearGen('download_project') command.
% Matlab 2016.
% In case this option doesn't work for you, you can also download the project as a bundle
% from this link: https://osf.io/zud6h/
%
%
% Usage:
% VARARGIN sets an action related to an analysis, figure preparation or
% simple house keeping routines, such as data getters. For example,
% FPSA_FearGen('get_subjects') would return the subjects used for this
% manuscript. By analyzing that piece of code you can see how the selection
% procedure detailed in the manuscript is actually implemented.
% Some actions require inputs, these can be provided with 2nd, 3rd, so
% forth VARARGIN. For example, 'get_fpsa_fair' action requires two
% additional input arguments i/ fixations and ii/ runs. FIXATIONS
% determines which fixations to use to compute a dissimilarity matrix and
% RUNS determine pre- or post-learning phases. By convention baseline phase
% is denoted with 2, whereas Generalization phase by 4.
%
% Examples:
% FPSA_FearGen('get_subjects') retuns the indices to participants, included
%   in the analysis. These correspond to subXXX folder numbers.
% FPSA_FearGen('get_fixmat')
%   Returns the fixation data (cached in the midlevel folder).
% FPSA_FearGen('get_trialcount',2)
%   Returns number of trials per each subject and condition for baseline
%   phase (2). Use 4 for test phase.
% FPSA_FearGen('fix_counts',fixmat) counts the fixation density on 4
%   different ROIs used in the manuscript. Get FIXMAT with FPSA_FearGen('get_fixmat')
% FPSA_FearGen('get_fixmap',fixmat,12,1) returns FDMs using FIXMAT, for
%   participant 12 and fixation index 1.
% FPSA_FearGen('get_fixmap',fixmat,subjects,1:100) returns FDMs for
%   SUBJECTS and fixations from 1 to 100 (well, basically all the
%   fixations). Each FDM is a column vector.
% FPSA_FearGen('plot_fdm',maps) plots FDMs similar to Fig 3A.
% FPSA_FearGen('plot_ROIs') plots the ROIs.
% FPSA_FearGen('get_fpsa',1:100) would return the dissimilarity matrix
%   computed with exploration patterns containing all the fixations
%   (1:100). This is not a recommended method for computing dissimilarity
%   matrices, rather use:
% FPSA_FearGen('get_fpsa_fair',1:100,1:3) which would compute a fair
%   dissimilarity matrix for each run separately. 3 runs in generalization
%   are averaged later.
%
%
% INITIAL SETUP:
% Use FPSA_FearGen('download_project') to give it a start with it.
% Remember you will need basic unix tools for that such tar, git,
% unzip. This will download the data and the necessary scripts to your local
% machine and will add paths. Before starting you need to set the
% PATH_PROJECT variable below for your own liking.
%
% Contact: sonat@uke.de; lkampermann@uke.de

%% Set the default parameters
path_project         = sprintf('%s%s',homedir,'/Documents/Experiments/project_FPSA_FearGen/');% location of the project folder (MUST END WITH A FILESEP);
condition_borders    = {'' 1:8 '' 9:16};                                    % baseline and test condition labels.
block_extract        = @(mat,y,x,z) mat((1:8)+(8*(y-1)),(1:8)+(8*(x-1)),z); % a little routing to extract blocks from RSA maps
tbootstrap           = 1000;                                                % number of bootstrap samples
method               = 'correlation';                                       % methods for RSA computation
current_subject_pool = 0;                                                   % which subject pool to use (see get_subjects)
runs                 = 1:3;                                                 % which runs of the test phase to be used
criterion            ='strain' ;                                            % criterion for the MDS analysis.
force                = 0;                                                   % force recaching of results.
kernel_fwhm          = Fixmat([],[]).PixelPerDegree*.8;                            % size of the smoothing window (.8 degrees by default);
url                  = 'https://www.dropbox.com/s/0wix64zy2dlwh8g/project_FPSA_FearGen.tar.gz?dl=1';
%% overwrite default parameters with the input
invalid_varargin = logical(zeros(1,length(varargin)));
for nf = 1:length(varargin)
    if strcmp(varargin{nf}     , 'tbootstrap')
        tbootstrap           = varargin{nf+1};
    elseif strcmp(varargin{nf} , 'method')
        method               = varargin{nf+1};
    elseif strcmp(varargin{nf} , 'current_subject_pool')
        current_subject_pool = varargin{nf+1};
    elseif strcmp(varargin{nf} , 'runs')
        runs                 = varargin{nf+1};
    elseif strcmp(varargin{nf} , 'criterion')
        criterion            = varargin{nf+1};
    elseif strcmp(varargin{nf} , 'force')
        force                = varargin{nf+1};
    elseif strcmp(varargin{nf} , 'kernel_fwhm')
        kernel_fwhm          = varargin{nf+1};
    else
        invalid_varargin(nf) = true;%detect which varargins modify the default values and delete them
    end
end
varargin([find(~invalid_varargin) find(~invalid_varargin)+1]) = [];%now we have clean varargin cellarray we can continue

%%
if strcmp(varargin{1},'download_project');
    if ~exist(path_project)
        mkdir(path_project);
        mkdir([path_project 'tmp/'])
    end
    cd(path_project)
    %downloads the data and stimuli, download the code from github, and add
    %them to matlab path.
    %download data
    fprintf('Downloading data (this will take a while)...\n');
    tarfile              = ['/tmp/dummy.tar.gz'];
    s                    = urlwrite(url,tarfile);%download the data
    fprintf('Untarring the data...\n');
    untar(tarfile,'/tmp/');%untar it to the same location
    fprintf('Moving data to PATH_PROJECT...\n');
    movefile('/tmp/project_FPSA_FearGen/*',regexprep(path_project,'/$',''));%copy the contents of the file to PATH_PROJECT
    
    
    %download 3 other repositories, which we depend on here.
    fprintf('Downloading the analysis code and adding it to path...\n');
    cd(path_project);
    mkdir code
    cd code
    system(['git clone https://github.com/selimonat/fancycarp.git']);
    cd('./fancycarp');
    system(['git checkout FPSA_FearGen']);
    addpath(pwd)
    cd('..');
    system(['git clone https://github.com/selimonat/globalfunctions.git']);
    cd('./globalfunctions');
    addpath(pwd)
    cd ..
    
elseif strcmp(varargin{1},'get_subjects');
    %% returns subject indices based on the CURRENT_SUBJECT_POOL variable.
    % For the paper we use current_pool = 1, which discards all subjects:
    % who are not calibrated good enough +
    % who did not get the CS+ - UCS association.
    % Results are cached, use FORCE = 1 to recache. Set
    % CURRENT_SUBJECT_POOL = 0 to not select participants.
    
    filename = sprintf('%s/data/midlevel/subjectpool_%03d.mat',path_project,current_subject_pool);
    if exist(filename) == 0 | force
        if current_subject_pool == 0;
            subjects = Project.subjects(Project.subjects_ET);
        elseif current_subject_pool == 1%find tuned people;
            
            fprintf('finding tuned subjects first...\n');
            p=[];sub=[];pval=[];
            for n = Project.subjects(Project.subjects_ET);
                s    = Subject(n);
                s.fit_method = 8;%mobile vonMises function;
                p    = [p    ; s.get_fit('rating',4).params];
                pval = [pval ; s.get_fit('rating',4).pval];
                sub  = [sub  ; n];
            end
            valid    = (abs(p(:,3)) < 45) & pval > -log10(.05);%selection criteria
            fprintf('Found %03d valid subjects...\n',sum(valid));
            subjects = sub(valid);
            %             save(filename,'subjects');
        end
    else
        load(filename);
    end
    varargout{1} = subjects;
    allsubs      = Project.subjects(Project.subjects_ET);
    varargout{2} = allsubs(~ismember(allsubs,subjects));
elseif strcmp(varargin{1},'get_trialcount')
    %% Sanity check for number of trials per subject.
    % goes through subjects and counts the number of trials in a PHASE. The
    % output is [subjects conditions].
    % VARARGIN{2}.
    %
    % Example: FPSA_FearGen('get_trialcount',4) for phase 4 (test phase).
    % Example: FPSA_FearGen('get_trialcount',2) for phase 4 (test phase).
    phase  = varargin{2};
    fixmat = FPSA_FearGen('get_fixmat');
    s      = 0;
    for ns = unique(fixmat.subject)
        s = s+1;
        c = 0;
        for nc = unique(fixmat.deltacsp)
            c = c+1;
            C(s,c) = (length(unique(fixmat.trialid(fixmat.subject == ns & fixmat.deltacsp == nc & fixmat.phase == phase))));
        end
    end
    varargout{1} = C;
    imagesc(C(:,1:8));
    colorbar;
elseif strcmp(varargin{1},'get_fixmat');
    %% load the fixation data in the form of a Fixmat.
    %   For more information on Fixmat structure refer to (1), where this
    %   data is also published.
    %   Will return fixmat for the baseline and test phases. Generalization
    %   phase has 3 runs, by default all are returned.
    %   Use force = 1 to recache (defined at the top).
    %
    %   Example: fixmat = FPSA_FearGen('get_fixmat')
    
    filename = sprintf('%s/data/midlevel/fixmat_subjectpool_%03d_runs_%03d_%03d.mat',path_project,current_subject_pool,runs(1),runs(end));
    fix      = [];
    if exist(filename) == 0 | force
        subjects = FPSA_FearGen('get_subjects',current_subject_pool);
        fix      = Fixmat(subjects,[2 4]);%all SUBJECTS, PHASES and RUNS
        %further split according to runs if wanted.
        valid    = zeros(1,length(fix.x));
        for run = runs(:)'
            valid = valid + ismember(fix.trialid , (1:120)+120*(run-1))&ismember(fix.phase,4);%run selection opeates only for phase 4
        end
        %we dont want to discard phase02 fixations
        valid    = valid + ismember(fix.phase,2);
        fix.replaceselection(valid);
        fix.ApplySelection;
        save(filename,'fix');
    else
        load(filename)
    end
    fix.kernel_fwhm = kernel_fwhm;
    varargout{1}    = fix;
elseif strcmp(varargin{1},'fix_counts')
    %% Sanity check: counts fixations in 5 different ROI
    %during baseline and generalization, returns [subjects, roi, phase].
    fixmat         = varargin{2};
    fixmat.unitize = 0;
    subjects       = unique(fixmat.subject);
    c = 0;
    for ns = subjects(:)'
        fprintf('Counting fixations in subject: %03d.\n',ns)
        c = c+1;
        p = 0;
        for phase = [2 4]
            p = p +1;
            fixmat.getmaps({'phase' phase 'subject' ns});
            dummy        = fixmat.maps;
            count(c,:,p) = fixmat.EyeNoseMouth(dummy,0);
        end
    end
    varargout{1} = count;
    
elseif  strcmp(varargin{1},'get_fixmap')
    %% General routine to compute a fixation map for a SUBJECT recorded at both phases for fixations FIX (vector) based on a FIXMAT.
    % maps are mean corrected for each phase separately.
    fixmat  = varargin{2};
    subject = varargin{3};
    fixs    = varargin{4};
    %create the query cell to talk to Fixmat object
    maps    = [];
    for phase = [2 4];
        v    = [];
        c    = 0;
        for cond = -135:45:180
            c    =  c+1;
            v{c} = {'phase', phase, 'deltacsp' cond 'subject' subject 'fix' fixs};
        end
        fixmat.getmaps(v{:});%real work done by the fixmat object.
        maps = cat(2,maps,demean(fixmat.vectorize_maps')');%within phase mean subtraction
    end
    varargout{1} = maps;
elseif  strcmp(varargin{1},'get_fixmap_oddeven')
    %% same as get_fixmat however returns 2 times many fixation maps separated by odd/even trial numbers.
    fixmat  = varargin{2};
    subject = varargin{3};
    fixs    = varargin{4};
    %creaete the query cell
    maps    = [];
    for phase = [2 4];
        v    = [];
        c    = 0;
        for cond = -135:45:180
            c    =  c+1;
            v{c} = {'phase', phase, 'deltacsp' cond 'subject' subject 'fix' fixs};
        end
        fixmat.getmaps_split(v{:});
        if ~isempty(fixmat.maps)
            maps = cat(2,maps,demean(fixmat.vectorize_maps')');
        else
            maps = [];
        end
    end
    varargout{1} = maps;
elseif strcmp(varargin{1},'plot_fdm');
    %% plot routine for FDMs used in the paper in a similar way to Figure 3A. Use the second VARARGIN to plot ROIs on top.
    % VARARGIN{1} contains fixation maps in the form of [x,y,condition].
    % The output of FPSA_FearGen('get_fixmap',...) has to be accordingly
    % reshaped. size(VARARGIN{1},3) must be a multiple of 8.
    maps          = varargin{2};
    tsubject      = size(maps,3)/8;
    contour_lines = 0;%FACIAL ROIs Plot or not.
    fs            = 18;%fontsize;
    if nargin == 3
        contour_lines = varargin{3};
    end
    %     grids         = [linspace(prctile(fixmat.maps(:),0),prctile(fixmat.maps(:),10),10) linspace(prctile(fixmat.maps(:),90),prctile(fixmat.maps(:),100),10)];
    %     [d u]         = GetColorMapLimits(maps(:),2.5);
    %     grids         = [linspace(d,u,5)];
    t             = repmat(circshift({'CS+' '+45' '+90' '+135' ['CS' char(8211)] '-135' '-90' '-45'},[1 3]),1,tsubject);
    colors        = GetFearGenColors;
    colors        = repmat(circshift(colors(1:8,:),0),tsubject,1);
    colormap jet;
    %%
    for n = 1:size(maps,3)
        if mod(n-1,8)+1 == 1;
            figure;set(gcf,'position',[1952 361 1743 714]);
        end
        hhhh(n)=subplot(1,8,mod(n-1,8)+1);
        imagesc(Fixmat([],[]).stimulus);
        hold on
        grids         = linspace(min(Vectorize(maps(:))),max(Vectorize(maps(:))),21);
        [a,h2]        = contourf(maps(:,:,n),grids,'color','none');
        caxis([grids(2) grids(end)]);
        %if n == 8
        %h4 = colorbar;
        %set(h4,'box','off','ticklength',0,'ticks',[[grids(4) grids(end-4)]],'fontsize',fs);
        %end
        hold off
        axis image;
        axis off;
        if strcmp(t{mod(n-1,8)+1}(1),'+') | strcmp(t{mod(n-1,8)+1}(1),'-')
            h= title(sprintf('%s%c',t{mod(n-1,8)+1},char(176)),'fontsize',fs,'fontweight','normal');
        else
            h= title(sprintf('%s',t{mod(n-1,8)+1}),'fontsize',fs,'fontweight','normal');
        end
        try
            h.Position = h.Position + [0 -50 0];
        end
        %
        drawnow;
        pause(.5);
        %
        try
            %%
            %             I      = find(ismember(a(1,:),h2.LevelList));
            %             [~,i]  = max(a(2,I));
            %             alphas = [repmat(.5,1,length(I))];
            %             alphas(i) = 0;
            contourf_transparency(h2,.75);
        end
        %%
        rectangle('position',[0 0 diff(xlim) diff(ylim)],'edgecolor',colors(mod(n-1,8)+1,:),'linewidth',7);
    end
    pause(1);
    for n = 1:size(maps,3);
        subplotChangeSize(hhhh(n),.01,.01);
    end
    
    if contour_lines
        hold on;
        rois = Fixmat([],[]).GetFaceROIs;
        for n = 1:4
            contour(rois(:,:,n),'k--','linewidth',1);
        end
    end
elseif strcmp(varargin{1},'plot_ROIs');
    rois = Fixmat([],[]).GetFaceROIs;
    for n = 1:4
        subplot(2,2,n)
        imagesc(Fixmat([],[]).stimulus);
        axis image;
        axis off;
        hold on;
        contour(rois(:,:,n),'k-','linewidth',1);
        hold off;
    end
    %     SaveFigure('~/Dropbox/feargen_lea/manuscript/figures/ROIs.png');
elseif strcmp(varargin{1},'get_fpsa')
    %% Routine to compute similarity matrices based on FDMs that are computed with FIXATIONS.
    % sim = FPSA_FearGen('get_fpsa',1:100) would compute a similarity matrix with all the
    % available fixations (given 1.5 presentation duration). To do it with
    % only the first fixation use 1 instead of 1:100.
    %
    % Example: sim = FPSA_FearGen('get_fpsa',1:100)
    fixations = varargin{2};
    %
    filename  = sprintf('%s/data/midlevel/fpsa_all_firstfix_%03d_lastfix_%03d_subjectpool_%03d_runs_%02d_%02d.mat',path_project,fixations(1),fixations(end),current_subject_pool,runs(1),runs(end));
    %
    if exist(filename) ==0 | force;
        fixmat   = FPSA_FearGen('get_fixmat');%returns by defaults all the 3 runs;
        subc     = 0;
        for subject = unique(fixmat.subject);
            subc                    = subc + 1;
            maps                    = FPSA_FearGen('get_fixmap',fixmat,subject,fixations);
            fprintf('Subject: %03d, Method: %s\n',subject,method);
            sim.(method)(subc,:)    = pdist(maps',method);%
        end
        save(filename,'sim');
    else
        load(filename);
    end
    varargout{1} = sim;
elseif strcmp(varargin{1},'get_fpsa2')
    %% same as get_fpsa but computes as many similarity matrices as time-windows (not used in this manuscript)
    %Time windows are computed based on WINDOW_SIZE and WINDOW_OVERLAP.
    
    force          = 0;
    %
    window_size    = varargin{2};
    window_overlap = varargin{3};
    t              = 0:1:(window_size-1);
    start_times    = 0:window_overlap:1500-window_size+1
    time           = repmat(start_times',1,length(t)) + repmat(t,length(start_times),1);
    %
    filename  = sprintf('%s/data/midlevel/fpsa2_all_windowsize_%03d_window_overlap_%03d_subjectpool_%03d_runs_%02d_%02d.mat',path_project,window_size,window_overlap,current_subject_pool,runs(1),runs(end));
    %
    if exist(filename) ==0 | force;
        tc = 0;
        for t = 1:size(time,1)-1
            tc = tc+1
            fixmat   = FPSA_FearGen('get_fixmat');
            fixmat.UpdateSelection('start',time(t,:),'stop',time(t+1,:));
            fixmat.ApplySelection;
            subc     = 0;
            for subject = unique(fixmat.subject);
                subc                    = subc + 1;
                maps                    = FPSA_FearGen('get_fixmap',fixmat,subject,1:100);
                fprintf('Subject: %03d, Method: %s\n',subject,method);
                sim.(method)(subc,:,tc)    = pdist(maps',method);%
            end
            try
                figure(1);imagesc(squareform(nanmean(sim.correlation(:,:,tc))));drawnow;
            end
        end
        save(filename,'sim');
    else
        load(filename);
    end
    varargout{1} = sim;
    
elseif strcmp(varargin{1},'get_fpsa_fair')
    %% gets an FPSA matrix per run to be fair to the baseline condition (main routine to get similarity matrices).
    % the FPSA for the 3 test-phase runs are individually computed and averaged.
    % Doing it the other way (i.e. average FDMs from the 3 phases and compute
    % FPSA as in get_fpsa) would have led to comparably less noisy FDMs for the test
    % phase and thus differences btw B and T simply because the number of
    % trials are different. See (4) for more information on how noise
    % affects similarity values
    %
    % Example usage:
    % sim = FPSA_FearGen('get_fpsa_fair',1:100,1:3);
    
    fixations = varargin{2};%which fixations
    runs      = varargin{3};%whichs runs would you like to have
    %
    for run = runs
        filename     = sprintf('%s/data/midlevel/fpsa_fair_firstfix_%03d_lastfix_%03d_subjectpool_%03d_run_%02d.mat',path_project,fixations(1),fixations(end),current_subject_pool,run(1));
        if exist(filename) ==0 | force;
            fixmat   = FPSA_FearGen('get_fixmat','runs',run);
            subc     = 0;
            for subject = unique(fixmat.subject);
                subc                    = subc + 1;
                maps                    = FPSA_FearGen('get_fixmap',fixmat,subject,fixations);
                fprintf('Subject: %03d, Run: %03d, Method: %s\n',subject,run,method);
                sim.(method)(subc,:,run)= pdist(maps',method);%
            end
            %average across runs
            sim.(method) = mean(sim.(method),3);
            save(filename,'sim');
        else
            load(filename);
        end
    end
    varargout{1} = sim;
    
elseif strcmp(varargin{1},'get_fpsa_oddeven')
    %% Computes FPSA with cross-validation based on odd and even trials.
    %  See also: get_fixmap_oddeven
    
    fixations = varargin{2};
    filename  = sprintf('%s/data/midlevel/fpsa_all_oddeven_firstfix_%03d_lastfix_%03d_subjectpool_%03d_runs_%02d_%02d.mat',path_project,fixations(1),fixations(end),current_subject_pool,runs(1),runs(end));
    force     = 0;
    %
    if exist(filename) ==0 | force
        fixmat   = FPSA_FearGen('get_fixmat');
        subc     = 0;
        for subject = setdiff(unique(fixmat.subject),58);
            subc                    = subc + 1;
            maps                    = FPSA_FearGen('get_fixmap_oddeven',fixmat,subject,fixations);
            if size(maps,2) == 32;
                fprintf('Subject: %03d, Method: %s\n',subject,method);
                sim.(method)(subc,:)    = pdist(maps',method);%
            else
                fprintf('no map found...\n');
                subc = subc - 1;
            end
        end
        save(filename,'sim');
    else
        load(filename);
    end
    varargout{1} = sim;
elseif strcmp(varargin{1},'plot_fpsa');
    %% A routine to plot similarity matrices
    figure;
    sim     = varargin{2};
    cormatz = squareform(nanmean(sim.correlation));
    cormatz = CancelDiagonals(cormatz,NaN);
    [d u]   = GetColorMapLimits(cormatz,2.5);
    imagesc(cormatz,[d u]);
    axis square;colorbar
    set(gca,'fontsize',15);
    axis off;
    
elseif strcmp(varargin{1},'get_block')
    %% will get the Yth, Xth block from similarity matrix SIM.
    % SQFM is the square_form of SIM.
    %
    % Example: fpsa = FPSA_FearGen('get_block',FPSA_FearGen('get_fpsa',1:100),2,2)
    sim  = varargin{2};
    y    = varargin{3};
    x    = varargin{4};
    r    = [];
    sqfm = [];
    for ns = 1:size(sim.correlation,1)
        dummy = squareform(sim.correlation(ns,:));
        B     = block_extract(dummy,y,x,1);%funhandle defined at the top.
        r     = cat(3,r,B);
        sqfm  = [sqfm;squareform(B)];
    end
    varargout{1} = r;
    varargout{2} = sqfm;
elseif strcmp(varargin{1},'get_mdscale')
    %% Routine to make MDS analysis using a SIMilarity matrix with NDIMENsions.
    %
    % Example: FPSA_FearGen('get_mdscale',mean(sim.correlation),2);
    sim                         = varargin{2};%sim is a valid similarity matrix;
    ndimen                      = varargin{3};
    viz                         = 1;
    [dummy stress disparities]  = mdscale(sim,ndimen,'Criterion',criterion,'start','cmdscale','options',statset('display','final','tolfun',10^-12,'tolx',10^-12));
    dummy                       = dummy';
    Y                           = dummy(:);
    varargout{1}                = Y;
    if viz
        FPSA_FearGen('plot_mdscale',Y);
    end
elseif strcmp(varargin{1},'plot_mdscale')
    %% Routine to plot the results of the get_mdscale
    Y      = varargin{2};
    ndimen = length(Y)./16;
    y      = reshape(Y,length(Y)/16,16)';%to make it easy plotting put coordinates to different columns;
    colors = GetFearGenColors;
    colors = [colors(1:8,:);colors(1:8,:)];
    if ndimen == 2
        plot(y([1:8 1],1),y([1:8 1],2),'.-.','linewidth',3,'color',[.6 .6 .6]);
        hold on;
        plot(y([1:8 1]+8,1),y([1:8 1]+8,2),'k.-.','linewidth',3);
        for nface = 1:16
            plot(y(nface,1),y(nface,2),'.','color',colors(nface,:),'markersize',120,'markerface',colors(nface,:));
        end
        hold off;
        %
        for n = 1:16;text(y(n,1),y(n,2),mat2str(mod(n-1,8)+1),'fontsize',25);end
    elseif ndimen == 3
        plot3(y([1:8 1],1),y([1:8 1],2),y([1:8 1],3),'o-','linewidth',3);
        hold on;
        plot3(y([1:8 1]+8,1),y([1:8 1]+8,2),y([1:8 1]+8,3),'ro-','linewidth',3);
        hold off;
        for n = 1:16;text(y(n,1),y(n,2),y(n,3),mat2str(mod(n-1,8)+1),'fontsize',25);end
    end
    
    
elseif strcmp(varargin{1},'FPSA_get_table')
    %% returns a table object for the FPSA modelling with fitlm, fitglm, etc.
    %
    % the table object contains the following variable names:
    % FPSA_B      : similarity matrices from baseline.
    % FPSA_G      : similarity matrices from test.
    % circle     : circular predictor consisting of a sum of specific and
    % unspecific components.
    % specific   : specific component based on quadrature decomposition
    % (the cosine factor).
    % unspecific : unspecific component based on quadrature decomposition
    % (the sine factor).
    % Gaussian : generalization of the univariate Gaussian component to
    % the similarity space.
    % subject    : indicator variable for subjects
    % phase      : indicator variable for baseline and generalizaation
    % phases.
    %
    % Example: FPSA_FearGen('FPSA_get_table',1:100)
    fixations = varargin{2};
    force = 1;
    runs      = 1:3;
    filename  = sprintf('%s/data/midlevel/fpsa_modelling_table_firstfix_%03d_lastfix_%03d_subjectpool_%03d_runs_%02d_%02d.mat',path_project,fixations(1),fixations(end),current_subject_pool,runs(1),runs(end));
    if ~exist(filename) | force
        %the full B and T similarity matrix which are jointly computed;
        sim       = FPSA_FearGen('get_fpsa_fair',fixations,runs);%
        %%we only want the B and T parts
        B         = FPSA_FearGen('get_block',sim,1,1);
        T         = FPSA_FearGen('get_block',sim,2,2);
        %once we have these, we go back to the compact form and concat the
        %stuff, now each column is a non-redundant FPSA per subject
        for n = 1:size(sim.correlation,1)
            BB(n,:) = squareform(B(:,:,n));
            TT(n,:) = squareform(T(:,:,n));
        end
        BB       = BB';
        TT       = TT';
        % some indicator variables for phase, subject identities.
        phase    = repmat([repmat(1,size(BB,1)/2,1); repmat(2,size(BB,1)/2,1)],1,size(BB,2));
        subject  = repmat(1:size(sim.correlation,1),size(BB,1),1);
        S        = subject(:);
        P        = phase(:);
        %% our models:
        %MODEL1: perfectly circular similarity model;
        %MODEL2: flexible circular similarity model;
        %MODEL3: Model2 + a Gaussian.
        % a circular FPSA matrix for B and T replicated by the number of subjects
        x          = [pi/4:pi/4:2*pi];
        w          = [cos(x);sin(x)];
        model1     = repmat(repmat(squareform_force(w'*w),1,1),1,size(subject,2));%we use squareform_force as the w'*w is not perfectly positive definite matrix due to rounding errors.
        %
        model2_c   = repmat(repmat(squareform_force(cos(x)'*cos(x)),1,1),1,size(subject,2));%
        model2_s   = repmat(repmat(squareform_force(sin(x)'*sin(x)),1,1),1,size(subject,2));%
        %
        %getcorrmat(amp_circ, amp_gau, amp_const, amp_diag, varargin)
        [cmat]     = getcorrmat(0,3,1,1);%see model_rsa_testgaussian_optimizer
        model3_g   = repmat(repmat(squareform_force(cmat),1,1),1,size(subject,2));%
        %% add all this to a TABLE object.
        t          = table(1-BB(:),1-TT(:),model1(:),model2_c(:),model2_s(:),model3_g(:),categorical(subject(:)),categorical(phase(:)),'variablenames',{'FPSA_B' 'FPSA_G' 'circle' 'specific' 'unspecific' 'Gaussian' 'subject' 'phase'});
        save(filename,'t');
    else
        load(filename);
    end
    varargout{1} = t;
elseif strcmp(varargin{1},'FPSA_model');
    %%
    
    fixations  = varargin{2};
    t          = FPSA_FearGen('FPSA_get_table',fixations);
    %% MIXED EFFECT MODEL
    % null model
    out.baseline.model_00_mixed          = fitlme(t,'FPSA_B ~ 1 + (1|subject)');
    out.generalization.model_00_mixed    = fitlme(t,'FPSA_G ~ 1 + (1|subject)');
    % FPSA_model_bottom-up model
    out.baseline.model_01_mixed          = fitlme(t,'FPSA_B ~ 1 + circle + (1 + circle|subject)');
    out.generalization.model_01_mixed    = fitlme(t,'FPSA_G ~ 1 + circle + (1 + circle|subject)');
    % FPSA_model_adversitycateg
    out.baseline.model_02_mixed          = fitlme(t,'FPSA_B ~ 1 + specific + unspecific +  (1 + specific + unspecific|subject)');
    out.generalization.model_02_mixed    = fitlme(t,'FPSA_G ~ 1 + specific + unspecific +  (1 + specific + unspecific|subject)');
    % FPSA_model_adversitytuning
    out.baseline.model_03_mixed          = fitlme(t,'FPSA_B ~ 1 + specific + unspecific + Gaussian + (1 + specific + unspecific + Gaussian|subject)');
    out.generalization.model_03_mixed    = fitlme(t,'FPSA_G ~ 1 + specific + unspecific + Gaussian + (1 + specific + unspecific + Gaussian|subject)');
    
    %% FIXED EFFECT MODEL
    % FPSA null model
    out.baseline.model_00_fixed          = fitlm(t,'FPSA_B ~ 1');
    out.generalization.model_00_fixed    = fitlm(t,'FPSA_G ~ 1');
    % FPSA_model_bottom-up model
    out.baseline.model_01_fixed          = fitlm(t,'FPSA_B ~ 1 + circle');
    out.generalization.model_01_fixed    = fitlm(t,'FPSA_G ~ 1 + circle');
    % FPSA_model_adversitycateg
    out.baseline.model_02_fixed          = fitlm(t,'FPSA_B ~ 1 + specific + unspecific');
    out.generalization.model_02_fixed    = fitlm(t,'FPSA_G ~ 1 + specific + unspecific');
    % FPSA_model_adversitytuning
    out.baseline.model_03_fixed          = fitlm(t,'FPSA_B ~ 1 + specific + unspecific + Gaussian');
    out.generalization.model_03_fixed    = fitlm(t,'FPSA_G ~ 1 + specific + unspecific + Gaussian');
    varargout{1}   = out;
    
elseif strcmp(varargin{1},'model2text')
    %handy function to dump model output to a text file that let you easily
    %paste it to the manuscript ;-\
    
    model = varargin{2};
    a     = evalc('disp(model)');
    fid   = fopen(sprintf('%s/data/midlevel/%s.txt',path_project,model.Formula),'w');
    fwrite(fid,a);
    fclose(fid);
    
    
elseif strcmp(varargin{1},'FPSA_model_singlesubject');
    %% Fits FPSA matrices decribed models.
    fixations  = varargin{2};
    t          = FPSA_FearGen('FPSA_get_table',fixations);
    %% test the model for B, T
    
    Model.model_01.w1 = [];
    Model.model_02.w1 = [];
    Model.model_02.w2 = [];
    Model.model_03.w1 = [];
    Model.model_03.w2 = [];
    Model.model_03.w3 = [];
    for ns = unique(t.subject)'
        fprintf('Fitting an circular and flexibile LM to subject %03d...\n',double(ns));
        t2                = t(ismember(t.subject,categorical(ns)),:);
        B                 = fitlm(t2,'FPSA_B ~ 1 + circle');
        T                 = fitlm(t2,'FPSA_G ~ 1 + circle');
        Model.model_01.w1 = [Model.model_01.w1; [B.Coefficients.Estimate(2) T.Coefficients.Estimate(2)]];
        %
        B                 = fitlm(t2,'FPSA_B ~ 1 + specific + unspecific');
        T                 = fitlm(t2,'FPSA_G ~ 1 + specific + unspecific');
        Model.model_02.w1 = [Model.model_02.w1; [B.Coefficients.Estimate(2) T.Coefficients.Estimate(2)]];
        Model.model_02.w2 = [Model.model_02.w2; [B.Coefficients.Estimate(3) T.Coefficients.Estimate(3)]];
        %
        B                 = fitlm(t2,'FPSA_B ~ 1 + specific + unspecific + Gaussian');
        T                 = fitlm(t2,'FPSA_G ~ 1 + specific + unspecific + Gaussian');
        Model.model_03.w1 = [Model.model_03.w1; [B.Coefficients.Estimate(2) T.Coefficients.Estimate(2)]];
        Model.model_03.w2 = [Model.model_03.w2; [B.Coefficients.Estimate(3) T.Coefficients.Estimate(3)]];
        Model.model_03.w3 = [Model.model_03.w3; [B.Coefficients.Estimate(4) T.Coefficients.Estimate(4)]];
    end
    varargout{1} = Model;
    
elseif strcmp(varargin{1},'figure_04C');
    %% plots the main model comparison figure;
    fixations = varargin{2};
    C         = FPSA_FearGen('FPSA_model_singlesubject',fixations);
    %%
    % circular model
    M         = mean(C.model_01.w1);
    SEM       = std(C.model_01.w1)./sqrt(61);
    %flexible model
    Mc        = mean(C.model_02.w1);
    SEMc      = std(C.model_02.w1)./sqrt(61);
    Ms        = mean(C.model_02.w2);
    SEMs      = std(C.model_02.w2)./sqrt(61);
    %gaussian model
    Mcg       = mean(C.model_03.w1);
    SEMcg     = std(C.model_03.w1)./sqrt(61);
    Msg       = mean(C.model_03.w2);
    SEMsg     = std(C.model_03.w2)./sqrt(61);
    Mg        = mean(C.model_03.w3);
    SEMg      = std(C.model_03.w3)./sqrt(61);
    
    
    %% get the p-values
    [H   P]     = ttest(C.model_01.w1(:,1)-C.model_01.w1(:,2));%compares baseline vs test the circular model parameters
    
    [Hc Pc]     = ttest(C.model_02.w1(:,1)-C.model_02.w1(:,2));%compares cosine before to after
    [Hs Ps]     = ttest(C.model_02.w2(:,1)-C.model_02.w2(:,2));%compares sine before to after
    [Hcs Pcs]   = ttest(C.model_02.w1(:,2)-C.model_02.w2(:,2));%compares cosine after to sine after
    % same as before
    [Hgc Pgc]   = ttest(C.model_03.w1(:,1)-C.model_03.w1(:,2));%compares cosine before to after
    [Hgcs Pgcs] = ttest(C.model_03.w1(:,2)-C.model_03.w2(:,2));%compares cosine after to sine after
    [Hgg Pgg]   = ttest(C.model_03.w3(:,1)-C.model_03.w3(:,2));%compares sine before to after
    
    %%
    %%
    figure(fixations(1));
    if ispc
        set(gcf,'position',[-200+500*fixations(1)        1400        898         604]);
    else
        set(gcf,'position',[2150         335         898         604]);
    end
    X    = [1 2 4 5 6 7  9 10 11 12 13 14]/1.5;
    Y    = [M Mc Ms Mcg Msg Mg];
    Y2   = [SEM SEMc SEMs SEMcg SEMsg SEMg];
    bw   = .5;
    hold off;
    for n = 1:size(Y,2)
        h       = bar(X(n),Y(n));
        legh(n) = h;
        hold on
        try %capsize is 2016b compatible.
            errorbar(X(n),Y(n),Y2(n),'k.','marker','none','linewidth',1.5,'capsize',10);
        catch
            errorbar(X(n),Y(n),Y2(n),'k.','marker','none','linewidth',1.5);
        end
        if ismember(n,[1 3 5 7 9 11])
            try %2016b compatibility.
                set(h,'FaceAlpha',.1,'FaceColor','w','EdgeAlpha',1,'EdgeColor',[0 0 0],'LineWidth',1.5,'BarWidth',bw,'LineStyle','-');
            catch
                set(h,'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',1.5,'BarWidth',bw,'LineStyle','-');
            end
        else
            try
                set(h,'FaceAlpha',.5,'FaceColor',[0 0 0],'EdgeAlpha',0,'EdgeColor',[.4 .4 .4],'LineWidth',1,'BarWidth',bw,'LineStyle','-');
            catch
                set(h,'FaceColor',[0 0 0],'EdgeColor',[.4 .4 .4],'LineWidth',1,'BarWidth',bw,'LineStyle','-');
            end
        end
    end
    box off;
    L          = legend(legh(1:2),{'Baseline' 'Generaliz.'},'box','off');
    try
        L.Position = L.Position + [0.1/2 0 0 0];
        L.FontSize = 12;
    end
    set(gca,'linewidth',1.8);
    % xticks
    xtick = [mean(X(1:2)) mean(X(3:4)) mean(X(5:6)) mean(X(7:8)) mean(X(9:10)) mean(X(11:12)) ];
    label = {'\itw_{\rmcircle}' '\itw_{\rmspec.}' '\itw_{\rmunspec.}' '\itw_{\rmspec.}' '\itw_{\rmunspec.}' '\itw_{\rmGaus.}' };
    for nt = 1:length(xtick)
        h = text(xtick(nt),-.02,label{nt},'horizontalalignment','center','fontsize',20,'rotation',45,'fontangle','italic','fontname','times new roman');
    end
    try
        set(gca,'xtick',[3 8]./1.5,'xcolor','none','color','none','XGrid','on','fontsize',16);
    catch
        set(gca,'xtick',[3 8]./1.5,'color','none','XGrid','on','fontsize',16);
    end
    %
    text(-.5,.215,'\beta','fontsize',28,'fontweight','bold');
    
    
    ylim([-.05 .2]);
    set(gca,'ytick',-.05:.05:.2,'yticklabel',{'-.05' '0' '.05' '.1' '.15' '.2'})
    axis normal
    % asteriks
    hold on
    ylim([min(ylim) .16]);
    h= line([X(1)-bw/2 X(2)+bw/2],repmat(max(ylim),1,2)-.01);set(h,'color','k','linewidth',1);
    h= line([X(3)-bw/2 X(4)+bw/2],repmat(max(ylim),1,2)-.01);set(h,'color','k','linewidth',1);
    h= line([X(4)-bw/2 X(6)+bw/2],repmat(max(ylim),1,2)-.0025);set(h,'color','k','linewidth',1);
    
    h= line([X(7)-bw/2 X(8)+bw/2],repmat(max(ylim),1,2)-.01);set(h,'color','k','linewidth',1);
%     h= line([X(8)-bw/2 X(10)+bw/2],repmat(max(ylim),1,2)-.0025);set(h,'color','k','linewidth',1);
%     h= line([X(11)-bw/2 X(12)+bw/2],repmat(max(ylim),1,2)-.1);set(h,'color','k','linewidth',1);
    %
    text(mean(X(1:2))  ,max(ylim)-.0075, pval2asterix(P),'HorizontalAlignment','center','fontsize',16);
    text(mean(X(3:4))  ,max(ylim)-.0075, pval2asterix(Pc),'HorizontalAlignment','center','fontsize',16);
    text(mean(X([4 6])),max(ylim)      , pval2asterix(Pcs),'HorizontalAlignment','center','fontsize',16);
    text(mean(X([7 8])),max(ylim)-.0075, pval2asterix(Pgc),'HorizontalAlignment','center','fontsize',16);
%     text(mean(X([8 10])),max(ylim)      , pval2asterix(Pgcs),'HorizontalAlignment','center','fontsize',16);
%     text(mean(X([11 12])),max(ylim)-.09      , pval2asterix(Pgg),'HorizontalAlignment','center','fontsize',12);
    % model names
    ylim([-.04 .2])
    %     h = line([X(1)-bw/2 X(2)+bw/2],[-.022 -.022],'linestyle','--');
    %     set(h(1),'color','k','linewidth',1,'clipping','off');
    text(mean(X(1:2)),.18,sprintf('Arousal\nmodel'),'Rotation',0,'HorizontalAlignment','center','FontWeight','normal','fontname','Helvetica','fontsize',14,'verticalalignment','bottom');
    %     h = line([X(3)-bw/2 X(6)+bw/2],[-.022 -.022],'linestyle','--');
    %     set(h(1),'color','k','linewidth',1,'clipping','off');
    text(mean(X(3:6)),.18,sprintf('Adversity\nCategorization\nmodel'),'Rotation',0,'HorizontalAlignment','center','FontWeight','normal','fontname','Helvetica','fontsize',14,'verticalalignment','bottom');
    
    %     h = line([X(7)-bw/2 X(end)+bw/2],[-.022 -.022],'linestyle','--');
    %     set(h(1),'color','k','linewidth',1,'clipping','off');
    text(mean(X(7:end)),.18,sprintf('Adversity\nTuning\nmodel'),'Rotation',0,'HorizontalAlignment','center','FontWeight','normal','fontname','Helvetica','fontsize',14,'verticalalignment','bottom');
    %%
    %     SaveFigure('~/Dropbox/feargen_lea/manuscript/figures/figure_03E.png');
    %
    
elseif strcmp(varargin{1},'figure_04C_BDNFcheck');
    %% plots the main model comparison figure;
    fixations = varargin{2};
    C         = FPSA_FearGen('FPSA_model_singlesubject',fixations);
    %%
    subgr{1} = ones(61,1);
    subgr{2} = Project.BDNF(ismember(Project.subjects_bdnf,FPSA_FearGen('get_subjects')))==1;
    subgr{3} = Project.BDNF(ismember(Project.subjects_bdnf,FPSA_FearGen('get_subjects')))==2;
    
    for n = 1:3
        inds = logical(subgr{n});
        M(n,:)         = mean(C.model_01.w1(inds,:));
        SEM(n,:)       = std(C.model_01.w1(inds,:))./sqrt(sum(subgr{n}));
        %flexible model
        Mc(n,:)        = mean(C.model_02.w1(inds,:));
        SEMc(n,:)      = std(C.model_02.w1(inds,:))./sqrt(sum(subgr{n}));
        Ms(n,:)        = mean(C.model_02.w2(inds,:));
        SEMs(n,:)      = std(C.model_02.w2(inds,:))./sqrt(sum(subgr{n}));
        %gaussian model
        Mcg(n,:)       = mean(C.model_03.w1(inds,:));
        SEMcg(n,:)     = std(C.model_03.w1(inds,:))./sqrt(sum(subgr{n}));
        Msg(n,:)       = mean(C.model_03.w2(inds,:));
        SEMsg(n,:)     = std(C.model_03.w2(inds,:))./sqrt(sum(subgr{n}));
        Mg(n,:)        = mean(C.model_03.w3(inds,:));
        SEMg(n,:)      = std(C.model_03.w3(inds,:))./sqrt(sum(subgr{n}));
        
        %% get the p-values
        [H(n)   P(n)]     = ttest(C.model_01.w1(inds,1)-C.model_01.w1(inds,2));%compares baseline vs test the circular model parameters
        
        [Hc(n) Pc(n)]     = ttest(C.model_02.w1(inds,1)-C.model_02.w1(inds,2));%compares cosine before to after
        [Hs(n) Ps(n)]     = ttest(C.model_02.w2(inds,1)-C.model_02.w2(inds,2));%compares sine before to after
        [Hcs(n) Pcs(n)]   = ttest(C.model_02.w1(inds,2)-C.model_02.w2(inds,2));%compares cosine after to sine after
        % same as before
        [Hgc(n) Pgc(n)]   = ttest(C.model_03.w1(inds,1)-C.model_03.w1(inds,2));%compares cosine before to after
        [Hgcs(n) Pgcs(n)] = ttest(C.model_03.w1(inds,2)-C.model_03.w2(inds,2));%compares cosine after to sine after
        [Hgg(n) Pgg(n)]   = ttest(C.model_03.w3(inds,1)-C.model_03.w3(inds,2));%compares sine before to after
    end
    %%
    %%
    figure;
    set(gcf,'position',[0        0        898         604]);
    X    = [1 2 4 5 6 7  9 10 11 12 13 14]/1.5;
    Y    = [M Mc Ms Mcg Msg Mg];
    Y2   = [SEM SEMc SEMs SEMcg SEMsg SEMg];
    bw   = .5;
    hold off;
    groupcolors = {[.4 .4 .4],[0 0 1],[1 0 0]};
    for subpl = 1:3
        subplot(3,1,subpl);
        for n = 1:size(Y,2)
            h       = bar(X(n),Y(subpl,n));
            legh(n) = h;
            hold on
            try %capsize is 2016b compatible.
                errorbar(X(n),Y(subpl,n),Y2(subpl,n),'k.','marker','none','linewidth',1.5,'capsize',10);
            catch
                errorbar(X(n),Y(subpl,n),Y2(subpl,n),'k.','marker','none','linewidth',1.5);
            end
            if ismember(n,[1 3 5 7 9 11])
                try %2016b compatibility.
                    set(h,'FaceAlpha',.1,'FaceColor','w','EdgeAlpha',1,'EdgeColor',[0 0 0],'LineWidth',1.5,'BarWidth',bw,'LineStyle','-');
                catch
                    set(h,'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',1.5,'BarWidth',bw,'LineStyle','-');
                end
            else
                try
                    set(h,'FaceAlpha',.5,'FaceColor',[0 0 0],'EdgeAlpha',0,'EdgeColor',groupcolors{subpl},'LineWidth',1,'BarWidth',bw,'LineStyle','-');
                catch
                    set(h,'FaceColor',[0 0 0],'EdgeColor',groupcolors{subpl},'LineWidth',1,'BarWidth',bw,'LineStyle','-');
                end
            end
        end
        box off;
        L          = legend(legh(1:2),{'Baseline' 'Generaliz.'},'box','off');
        try
            L.Position = L.Position + [0.1/2 0 0 0];
            L.FontSize = 12;
        end
        set(gca,'linewidth',1.8);
        % xticks
        xtick = [mean(X(1:2)) mean(X(3:4)) mean(X(5:6)) mean(X(7:8)) mean(X(9:10)) mean(X(11:12)) ];
        label = {'\itw_{\rmcircle}' '\itw_{\rmspec.}' '\itw_{\rmunspec.}' '\itw_{\rmspec.}' '\itw_{\rmunspec.}' '\itw_{\rmGaus.}' };
        for nt = 1:length(xtick)
            h = text(xtick(nt),-.02,label{nt},'horizontalalignment','center','fontsize',20,'rotation',45,'fontangle','italic','fontname','times new roman');
        end
        try
            set(gca,'xtick',[3 8]./1.5,'xcolor','none','color','none','XGrid','on','fontsize',16);
        catch
            set(gca,'xtick',[3 8]./1.5,'color','none','XGrid','on','fontsize',16);
        end
        %
        text(-.5,.215,'\beta','fontsize',28,'fontweight','bold');
        
        
        ylim([-.05 .2]);
        set(gca,'ytick',-.05:.05:.2,'yticklabel',{'-.05' '0' '.05' '.1' '.15' '.2'})
        axis normal
        % asteriks
        hold on
        ylim([min(ylim) .16]);
        h= line([X(1)-bw/2 X(2)+bw/2],repmat(max(ylim),1,2)-.01);set(h,'color','k','linewidth',1);
        h= line([X(3)-bw/2 X(4)+bw/2],repmat(max(ylim),1,2)-.01);set(h,'color','k','linewidth',1);
        h= line([X(4)-bw/2 X(6)+bw/2],repmat(max(ylim),1,2)-.0025);set(h,'color','k','linewidth',1);
        
        h= line([X(7)-bw/2 X(8)+bw/2],repmat(max(ylim),1,2)-.01);set(h,'color','k','linewidth',1);
        h= line([X(8)-bw/2 X(10)+bw/2],repmat(max(ylim),1,2)-.0025);set(h,'color','k','linewidth',1);
        h= line([X(11)-bw/2 X(12)+bw/2],repmat(max(ylim),1,2)-.1);set(h,'color','k','linewidth',1);
        %
        text(mean(X(1:2))  ,max(ylim)-.0075, pval2asterix(P(subpl)),'HorizontalAlignment','center','fontsize',12);
        text(mean(X(3:4))  ,max(ylim)-.0075, pval2asterix(Pc(subpl)),'HorizontalAlignment','center','fontsize',12);
        text(mean(X([4 6])),max(ylim)      , pval2asterix(Pcs(subpl)),'HorizontalAlignment','center','fontsize',12);
        text(mean(X([7 8])),max(ylim)-.0075, pval2asterix(Pgc(subpl)),'HorizontalAlignment','center','fontsize',12);
        text(mean(X([8 10])),max(ylim)      , pval2asterix(Pgcs(subpl)),'HorizontalAlignment','center','fontsize',12);
        text(mean(X([11 12])),max(ylim)-.09      , pval2asterix(Pgg(subpl)),'HorizontalAlignment','center','fontsize',12);
        % model names
        ylim([-.04 .2])
        %     h = line([X(1)-bw/2 X(2)+bw/2],[-.022 -.022],'linestyle','--');
        %     set(h(1),'color','k','linewidth',1,'clipping','off');
        text(mean(X(1:2)),.18,sprintf('Arousal\nmodel'),'Rotation',0,'HorizontalAlignment','center','FontWeight','normal','fontname','Helvetica','fontsize',14,'verticalalignment','bottom');
        %     h = line([X(3)-bw/2 X(6)+bw/2],[-.022 -.022],'linestyle','--');
        %     set(h(1),'color','k','linewidth',1,'clipping','off');
        text(mean(X(3:6)),.18,sprintf('Adversity\nCategorization\nmodel'),'Rotation',0,'HorizontalAlignment','center','FontWeight','normal','fontname','Helvetica','fontsize',14,'verticalalignment','bottom');
        
        %     h = line([X(7)-bw/2 X(end)+bw/2],[-.022 -.022],'linestyle','--');
        %     set(h(1),'color','k','linewidth',1,'clipping','off');
        text(mean(X(7:end)),.18,sprintf('Univariate\nGeneralization\nmodel'),'Rotation',0,'HorizontalAlignment','center','FontWeight','normal','fontname','Helvetica','fontsize',14,'verticalalignment','bottom');
        %%
        %     SaveFigure('~/Dropbox/feargen_lea/manuscript/figures/figure_03E.png');
        %
    end
    
    
elseif strcmp(varargin{1},'model2behavior')
    C         = FPSA_FearGen('FPSA_model_singlesubject',1:100);
    subs      = FPSA_FearGen('get_subjects'); %be sure to make current_subject_pool =0
    %
    %normalize ellipsoid parameter by other parameter.
    beta1 = C.model_02.w1(:,2);
    beta2 = C.model_02.w2(:,2);
    param = (beta1-beta2)./(beta1+beta2); %% ellipsoidness is specific vs unspecific beta
    
    %if baseline corrected:
    beta1bc = C.model_02.w1(:,2)-C.model_02.w1(:,1);
    beta2bc = C.model_02.w2(:,2)-C.model_02.w2(:,1);
    parambc = (beta1bc-beta2bc)./(beta1bc+beta2bc);
    
    
    % collect rating amplitudes
    ns = 0;
    for sub = subs(:)'
        ns = ns+1;
        s = Subject(sub);
        amp_cond(ns) = s.get_fit('rating',3).params(1);
        amp_test(ns) = s.get_fit('rating',4).params(1);
    end
    
    % get SCR ampl CS+ > CS-
    g   = Group(subs);
    out = g.getSCR(2.5:5.5);
    scr_cond = out.y(:,13)-out.y(:,17);
    scr_test = out.y(:,22)-out.y(:,26);
    
    selector = ismember(subs,g.ids);
    %%
    figure(1);
    subplot(2,2,1);plot(param(selector),scr_cond,'bo');lsline;ylabel('SCR CS+>CS- Cond'); xlabel('(b1-b2)/sum(b)');box off
    subplot(2,2,2);plot(param(selector),scr_test,'bo');lsline;ylabel('SCR CS+>CS- Test'); xlabel('(b1-b2)/sum(b)');box off
    subplot(2,2,3);plot(param,amp_cond,'ro');lsline;ylabel('Ampl. Rating Cond'); xlabel('(b1-b2)/sum(b)');box off
    subplot(2,2,4);plot(param,amp_test,'ro');lsline;ylabel('Ampl. Rating Test'); xlabel('(b1-b2)/sum(b)');box off
    st=supertitle('modelparam vs behavior');set(st,'fontsize',16,'position',[0 .5]);
    figure(2);
    subplot(2,2,1);plot(parambc(selector),scr_cond,'bo');lsline;ylabel('SCR CS+>CS- Cond'); xlabel('(b1-b2)/sum(b)');box off
    subplot(2,2,2);plot(parambc(selector),scr_test,'bo');lsline;ylabel('SCR CS+>CS- Test'); xlabel('(b1-b2)/sum(b)');box off
    subplot(2,2,3);plot(parambc,amp_cond,'ro');lsline;ylabel('Ampl. Rating Cond'); xlabel('(b1-b2)/sum(b)');box off
    subplot(2,2,4);plot(parambc,amp_test,'ro');lsline;ylabel('Ampl. Rating Test'); xlabel('(b1-b2)/sum(b)');box off
    st=supertitle('modelparam vs behavior, baseline corrected');set(st,'fontsize',16,'position',[0 .5]);
    figure(3);
    subplot(2,2,1);plot(beta1(selector),scr_cond,'bo');lsline;ylabel('SCR CS+>CS- Cond'); xlabel('beta1');box off
    subplot(2,2,2);plot(beta1(selector),scr_test,'bo');lsline;ylabel('SCR CS+>CS- Test'); xlabel('beta1');box off
    subplot(2,2,3);plot(beta1,amp_cond,'ro');lsline;ylabel('Ampl. Rating Cond'); xlabel('beta1');box off
    subplot(2,2,4);plot(beta1,amp_test,'ro');lsline;ylabel('Ampl. Rating Test'); xlabel('beta1');box off
    st=supertitle('beta1 vs behavior');set(st,'fontsize',16,'position',[0 .5]);
    figure(4);
    subplot(2,2,1);plot(beta1bc(selector),scr_cond,'bo');lsline;ylabel('SCR CS+>CS- Cond'); xlabel('beta1bc');box off
    subplot(2,2,2);plot(beta1bc(selector),scr_test,'bo');lsline;ylabel('SCR CS+>CS- Test'); xlabel('beta1bc');box off
    subplot(2,2,3);plot(beta1bc,amp_cond,'ro');lsline;ylabel('Ampl. Rating Cond'); xlabel('beta1bc');box off
    subplot(2,2,4);plot(beta1bc,amp_test,'ro');lsline;ylabel('Ampl. Rating Test'); xlabel('beta1bc');box off
    st=supertitle('beta1 vs behavior, baseline corrected');set(st,'fontsize',16,'position',[0 .5]);
    figure(5);
    subplot(2,2,1);plot(beta2(selector),scr_cond,'bo');lsline;ylabel('SCR CS+>CS- Cond'); xlabel('beta2');box off
    subplot(2,2,2);plot(beta2(selector),scr_test,'bo');lsline;ylabel('SCR CS+>CS- Test'); xlabel('beta2');box off
    subplot(2,2,3);plot(beta2,amp_cond,'ro');lsline;ylabel('Ampl. Rating Cond'); xlabel('beta2');box off
    subplot(2,2,4);plot(beta2,amp_test,'ro');lsline;ylabel('Ampl. Rating Test'); xlabel('beta2');box off
    st=supertitle('beta2 vs behavior');set(st,'fontsize',16,'position',[0 .5]);
    figure(6);
    subplot(2,2,1);plot(beta2bc(selector),scr_cond,'bo');lsline;ylabel('SCR CS+>CS- Cond'); xlabel('beta2bc');box off
    subplot(2,2,2);plot(beta2bc(selector),scr_test,'bo');lsline;ylabel('SCR CS+>CS- Test'); xlabel('beta2bc');box off
    subplot(2,2,3);plot(beta2bc,amp_cond,'ro');lsline;ylabel('Ampl. Rating Cond'); xlabel('beta2bc');box off
    subplot(2,2,4);plot(beta2bc,amp_test,'ro');lsline;ylabel('Ampl. Rating Test'); xlabel('beta2bc');box off
    st=supertitle('beta2 vs behavior, baseline corrected');set(st,'fontsize',16,'position',[0 .5]);
    
    names = {'beta1-2' 'beta1-2BC' 'beta1','beta1BC','beta2','beta2BC'};
    for a = 1:6;
        export_fig(figure(a),[path_project 'data/midlevel/' names{a} '.png']);
    end
    
elseif strcmp(varargin{1},'FPSA_get_table_behavior');
    %%
    % Target: relate model betas (representing ellipsoidness) to subject's ratings and scr 'behavior'.
    % Steps:
    % collect necessary data
    % set up table
    force    = 0;
    p        = Project;
    subs     = FPSA_FearGen('get_subjects');
    path2table = sprintf('%sdata/midlevel/table_predict_behavior_N%d.mat',path_project,length(subs));
    
    if ~exist(path2table)||force == 1
        %% prepare scr data
        scrsubs  = ismember(subs,p.subjects(p.subjects_scr));
        scrpath  = sprintf('%sdata/midlevel/SCR_N%d.mat',path_project,sum(scrsubs));
        %get scr data
        if ~exist(scrpath)
            g        = Group(subs(scrsubs));
            out      = g.getSCR(2.5:5.5);
            save(scrpath,'out');
            clear g
        else
            load(scrpath)
        end
        scr_test_nonparam            = nan(length(subs),1); % the table needs same number of rows, so we just fill a column of nans with scr params.
        scr_test_nonparam(scrsubs,:) = mean(out.y(:,[21 22 23]),2)-mean(out.y(:,[25 26 19]),2); %% diff between CSP and CSN (IS THIS ZSCORED?)
        
        %% scr data with fits        
        ns = 0;
        scr_test_parametric            = nan(length(subs),1);
        for sub = subs(scrsubs(:)');
            ns                          = ns + 1;
            s                           = Subject(sub);
            bla(ns,1)                   = s.get_fit('scr',4).params(1);
        end
        scr_test_parametric(scrsubs,1) = bla;
        %% prepare rating data
        % collect rating amplitudes
        ns = 0;
        for sub = subs(:)'
            ns     = ns+1;
            s      = Subject(sub);
            rating = s.get_rating(4);
%             Y      = zscore(rating.y);%zscored rating
%             Y      = accumarray(rating.x'/45+4,Y,[8 1],@mean)';
            Y = rating.y_mean;
%             amp_test_nonparam(ns,1)   = 
            rating_test_parametric(ns,1) = s.get_fit('rating',4).params(1);
            rating_test_nonparam(ns,1)   = mean(Y([3 4 5]))-mean(Y([1 7 8]));
        end
        %% get model parameters
        C          = FPSA_FearGen('FPSA_model_singlesubject',1:100);
        beta1         = C.model_02.w1(:,2);
        beta2         = C.model_02.w2(:,2);        
        %% concatenate everything in the table
        t = table(subs(:),rating_test_parametric,rating_test_nonparam,scr_test_parametric,scr_test_nonparam,beta1,beta2,'variablenames',{'subject_id' 'rating_test_parametric','rating_test_nonparam','scr_test_parametric','scr_test_nonparam','beta1','beta2'});
        save(path2table,'t');
    else
        fprintf('Found table at %s, loading it.\n',path2table)
        load(path2table);
    end
    %%
    varargout{1} = t;
elseif strcmp(varargin{1},'model_fpsa_testgaussian_optimizer');
    %% create Gaussian models with different parameters to find the best one to compare against the flexible model
    t           = FPSA_FearGen('FPSA_get_table',1:100);
    tsubject    = length(unique(t.subject));
    res         = 50;
    amps        = linspace(0.25,2,res);
    sds         = linspace(0.25,3,res);
    c           = 0;
    
    for amp = amps
        for sd = sds
            c          = c + 1;
            %
            [cmat]     = getcorrmat(0,amp,0,1,sd);
            imagesc(cmat,[-1 1]);
            colorbar;title('Currently fitted Gau component');drawnow;
            model3_g   = Vectorize(repmat(repmat(squareform_force(cmat),1,1),1,tsubject));%
            %
            t.gau      = model3_g(:);
            a          = fitlm(t,'FPSA_G ~ 1 + specific + unspecific + Gaussian');
            BIC2(c)    = a.ModelCriterion.BIC;
        end
    end
    BIC2         = reshape(BIC2,res,res);
    varargout{1} = BIC2;
    
    %% prepare the output;
    [y x]  = find(BIC2 == min(BIC2(:)));
    amp    = amps(x);
    sd     = sds(y);
    clf
    imagesc(amps,sds,reshape(BIC2,res,res));
    hold on
    plot(amp,sd,'ko','markersize',25);
    hold off;
    title( 'BIC = f(amp,sd)');
    xlabel('amp');
    ylabel('sd');
    fprintf('The best amplitude and sd parameters are as follows:\n AMP: %3.5g, SD: %3.5g\n',amp,sd);
    
elseif strcmp(varargin{1},'model_fpsa_parameter_timecourse')
    for nfix = 1:5;
        a{nfix} = FPSA_FearGen('model_fpsa_testflexible',nfix);
    end;
    e=[];
    for nfix = 1:5;
        e(nfix,:) = a{nfix}.Coefficients.Estimate(2:end);
    end;
    bar(e)
    
    
elseif strcmp(varargin{1},'NvsO')
    %% compares neighboring correlation to opposing correlations
    % BLOCK = 1 | 2 for baseline | test, respectively
    sim   = varargin{2};
    block = varargin{3};
    r     = FPSA_FearGen('get_block',sim,block,block);
    for ns = 1:size(r,3)
        c.N(:,ns) = diag(r(:,:,ns),1);
        c.O(:,ns) = diag(r(:,:,ns),4);
    end
    varargout{1} = c;
    [h p stats bla] = ttest(fisherz(1-mean(c.N))-fisherz(1-mean(c.O)))
    
    
elseif strcmp(varargin{1},'CompareB2T_RSA')
    %% returns the coordinates and pvalue of the test comparing corresponding similarity entries of between baseline and generalization.
    %the full B and T similarity matrix;
    sim = FPSA_FearGen('get_fpsa_fair',1:100,1:3);%
    %%we only want the B and T parts
    [~,B] = FPSA_FearGen('get_block',sim,1,1);
    [~,T] = FPSA_FearGen('get_block',sim,2,2);
    %fisher transform and make a ttest
    [h p ] = ttest(fisherz(B)-fisherz(T));
    h      = squareform_force(h);
    p      = squareform_force(p);
    [i]    = find(p < .05);
    p      = p(i);
    varargout{1} = [i p];
    
    
elseif strcmp(varargin{1},'get_betas')
    %% compute loadings on these
    sim    = varargin{2};
    tsub   = size(sim.correlation,1);
    tblock = length(squareform(sim.correlation(1,:)))/8;
    fprintf('Found %02d blocks\n',tblock);
    betas  = [];
    X      = FPSA_FearGen('get_design_matrix');
    for nblock = 1:tblock
        n      = 0;
        data   = FPSA_FearGen('get_block',sim,nblock,nblock);
        while n < tbootstrap
            n                  = n +1;
            i                  = randsample(1:tsub,tsub,1);
            Y                  = ( 1-squareform(nanmean(data(:,:,i),3)) );
            betas(n,:,nblock)  = X\Y';
        end
    end
    % get errorbars for that
    ci           = prctile(betas,[2.5 97.5]);
    varargout{1} = squeeze(mean(betas));
    varargout{2} = ci;
    %
elseif strcmp(varargin{1},'test_betas')
    %%
    sim     = FPSA_FearGen('get_fpsa',1:100);
    [a b c] = FPSA_FearGen('get_betas_singlesubject',sim);
    [h p]   = ttest(c(:,2,1)-c(:,2,2));
    fprintf('t-test for physical similarity H: %d, p-value:%3.5g\n',h,p);
    [h p]   = ttest(c(:,3,1)-c(:,3,2));
    fprintf('t-test for cs+ similarity H     : %d, p-value:%3.5g\n',h,p);
elseif strcmp(varargin{1},'get_betas_singlesubject')
    %% compute loadings on these
    sim    = varargin{2};
    tsub   = size(sim.correlation,1);
    tblock = length(squareform(sim.correlation(1,:)))/8;
    fprintf('Found %02d blocks\n',tblock);
    betas  = [];
    X      = FPSA_FearGen('get_design_matrix');
    for nblock = 1:tblock
        data   = FPSA_FearGen('get_block',sim,nblock,nblock);
        for n = 1:tsub
            Y                  = ( 1-squareform(mean(data(:,:,n),3)) );
            betas(n,:,nblock)  = X\Y';
        end
    end
    % get errorbars for that
    ci           = nanstd(betas)./sqrt(tsub);
    varargout{1} = squeeze(nanmean(betas));
    varargout{2} = [nanmean(betas)-ci/2 ;nanmean(betas)+ci/2];
    varargout{3} = betas;
    %
    
elseif strcmp(varargin{1},'searchlight')
    %% conducts a searchlight analysis on the FDMs using a moving window of about 1 degrees
    % Default window parameters B1 and B2 are 1, and 15; (1 degree running
    % average windows with full overlap).
    % At each searchlight position the flexible model is fit.
    b1                = varargin{2};
    b2                = varargin{3};
    fixations         = 1:100;
    runs_per_phase{2} = 1;
    runs_per_phase{4} = runs;
    fun               = @(block_data) FPSA_FearGen('fun_handle',block_data.data);%what we will do in every block
    
    runc             = 0;%1 run from B + 3 runs from T.
    for phase = [2 4];
        conds = condition_borders{phase};%default parameter
        for run = runs_per_phase{phase}
            runc             = runc + 1;
            fixmat           = FPSA_FearGen('get_fixmat','runs',run);%get the fixmat for this run
            filename         = DataHash({fixmat.kernel_fwhm,b1,b2,phase,run,fixations});
            subc = 0;
            for subject = unique(fixmat.subject);
                subc                 = subc + 1;%subject counter
                path_write           = sprintf('%ssub%03d/p%02d/midlevel/%s.mat',path_project,subject,phase,filename);
                cprintf([1 0 0],'Processing subject %03d\ncache name: %s\n',subject,path_write);
                if exist(fileparts(path_write)) == 0;mkdir(fileparts(path_write));end;%create midlevel folder if not there.
                %analysis proper
                if exist(path_write) == 0 | force
                    % create the query cell
                    maps             = FPSA_FearGen('get_fixmap',fixmat,subject,fixations);
                    maps             = reshape(maps(:,conds),[500 500 length(conds)]);
                    out              = blockproc(maps,[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true,'UseParallel',true,'DisplayWaitbar',false);
                    save(path_write,'out');
                else
                    cprintf([0 1 0],'Already cached...\n');
                    load(path_write);
                end
                B1(:,:,:,subc,runc) = out;
            end
        end
    end
    varargout{1} = B1;
elseif strcmp(varargin{1},'fun_handle')
    %% This is the function kernel executed for each seachlight position.
    
    maps = varargin{2};
    maps = reshape(maps,[size(maps,1)*size(maps,2) size(maps,3)]);
    if all(sum(abs(maps)))
        Y            = 1-pdist(maps','correlation');
        X            = FPSA_FearGen('get_design_matrix');%returns the design matrix for the flexible ellipsoid model
        betas(1,1,:) = X\Y';
    else
        betas(1,1,:)= [NaN NaN NaN];
    end
    varargout{1} = betas;
    
elseif strcmp(varargin{1},'get_design_matrix');
    %% Design matrix for the ellipsoide model.
    
    x          = [pi/4:pi/4:2*pi];
    w          = [cos(x);sin(x)];
    %
    model2_c   = squareform_force(cos(x)'*cos(x));
    model2_s   = squareform_force(sin(x)'*sin(x));    
    X          = [ones(length(model2_c(:)),1) model2_c model2_s];        
    varargout{1}  = X;
    
    
elseif strcmp(varargin{1},'plot_searchlight')
    %% will spatially plot the ellipses area using the sqrt(cosine^2 +
    %sine^2) weight combination term.
    
    Mori            = FPSA_FearGen('searchlight',1,15);
    %MORI(pixel,pixle,regressor,subject,phase)
    %%
    M                = nanmean(Mori,4);%average across subjects.
    M(:,:,1,:,:)     = sqrt(M(:,:,2,:,:).^2+M(:,:,3,:,:).^2);%average across sine and cosine
    M(:,:,2:end,:,:) = [];%delete cos and sine;
    M(:,:,:,:,2)     = mean(M(:,:,:,:,2:end),5);%merge the 3 test runs.
    M(:,:,:,:,3:end) = [];
    %%
    fixmat           = Fixmat([],[]);
    fixmat.maps      = squeeze(M);
    fixmat.plot;    
    %%        
    figure;
    B = fixmat.EyeNoseMouth(M(:,:,1,1,2))-fixmat.EyeNoseMouth(M(:,:,1,1,1));bar(B(1:4))
    
    
elseif strcmp(varargin{1},'searchlight_stimulus')
    %% applies the search light analysis to the V1 representations.
    b1          = varargin{2};
    b2          = varargin{3};
    noise_level = varargin{4};
    filename    = 'stimulus_searchlight';
    path_write  = sprintf('%smidlevel/%s_noiselevel_%02d.mat',path_project,filename,noise_level);
    fun         = @(block_data) FPSA_FearGen('fun_handle',block_data.data);%what we will do in every block
    maps        = [];
    for n = 1:8
        maps(:,:,n) = imread(sprintf('%sstimuli/%02d.bmp',path_project,n));
    end
    obj  = Fixmat([],[]);
    maps = maps + rand(size(maps))*noise_level;
    maps = maps( obj.rect(1):obj.rect(1)+obj.rect(3)-1,  obj.rect(2):obj.rect(2)+obj.rect(4)-1,:);
    
    if exist(path_write) == 0
        % create the query cell
        out              = blockproc(maps,[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true,'UseParallel',true);
        save(path_write,'out');
    else
        cprintf([0 1 0],'Already cached...\n');
        load(path_write);
    end
    varargout{1} = out;
    %
    %     subplot(1,2,1);
    figure;
    imagesc(out(:,:,2));
    %     hold on
    %     f   = Fixmat([],[]);
    %     roi = f.GetFaceROIs;
    %     [~,h] = contourf(mean(roi(:,:,1:4),3));
    %     h.Fill = 'off';
    %     axis image;
    %     hold off;
    %     subplot(1,2,2);
    %     b = f.EyeNoseMouth(out(:,:,2),0)
    %     bar(b(1:4));
elseif strcmp(varargin{1},'searchlight_bs')
    %%
    fixmat   = varargin{2}
    b1       = varargin{3};%1
    b2       = varargin{4};%15
    %
    tsub     = length(unique(fixmat.subject));
    fun      = @(block_data) FPSA_FearGen('fun_handle',block_data.data);%what we will do in every block
    bs       = 0;
    while bs < 1000
        bs                = bs+1;
        fprintf('Processing bs %03d\n',bs);
        % craete the query cell
        subject          = randsample(1:tsub,tsub,1);
        maps             = FPSA_FearGen('get_fixmap',fixmat,subject,1:100);
        maps             = reshape(maps,[500 500 16]);
        B1(:,:,:,bs,1)   = blockproc(maps(:,:,1:8),[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true,'UseParallel',true);
        B1(:,:,:,bs,2) = blockproc(maps(:,:,9:16),[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true,'UseParallel',true);
        c = 0;
        for m = 1:2
            for n = 1:3
                c = c +1;
                subplot(2,3,c)
                imagesc(nanmean(B1(:,:,n,:,m),4));colorbar;
                drawnow;
            end
        end
    end
    varargout{1} = B1;
    
elseif strcmp(varargin{1},'beta_counts')
    %% counts searchlight values in face rois
    fixmat      = FPSA_FearGen('get_fixmat');
    b1          = 1;
    b2          = 15;
    out         = FPSA_FearGen('searchlight',fixmat,b1,b2);
    for np = 1:2
        for ns = 1:size(out,4);
            for beta = 2%1:size(out,3)
                map            = out(:,:,beta,ns,np);
                count(ns,:,np) = fixmat.EyeNoseMouth(map,0);
            end
        end
    end
    varargout{1} = count;
    cr           = count;
    mean((cr(:,1,1)-cr(:,2,1))-(cr(:,1,2)-cr(:,2,2)))
    
elseif strcmp(varargin{1},'anova')
    %%
    fixmat   = FPSA_FearGen('get_fixmat');
    cr       = FPSA_FearGen('beta_counts',fixmat,1,15);
    tsubject = length(unique(fixmat.subject))
    
    Y        = [cr(:,1,1) cr(:,2,1) cr(:,1,2) cr(:,2,2)];
    figure;
    errorbar(mean(Y),std(Y)./sqrt(size(Y,1)));
    y        = [cr(:,1,1);cr(:,2,1);cr(:,1,2);cr(:,2,2)];
    side     = [ones(tsubject,1);ones(tsubject,1)*2;ones(tsubject,1);ones(tsubject,1)*2];
    phase    = [ones(tsubject,1);ones(tsubject,1);ones(tsubject,1)*2;ones(tsubject,1)*2];
    anovan(y,{side(:) phase(:)},'model','full')
    
elseif strcmp(varargin{1},'eyebehave_params')
    
    savepath = sprintf('%s/data/midlevel/',path_project);
    filename = 'eyebehave_params.mat';
    subs = FPSA_FearGen('get_subjects');
    
    visualization = 1; %if you want all plots to be created
    
    if nargin > 1
        force = varargin{2};
    end
    
    if ~exist([savepath filename]) || force == 1
        
        fix = Fixmat(subs,[2 3 4]);
        
        %mean number of fixations for this combination
        [d.fixN.data, d.fixN.info] = fix.histogram;
        dummy = [];
        dummy2 = [];
        
        sc = 0;
        for sub = unique(fix.subject)
            fprintf('\nWorking on sub %02d, ',sub)
            sc= sc+1;
            pc = 0;
            for ph = 2:4
                fprintf('phase %d. ',ph)
                pc = pc+1;
                cc=0;
                for cond = unique(fix.deltacsp)
                    cc=cc+1;
                    ind = logical((fix.subject==sub).*(fix.phase == ph).* (fix.deltacsp == cond));
                    %mean duration of fixations for this phase/sub/cond
                    d.fixdur.m(sc,pc,cc) = mean(fix.stop(ind)-fix.start(ind));
                    %% for entropy, we need single trials, otherwise the trial number contributing to mean FDM (for this cond-phase-sub) biases the entropy computation
                    %mean entropy for this combination
                    fix.unitize = 0;
                    tc = 0;
                    for tr = unique(fix.trialid(ind)) %loop through trials of this cond-phase-sub
                        tc = tc+1;
                        fix.getmaps({'trialid' tr 'phase' ph 'subject' sub 'deltacsp' cond});
                        dummy(tc) = FDMentropy(fix.vectorize_maps);
                        dummy2(tc) = entropy(fix.vectorize_maps);
                    end
                    d.FDMentropy.m(sc,pc,cc) = mean(dummy);
                    d.entropy.m(sc,pc,cc) = mean(dummy2);
                    fix.unitize = 1;
                    dummy = [];
                    dummy2 = [];
                end
            end
        end
        save([savepath filename],'d','subs');
    else
        load([savepath filename]);
    end
    
    varargout{1} = d;
    
    if visualization ==1
        subtitles = {'base','cond','test'};
        %% Number of fixations per condition and phase
        figure;
        meancorr = 1;
        for n=1:3;
            if meancorr  == 1
                data = squeeze(d.fixN.data(:,:,n)) - repmat(nanmean(d.fixN.data(:,:,n),2),1,size(d.fixN.data,2));
            else
                data = squeeze(d.fixN.data(:,:,n));
            end
            
            subplot(2,3,n);
            boxplot(data);
            set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',60,'fontsize',14);
            title(subtitles{n})
            box off
            if meancorr ==1
                ylabel('N fixations (subj.mean corr.)')
            else
                ylabel('N fixations')
            end
            axis square
            subplot(2,3,n+3)
            Project.plot_bar(-135:45:315,nanmean(data),nanstd(data)./sqrt(length(unique(fix.subject))));
            set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',60,'fontsize',14);
            
            box off
            if meancorr ==1
                ylabel('N fixations (subj.mean corr.)')
            else
                ylabel('N fixations')
            end
            axis square
        end
        for n = 4:6
            subplot(2,3,n)
            ylim([min(mean(data))-.5 max(mean(data))+.7])
            %             ylim([nanmean(d.fixN.data(:))-nanstd(d.fixN.data(:)) nanmean(d.fixN.data(:))+nanstd(d.fixN.data(:))])
        end
        %% Fixation durations per condition
        figure;
        meancorr = 1;
        for n=1:3;
            if meancorr  ==1
                data = squeeze(d.fixdur.m(:,n,:)) - repmat(nanmean(d.fixdur.m(:,n,:),3),1,size(d.fixdur.m,3));
            else
                data = squeeze(d.fixdur.m(:,n,:));
            end
            subplot(2,3,n);
            boxplot(data);
            set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',60,'fontsize',14);
            title(subtitles{n})
            box off
            ylabel('fix duration [ms]')
            axis square
            subplot(2,3,n+3)
            Project.plot_bar(-135:45:315,nanmean(data),nanstd(data)./sqrt(length(subs)));
            set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',60,'fontsize',14);
            box off
            ylabel('fix duration [ms]')
            axis square
        end
        for n = 4:6
            subplot(2,3,n)
            ylim([nanmean(data(:))-nanstd(data(:)) nanmean(data(:))+nanstd(data(:))])
        end
        % difference of fixation durations between test and baseline
        figure;
        diffTB = squeeze(d.fixdur.m(:,3,:)-d.fixdur.m(:,1,:)); % this is practically meancorrected for the sub then.. we just look at the difference.
        subplot(1,2,1)
        boxplot(diffTB);
        set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',60,'fontsize',14);
        hold on;
        line(xlim,[0 0])
        box off
        ylabel('duration diff Test-Base [ms]')
        axis square
        subplot(1,2,2)
        Project.plot_bar(-135:45:315,nanmean(diffTB),nanstd(diffTB)./sqrt(length(subs)));
        set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',60,'fontsize',14);
        hold on;
        line(xlim,[0 0])
        box off
        ylabel('duration diff Test-Base [ms]')
        axis square
        
        %% Entropy of fixmaps per condition, sub and phase (matlab entropy)
        figure;
        meancorr = 1;
        for n = 1:3;
            if meancorr  ==1
                data = squeeze(d.FDMentropy.m(:,n,:)) - repmat(nanmean(d.FDMentropy.m(:,n,:),3),1,size(d.FDMentropy.m,3));
            else
                data = squeeze(d.FDMentropy.m(:,n,:));
            end
            subplot(2,3,n);
            boxplot(data);
            set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',45,'fontsize',14);
            title(subtitles{n})
            box off
            ylabel('FDM entropy [a.u.]')
            axis square
            subplot(2,3,n+3);
            Project.plot_bar(-135:45:315,nanmean(data),nanstd(data)./sqrt(length(subs)));
            set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',45,'fontsize',14);
            box off
            ylabel('FDM entropy [a.u.]')
            axis square
            if meancorr == 1 && n==3
                subplot(2,3,6)
                dada.x = repmat(-135:45:180,74,1);
                dada.y = data(:,1:8)*1000; % Tuning obj not working for such small numbers.
                dada.ids = 1:74;
                t = Tuning(dada);t.GroupFit(3);
                plot(t.groupfit.x_HD,[t.groupfit.fit_HD+mean(t.y(:))]./1000,'k','LineWidth',2)
            end
            
        end
        for n = 4:6
            subplot(2,3,n)
            ylim([nanmean(data(:))-nanstd(data(:)) nanmean(data(:))+nanstd(data(:))])
        end
        
        
        % %        %% Entropy of fixmaps per condition, sub and phase (own entropy)
        % %        function E = FDMentropy(fdm)
        % %        % computes entropy of a fixation density map.
        % %        % Map should be normalized anyway. If not, this function does it.
        % %
        % %        % remove zero entries in p
        % %        fdm(fdm==0) = [];
        % %
        % %        if sum(fdm) ~=0
        % %            % normalize p so that sum(p) is one.
        % %            fdm = fdm ./ numel(fdm);
        % %        end
        % %
        % %        E = -sum(fdm.*log2(fdm));
        % %        end
        %
        %     figure;
        %     for n=1:3;
        %         subplot(2,3,n);
        %         boxplot(squeeze(d.FDMentropy.m(:,n,:)));
        %         set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',45,'fontsize',14);
        %         title(subtitles{n})
        %         box off
        %         ylabel('FDM entropy [a.u.]')
        %         axis square
        %         subplot(2,3,n+3);
        %         Project.plot_bar(-135:45:315,nanmean(squeeze(d.FDMentropy.m(:,n,:))),nanstd(squeeze(d.FDMentropy.m(:,n,:)))./sqrt(length(unique(fix.subject))));
        %         set(gca,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 'Odd' 'Null'},'XTickLabelRotation',45,'fontsize',14);
        %         box off
        %         ylabel('FDM entropy [a.u.]')
        %         axis square
        %     end
        %     for n = 4:6
        %         subplot(2,3,n)
        %         ylim([nanmean(d.FDMentropy.m(:))-nanstd(d.FDMentropy.m(:)) nanmean(d.FDMentropy.m(:))+nanstd(d.FDMentropy.m(:))])
        %     end
        
        %%
        figure;
        colors = GetFearGenColors;
        colors = [colors; .8 .8 .8];
        set(0,'DefaultAxesColorOrder',colors);
        subplot(1,3,1)
        title('Num Fixations','fontsize',14);hold on
        errorbar([repmat(1:3,11,1)'+rand(3,11)*.2],squeeze(mean(d.fixN.data))',squeeze(nanstd(d.fixN.data))'./sqrt(length(d.fixN.data)),'LineWidth',1.5,'LineStyle',':')
        box off;
        axis square;
        ylabel('Num Fixations M/SEM')
        subplot(1,3,2)
        title('Fixation duration','fontsize',14);hold on
        errorbar([repmat(1:3,11,1)+rand(11,3)*.2]',squeeze(mean(d.fixdur.m)),squeeze(nanstd(d.fixdur.m))./sqrt(length(d.fixdur.m)),'LineWidth',1.5,'LineStyle',':')
        box off;
        axis square;
        ylabel('Fixation duration[ms] M/SEM')
        subplot(1,3,3)
        title('Entropy','fontsize',14);hold on
        errorbar([repmat(1:3,11,1)+rand(11,3)*.2]',squeeze(mean(d.entropy.m)),squeeze(nanstd(d.entropy.m))./sqrt(length(d.entropy.m)),'LineWidth',1.5,'LineStyle',':')
        box off;
        axis square;
        ylabel('FDM entropy [a.u.] M/SEM')
        legend({'-135','-90','-45','CS+','45','90','135','CS-','UCS','Odd','Null'},'fontsize',14,'location','best')
    end
elseif strcmp(varargin{1},'inter-subject-variance')
    
    subs = FPSA_FearGen('get_subjects');
    
    pc = 0;
    for ph = 2:4
        pc = pc+1;
        sc = 0;
        for sub = unique(subs)
            sc = sc+1;
            fix = Fixmat(sub,ph);
            cc = 0;
            for cond = fix.realcond
                cc = cc+1;
                ind = fix.deltacsp == cond;
                tc = 0;
                for tr = unique(fix.trialid(ind)) %loop through trials of this cond-phase-sub
                    tc = tc+1;
                    fix.getmaps({'trialid' tr 'phase' ph 'subject' sub 'deltacsp' cond});
                    dummy(tc) = FDMentropy(fix.vectorize_maps);
                    dummy2(tc) = entropy(fix.vectorize_maps);
                end
            end
        end
        fix.getsubmaps;
        submaps  = fix.maps;
        groupave = mean(submaps,3);
        fix.maps = [submaps - repmat(groupave,1,1,length(subs))];
        fix.maps = imresize(fix.maps,.1,'method','bilinear');
        
        ISV(pc,:) = var(fix.vectorize_maps);
        
    end
    
elseif strcmp(varargin{1},'SVM')
    %This script trains a linear SVM training CS+ vs. CS- for phases 2 and 4.
    %It collects the data and computes the eigenvalues on the
    %fly for chosen(or a range of parameters) (kernel_fwhm, number of
    %eigenvalues). As the number of trials are lower in the
    %baseline, all the test-training sessions should use the same number of
    %trials. For example to keep the comparisons comparable, in baseline 11
    %trials in the baseline, with .5 hold-out, one needs to sample the same
    %number of trials from the testphase before training.
    %the option random = 1 randomizes labels to determine the chance classification
    %performance level.
    random = 0;
    exclmouth = 1;
    tbootstrap       = 1000; %number of bootstraps
    phase            = [2 4 4 4];%baseline = 2, test = 4
    holdout_ratio    = .5; %holdout_ratio for training vs. test set
    teig             = 100; %up to how many eigenvalues should be included for tuning SVM?
    crit             = 'var';%choose 'ellbow' classification or 'var' 90% variance explained.
    cutoffcrit       = .9;
    R                = [];%result storage for classification performance
    HP               = [];%result storage for single subject hyperplanes.
    AVEHP            = [];%result storate for average hyperplane
    
    eigval           = [];
    trialselect      = {1:120 1:120 121:240 241:360};
    
    
    subjects = FPSA_FearGen('get_subjects');
    o = 0;
    for run = 1:4 % phase 2, phase 4.1 phase 4.2 phase 4.3
        o = o+1;
        fix             = Fixmat(subjects,phase(run));%get the data
        if exclmouth == 1
            roi = fix.GetFaceROIs;
        end
        fix.unitize     = 1;%unitize fixmaps or not (sum(fixmap(:))=0 or not).
        %% get number of trials per condition and subject: Sanity check...
        M               = [];%this will contain number of trials per subject and per condition. some few people have 10 trials (not 11) but that is ok. This is how I formed the subject_exlusion variable.
        sub_c           = 0;
        for ns = subjects(:)'
            sub_c = sub_c + 1;
            nc = 0;
            for cond = unique(fix.deltacsp)
                nc          = nc +1;
                i           = ismember(fix.phase,phase(run)).*ismember(fix.subject,ns).*ismember(fix.deltacsp,cond);%this is on fixation logic, not trials.
                i           = logical(i);
                M(sub_c,nc,o) = length(unique(fix.trialid(i)));
            end
        end
        %% get all the single trials in a huge matrix D together with labels.
        global_counter = 0;
        clear D;%the giant data matrix
        clear labels;%and associated labels.
        for ns = subjects(:)'
            for deltacsp = -135:45:180;
                i              = ismember(fix.trialid,trialselect{run}).*(fix.subject == ns).*(fix.deltacsp == deltacsp);
                trials         = unique(fix.trialid(i == 1));
                trial_counter  = 0;
                for trialid = trials
                    trial_counter       = trial_counter + 1;
                    global_counter      = global_counter +1;
                    c                   = global_counter;
                    v                   = {'subject' ns 'deltacsp' deltacsp 'trialid' trialid};
                    fix.getmaps(v);
                    if exclmouth == 1
                        fix.maps(roi(:,:,4)) = 0;
                    end
                    D(:,c)              = Vectorize(imresize(fix.maps,.1));
                    labels.sub(c)       = ns;
                    labels.phase(c)     = phase(run);
                    labels.trial(c)     = trial_counter;%some people have less trials check it out with plot(labels.trial)
                    labels.cond(c)      = deltacsp;
                end
            end
        end
        %% DATA2LOAD get the eigen decomposition: D is transformed to TRIALLOAD
        fprintf('starting covariance computation\n')
        covmat    = cov(D');
        fprintf('done\n')
        fprintf('starting eigenvector computation\n')
        [e dv]    = eig(covmat);
        fprintf('done\n')
        dv        = sort(diag(dv),'descend');
        eigval(:,run) = dv;
        %     figure(100);
        %     plot(cumsum(dv)./sum(dv),'o-');xlim([0 200]);drawnow
        eigen     = fliplr(e);
        %collect loadings of every trial
        trialload = D'*eigen(:,1:teig)*diag(dv(1:teig))^-.5;%dewhitened
        %% LIBSVM business
        neigs = [7 12 8 10]; %check eigenvalues and put numbers of EV here, based on ellbow criterion.
        if strcmp(crit,'ellbow')
            neig = neigs(run);
        elseif strcmp(crit,'var')
            neig = find(cumsum(dv)./sum(dv)>cutoffcrit,1,'first');
        end
        sub_counter = 0;
        result      = [];
        w           = [];
        for sub = unique(labels.sub)%go subject by subject
            fprintf('run:%d-Eig:%d-Sub:%d\n',run,neig,sub);
            if random == 1
                warning('Randomizing labels as wanted. \n')
            end
            sub_counter = sub_counter + 1;
            ind_all     = ismember(labels.sub,sub);%this subject, this phase.
            %
            for n = 1:tbootstrap%
                Ycond   = double(labels.cond(ind_all))';%labels of the fixation maps for this subject in this phase.
                X       = trialload(ind_all,1:neig);%fixation maps of this subject in this phase.
                % now normal Holdout for every phase (which should all have the
                % same number of trials now)
                P       = cvpartition(Ycond,'Holdout',holdout_ratio); % divide training and test datasets respecting conditions
                i       = logical(P.training.*ismember(Ycond,[0 180]));%train using only the CS+ and CS? conditions.
                if random ==1
                    model   = svmtrain(Shuffle(Ycond(i)), X(i,1:neig), '-t 0 -c 1 -q'); %t 0: linear, -c 1: criterion, -q: quiet
                else
                    model   = svmtrain(Ycond(i), X(i,1:neig), '-t 0 -c 1 -q'); %t 0: linear, -c 1: criterion, -q: quiet
                end
                % get the hyperplane
                try
                    w(:,sub_counter,n)          = model.SVs'*model.sv_coef;
                catch
                    keyboard%sanity check: stop if something is wrong
                end
                %%
                cc=0;
                for cond = unique(Ycond)'
                    cc                          = cc+1;
                    i                           = logical(P.test.*ismember(Ycond,cond));%find all indices that were not used for training belonging to COND.
                    [~, dummy]                  = evalc('svmpredict(Ycond(i), X(i,:), model);');%doing it like this supresses outputs.
                    dummy                       = dummy == 0;%binarize: 1=CS+, 0=Not CS+
                    result(cc,n,sub_counter)    = sum(dummy)/length(dummy);%get the percentage of CS+ responses for each CONDITION,BOOTSTR,SUBJECT
                end
            end
        end
        %once the data is there compute relevant output metrics:
        R(:,run)      = mean(mean(result,2),3);%average across bootstaps the classification results
        AVEHP(:,:,run) = reshape(mean(eigen(:,1:neig)*mean(w,3),2),[50 50 1]);%average hyperplane across subjects
        HP(:,:,:,run) = reshape(eigen(:,1:neig)*mean(w,3),[50 50 size(eigen(:,1:neig)*mean(w,3),2)]); %single hyperplanes
        
        savepath = sprintf('%s/data/midlevel/SVM/',path_project);
        filename = sprintf('/SVM_NEV%d_FWHM30_r%d_run%d_crit%s_exclmouth_%d.mat',neig,random,run,crit,exclmouth);
        if exist(savepath)
            save([savepath filename],'neig','R','eigval','result','HP','AVEHP');
        else
            fprintf('Creating SVM results folder...\n')
            mkdir(savepath);
            save([savepath filename],'neig','R','eigval','result','HP','AVEHP');
        end
    end
    %% check mouth exclusion
    savepath = sprintf('%s/data/midlevel/SVM/',path_project);
    files = cellstr(ls(savepath));
    %     neigs = [63 69 74 75 14 19 17 19];
    mouthex = [0 0 0 0  1 1 1 1];
    run = [1:4 1:4];
    for c = 1:8
        expr = sprintf('r0_run%d_critvar_exclmouth_%d.mat',run(c),mouthex(c));
        findfile = regexp(files,expr,'match');
        ind = find(~cellfun(@isempty,findfile));
        load([savepath files{ind}],'result');
        fprintf('Loading file %s\n',[savepath files{ind}])
        results(:,:,c) = squeeze(mean(result,2));
    end
    M = squeeze(mean(results,2));
    SE = squeeze(std(results,[],2)./sqrt(size(results,2)));
    
    % average the three test runs
    M = cat(2,M(:,1:4),mean(M(:,2:4),2),M(:,5:8),mean(M(:,6:8),2));
    SE = cat(2,SE(:,1:4),mean(SE(:,2:4),2),SE(:,5:8),mean(SE(:,6:8),2));
    
    ind = [1 5 6 10];
    ylab = {'Mouth INCL' '' 'Mouth EXCL' ''};
    for n = 1:4
        subplot(2,2,n)
        xlim([-170 215]);l=line(xlim,[.5 .5]); hold on;set(l,'Color','k','LineStyle',':');
    end
    for n = 1:4
        subplot(2,2,n)
        Project.plot_bar(-135:45:180,M(:,ind(n)),SE(:,ind(n)));
        hold on;
        ylim([.3 .7])
        set(gca,'YTick',.3:.1:.7,'XTick',[0 180],'XTickLabel',{'CS+' 'CS-'},'FontSize',14);
        box off
        axis square
        ylabel(ylab{n})
        set(gca,'YTick',[.3 .5 .7])
        xlim([-170 215])
    end
elseif strcmp(varargin{1},'get_trials_for_svm_idiosync')
    % here we classify subjects to show their idiosyncratic patterns.
    % we classify them within each phase seperately to account for diff.
    % number of trials, i.e. baseline, test1,test2,test3
    r             = 0; % for label randomization check
    name_analysis = 'subjects_inphase'; %classify subjects, respect phases
    savepath      = fullfile(path_project,['data\midlevel\SVM_idiosyncrasy\']);
    filename      = 'trialdata.mat';
    savedfile     = [savepath filename];
    %   savepath      = fullfile(path_project,['data\midlevel\SVM_idiosyncrasy_' name_analysis '_rand' num2str(r) '.mat']);
    if nargin >1
        force = varargin{2};
    end
    
    if ~exist(savepath)
        mkdir(savepath)
    end
    
    if ~exist(savedfile) || force==1
        if ~exist(savedfile)
            fprintf('Trial data not yet found at %s, \ncollecting it now. \n',savepath)
            WaitSecs(2); %make it readable.
        elseif force == 1
            fprintf('Forced to compute this file anew: \n %s',savepath)
            WaitSecs(2); %make it readable.
        end
        
        scale         = .1;
        run2ph        = [2 4 4 4];
        run2trial     = {1:120 1:120 121:240 241:360};
        subs          = FPSA_FearGen('get_subjects');
        nsubs         = length(subs);
        
        cutoffcrit     = .9; %percentage of variance that should be explained.
        
        fprintf('Started analysis (%s): %s\n',datestr(now,'hh:mm:ss'),name_analysis);
        
        ttrial           = nsubs*length(cell2mat(run2trial));
        D                = NaN(2500,ttrial);
        labels.sub       = NaN(1,ttrial);
        labels.run       = NaN(1,ttrial);
        labels.trial     = NaN(1,ttrial);
        labels.cond      = NaN(1,ttrial);
        v = [];
        c = 0;
        for sub = subs(:)'
            for run = 1:4
                fix = Fixmat(sub,run2ph(run));
                for tr = run2trial{run}
                    v = {'subject' sub, 'trialid' tr 'deltacsp' fix.realcond};
                    fprintf('Subject %d run %d trial %d\n',sub,run,tr);
                    fix.getmaps(v);
                    if ~any(isnan(fix.maps(:)))
                        c                   = c+1;
                        %scale it if necessary
                        if scale ~= 1
                            fix.maps        = imresize(fix.maps,scale,'method','bilinear');
                        end
                        D(:,c)              = fix.vectorize_maps;
                        labels.sub(c)       = sub;
                        labels.run(c)       = run;
                        labels.trial(c)     = tr;
                        labels.cond(c)      = unique(fix.deltacsp(fix.selection));
                    end
                end
            end
        end
        
        %cut the nans
        todelete = isnan(labels.sub); %could be derived from labels.anything, means empty fixmap
        fprintf('Will delete %g trials...\n',sum(todelete));
        D(:,todelete)=[];
        labels.sub(:,todelete)=[];
        labels.run(:,todelete)=[];
        labels.trial(:,todelete)=[];
        labels.cond(:,todelete)=[];
        
        c = 0;
        for l = unique(labels.sub)
            c = c + 1;
            labels.easy_sub(labels.sub == l) = c;
        end
        
        %% DATA2LOAD get the eigen decomposition: D is transformed to TRIALLOAD
        fprintf('starting covariance computation\n')
        covmat    = cov(D');
        fprintf('done\n')
        fprintf('starting eigenvector computation\n')
        [e dv]    = eig(covmat);
        fprintf('done\n')
        dv        = sort(diag(dv),'descend');
        eigval(:,run) = dv;
        %     figure(100);
        %     plot(cumsum(dv)./sum(dv),'o-');xlim([0 200]);drawnow
        eigen     = fliplr(e);
        neig = find(cumsum(dv)./sum(dv)>cutoffcrit,1,'first');
        
        fprintf('%002d Eigenvectors explain 90% of variance. \n',neig)
        
        %collect loadings of every trial
        trialload = D'*eigen(:,1:neig)*diag(dv(1:neig))^-.5;%dewhitened
        
        
        save(savedfile)
    else
        load(savedfile);
        fprintf('Was already saved, loading it from \n%s\n',savedfile)
        Ndim = size(trialload);
        fprintf('Found a Matrix of dim = %03d (trials) x %03d (EV loadings).\n',Ndim(1),Ndim(2))
        fprintf('Labels show that we have %02d subjects, which are the following:\n',length(unique(labels.sub)));
        unique(labels.sub)
        WaitSecs(1); %make it readable.
    end
    varargout{1} = trialload;
    varargout{2} = labels;
    
elseif strcmp(varargin{1},'svm_howmanytrialspersub')
    
    [~,labels] = FPSA_FearGen('get_trials_for_svm_idiosync');
    for run = 1:4
        ind = labels.run == run;
        ntrials(:,run)    = histc(labels.easy_sub(ind),unique(labels.easy_sub));
    end
    disp([unique(labels.easy_sub)' ntrials])
    
    fprintf('Minimum number of trials is (run 1 2 3 4):\n')
    [minN, minisub] = min(ntrials);
    disp(minN)
    
    fprintf('This should be fed into SVM classifying subjects to make it fairer.\n')
    
    
    varargout{1} = minN;
    
elseif strcmp(varargin{1},'svm_classify_subs_1vs1')
    
    nbootstrap = 1000;
    randomize  = 0;
    PHoldout    = .2;
    
    savepath      = fullfile(path_project,'data\midlevel\SVM_idiosyncrasy\');
    filename      = sprintf('performance_1vs1_r%d.mat',randomize);
    savedfile     = [savepath filename];
    %   savepath      = fullfile(path_project,['data\midlevel\SVM_idiosyncrasy_' name_analysis '_rand' num2str(r) '.mat']);
    if nargin > 1
        force = varargin{2};
    end
    
    if ~exist(savepath)
        mkdir(savepath)
    end
    
    if ~exist(savedfile) || force==1
        
        if ~exist(savedfile)
            fprintf('Trial data not yet found at %s, \ncollecting it now. \n',savepath)
            WaitSecs(2); %make it readable.
        elseif force == 1
            fprintf('Forced to compute this file anew: \n %s',savepath)
            WaitSecs(2); %make it readable.
        end
        
        [data,labels] = FPSA_FearGen('get_trials_for_svm_idiosync');
        nsub = max(unique(labels.easy_sub));
        ntrialspersub = FPSA_FearGen('svm_howmanytrialspersub');
        
        
        for run = 1:4
            runind = labels.run == run;
            for s1 = 1:nsub
                for s2 = 1:nsub
                    if s1 < s2;
                        fprintf('Run: %02d. Classifying sub %02d vs %02d.\n',run,s1,s2)
                        for n = 1:nbootstrap
                            select    = logical(ismember(labels.easy_sub,[s1 s2]).*runind);
                            Y         = labels.easy_sub(select)';
                            X         = data(select,:);
                            Y1        = randsample(find(Y == s1),ntrialspersub(run));
                            Y2        = randsample(find(Y == s2),ntrialspersub(run));
                            Y         = Y([Y1;Y2]);
                            X         = X([Y1;Y2],:);
                            if randomize == 1
                                Y = Shuffle(Y);
                            end
                            
                            P         = cvpartition(Y,'Holdout',PHoldout); %prepares trainings vs testset
                            cmq     = sprintf('-t 0 -c 1 -q');
                            ind     = logical(P.training);
                            model   = svmtrain(Y(ind), X(ind,:), cmq);
                            ind     = logical(P.test);
                            [~,predicted]               = evalc('svmpredict(Y(ind), X(ind,:), model);');%doing it like this supresses outputs.
                            %                         confmats(:,:,n)   = confusionmat(Y(ind),predicted,'order',[s1 s2]);
                            result(n)         = sum(Y(ind)==predicted)./length(predicted);
                        end
                        performance(s1,s2,run,:) = result;
                        result = [];
                    end
                end
            end
        end
        save(savedfile,'performance','nsub')
    else
        load(savedfile)
        fprintf('Was already saved, loading it from \n%s\n',savedfile)
    end
    
    varargout{1} = performance;
    nsub = size(performance,2);
    keyboard
    av_perf = mean(performance,3);
    av_perf = [av_perf;zeros(1,nsub)];
    perf = av_perf(logical(triu(ones(nsub,nsub))));
    
    average_performance = mean(perf);
    std_average_performance     = std(perf);
    
elseif strcmp(varargin{1},'svm_classify_subs_1vsrest')
    
    
    nbootstrap = 1000;
    randomize  = 0;
    PHoldout    = .2;
    
    
    savepath      = fullfile(path_project,'data\midlevel\SVM_idiosyncrasy\');
    filename      = sprintf('performance_1vsrest_r%d.mat',randomize);
    savedfile     = [savepath filename];
    %   savepath      = fullfile(path_project,['data\midlevel\SVM_idiosyncrasy_' name_analysis '_rand' num2str(r) '.mat']);
    if nargin > 1
        force = varargin{2};
    end
    
    if ~exist(savepath)
        mkdir(savepath)
    end
    
    if ~exist(savedfile) || force==1
        
        if ~exist(savedfile)
            fprintf('Trial data not yet found at %s, \ncollecting it now. \n',savepath)
            WaitSecs(2); %make it readable.
        elseif force == 1
            fprintf('Forced to compute this file anew: \n %s',savepath)
            WaitSecs(2); %make it readable.
        end
        
        [data,labels] = FPSA_FearGen('get_trials_for_svm_idiosync');
        nsub = max(unique(labels.easy_sub));
        performance = [];
        
        for run = 1:4
            fprintf('Starting 1vsrest classification for run = %d.\n',run)
            
            fprintf('Bootstrap count: N = 0... ')
            for n = 1:nbootstrap
                if mod(n,200) == 0
                    fprintf('%03d... ',n)
                end
                
                select    = logical(labels.run == run);
                Y         = labels.easy_sub(select)';
                X         = data(select,:);
                if randomize == 1
                    Y = Shuffle(Y);
                end
                
                P       = cvpartition(Y,'Holdout',PHoldout); %prepares trainings vs testset
                cmq     = sprintf('-t 0 -c 1 -q');
                ind     = logical(P.training);
                model   = ovrtrain(Y(ind), X(ind,:), cmq);
                ind     = logical(P.test);
                [~,predicted]               = evalc('ovrpredict(Y(ind), X(ind,:), model);');%doing it like this supresses outputs.
                
                confmats(:,:,run,n)   = confusionmat(Y(ind),predicted,'order',unique(labels.easy_sub));
                result(n)         = sum(Y(ind)==predicted)./length(predicted);
            end
            fprintf('\n')
            performance(run,:) = result;
        end
        save(savedfile,'confmats','performance','nsub')
    else
        load(savedfile)
        fprintf('Was already saved, loading it from \n%s\n',savedfile)
    end
    
    varargout{1} = performance;
    varargout{2} = confmats;
    
elseif strcmp(varargin{1},'SFig_03')
    %% plot
    savepath = sprintf('%s/data/midlevel/SVM/',path_project);
    files = cellstr(ls(savepath));
    %     neigs = [63 69 74 75 14 19 17 19];
    crit = {'ellbow','ellbow','ellbow','ellbow','var','var','var','var'};
    run = [1:4 1:4];
    for c = 1:8
        expr = sprintf('r0_run%d_crit%s_exclmouth_0.mat',run(c),crit{c});
        findfile = regexp(files,expr,'match');
        ind = find(~cellfun(@isempty,findfile));
        load([savepath files{ind}],'result');
        results(:,:,c) = squeeze(mean(result,2));
    end
    
    M = squeeze(mean(results,2));
    SE = squeeze(std(results,[],2)./sqrt(size(results,2)));
    
    % average the three test runs
    M = cat(2,M(:,1:4),mean(M(:,2:4),2),M(:,5:8),mean(M(:,6:8),2));
    SE = cat(2,SE(:,1:4),mean(SE(:,2:4),2),SE(:,5:8),mean(SE(:,6:8),2));
    
    
    figure(1001)
    clf
    for n = 1:10;subplot(2,5,n);xlim([-170 215]);l=line(xlim,[.5 .5]); hold on;set(l,'Color','k','LineStyle',':');end
    hold on
    labels = {'Baseline' 'Test_1' 'Test_2' 'Test_3' 'Test_M' };
    for n = 1:10
        subplot(2,5,n);
        Project.plot_bar(-135:45:180,M(:,n));
        hold on;
        errorbar(-135:45:180,M(:,n),SE(:,n),'k.','LineWidth',2)
        ylim([.3 .7])
        set(gca,'YTick',.3:.1:.7,'XTick',[0 180],'XTickLabel',{'CS+' 'CS-'},'FontSize',14);
        box off
        axis square
        if ismember(n,[1 6])
            ylabel('Classified as CS+')
        end
        set(gca,'YTick',[.3 .5 .7])
        xlim([-170 215])
        if n<6
            title(labels{n})
        end
    end
    %% partial figure for figure_04
    figure(1002)
    clf
    for n = 1:3;subplot(1,3,n);xlim([-170 215]);l=line(xlim,[.5 .5]); hold on;set(l,'linewidth',1.5,'Color','k','LineStyle',':');end % so the lines are behind bars
    hold on
    labels = {'Baseline' 'Test' 'Test'};
    spc = 0;
    for n = [6 10]
        spc = spc + 1;
        subplot(1,3,spc);
        Project.plot_bar(-135:45:180,M(:,n));
        e = errorbar(-135:45:180,M(:,n),SE(:,n),'k.');
        set(e,'LineWidth',2,'Color','k')
        hold on;
        ylim([.3 .7])
        set(gca,'YTick',.3:.1:.7,'XTick',[0 180],'XTickLabel',{'CS+' 'CS-'},'FontSize',14);
        box off
        axis square
        ylabel('Classified as CS+')
        set(gca,'YTick',[.3 .5 .7])
        xlim([-170 215])
        title(labels{spc})
    end
    subplot(1,3,3)
    clear results
    savepath = sprintf('%s/data/midlevel/SVM/',path_project);
    files = cellstr(ls(savepath));
    %     neigs = [63 69 74 75 14 19 17 19];
    crit = {'var','var','var','var'};
    for c = 1:4
        expr = sprintf('r0_run%d_crit%s_exclmouth_1.mat',c,crit{c});
        findfile = regexp(files,expr,'match');
        ind = find(~cellfun(@isempty,findfile));
        load([savepath files{ind}],'result');
        results(:,:,c) = squeeze(mean(result,2));
    end
    
    M = squeeze(mean(results,2));
    SE = squeeze(std(results,[],2)./sqrt(size(results,2)));
    M = mean(M(:,2:4),2);
    SE = mean(SE(:,2:4),2);
    Project.plot_bar(-135:45:180,M);
    hold on
    e = errorbar(-135:45:180,M,SE,'k.');
    set(e,'LineWidth',2,'Color','k')
    ylim([.3 .7])
    set(gca,'YTick',.3:.1:.7,'XTick',[0 180],'XTickLabel',{'CS+' 'CS-'},'FontSize',14);
    box off
    axis square
    ylabel('Classified as CS+')
    set(gca,'YTick',[.3 .5 .7])
    xlim([-170 215])
    title(labels{spc})
    
    for n = 1:3
        subplot(1,3,n);
        set(gca,'LineWidth',1.5,'FontSize',16)
        
    end
    
elseif strcmp(varargin{1},'figure_01A');
    %% this is the figure of faces, FDMs and dissimilarity matrices
    % get V1 dissimilarity
    p = Project;
    path2v1 = [strrep(p.path_project,'data','stimuli') 'V1\V1_*'];
    dummy           = dir([path2v1 '*.mat']);
    v1files         = [repmat([fileparts(path2v1) filesep],length(dummy),1) vertcat(dummy(:).name)];
    tfiles          = size(v1files,1);
    im0              = [];
    im               = [];
    c =0;
    for i = 1:tfiles
        c = c+1;
        dummy       = load(v1files(i,:));
        im0(:,:,c)   = dummy.v1;
    end
    im = im0  - repmat(mean(im0,3),[1 1 8]); % mean correction
    dissv1 = 1-corr(reshape(im,[400*400,8])); %dissimilarity v1
    
    figure(1)
    subplot(2,1,1)
    dissv1 = CancelDiagonals(dissv1,NaN);
    [d u]   = GetColorMapLimits(dissv1,2.5);
    imagesc(dissv1,[0 1.8]);
    axis square;c = colorbar;
    set(c,'Ytick',caxis)
    set(gca,'fontsize',15);
    axis off
    
    %chose exemplary subject
    %     figure;
    %     subs = FPSA_FearGen('get_subjects');
    %     sim = FPSA_FearGen('get_rsa',1:100);
    %     for n = 1:61;
    %         subplot(8,8,n);
    %         dissmatz = squareform(sim.correlation(n,:));
    %         dissmatz = CancelDiagonals(dissmatz,0);
    %         dissmatz = dissmatz(9:end,9:end);
    %         [d u]   = GetColorMapLimits(dissmatz,2.5);
    %         imagesc(dissmatz,[0 2]);
    %             axis square;
    %     end
    
    
    exemplsub = 45;
    
    
    subplot(2,1,2);
    subs = FPSA_FearGen('get_subjects');
    mat = FPSA_FearGen('get_fpsa',1:100);
    dissmat = squareform(mat.correlation(subs==exemplsub,:)); %get_rsa uses pdist, which is 1-r already, so no 1-squareform necessary.
    dissmat = dissmat(9:end,9:end);
    dissmat = CancelDiagonals(dissmat,NaN);
    %     [d2 u2]   = GetColorMapLimits(dissmat,2.5);
    imagesc(dissmat,[0 1.8]);
    axis square;c = colorbar;
    set(c,'Ytick',caxis)
    set(gca,'fontsize',15);
    axis off;
    
    
    fix = Fixmat(exemplsub,4);
    s   = Subject(exemplsub);
    figure(exemplsub);
    fix.contourplot(4,[.2 .3 .4 .5 .6]);
    colors        = GetFearGenColors;
    
    colors        = circshift(colors(1:8,:),s.csp-4); %need to resort them so that csp is red. (check again)
    
    for n = 1:8
        subplot(3,3,n);
        h = gca;
        x = xlim;
        y = ylim;
        rectangle('position',[x(1) y(1) diff(xlim) diff(ylim)],'edgecolor',colors(n,:),'linewidth',7);
    end
elseif strcmp(varargin{1},'figure_01B');
    %% this is the cartoon figure where hypotheses are presented;
    figure
    set(gcf,'position',[1958         247        1443         740]);
    %few fun definition to write axis labels
    small_texth = @(h) evalc('h = text(4,9,''CS+'');set(h,''HorizontalAlignment'',''center'',''fontsize'',6);h = text(8,9,''CS-'');set(h,''HorizontalAlignment'',''center'',''fontsize'',6);hold on;plot([4 4],[ylim],''k--'',''color'',[0 0 0 .4]);plot([xlim],[4 4],''k--'',''color'',[0 0 0 .4])');
    small_textv = @(h) evalc('h = text(.5,4,''CS+'');set(h,''HorizontalAlignment'',''right'',''fontsize'',6);h = text(.5,8,''CS-'');set(h,''HorizontalAlignment'',''right'',''fontsize'',6);');
    params = {{[.5 .5] 0} {[4.5 4.5] 0} {[4.5 2.5] 0} {[4.5 4.5] 4.5}};
    titles = {sprintf('Bottom-up\nSaliency') sprintf('Arousal') sprintf('Tuned\nExploration') sprintf('Univariate\nGeneralization')};
    width  = 2.3;
    d      = [-.8 -20 -20 -20];
    u      = [ .8  20  20  20];
    %
    spi = {[1 2 13 14] 25 26 [37 38 49 50]};
    for ncol = 1:4
        
        [model w] = getcorrmat(params{ncol}{1},params{ncol}{2},0,1,width);
        model     = 1-corrcov(model);
        
        % ori   = mdscale(1-model,2,'Criterion','strain','start','cmdscale','options',statset('display','final','tolfun',10^-12,'tolx',10^-12));
        colors    = GetFearGenColors;
        colors    = [colors(1:8,:);colors(1:8,:)];
        %
        [y]       = mdscale(model,2,'Criterion',criterion,'start','cmdscale','options',statset('display','final','tolfun',10^-12,'tolx',10^-12));
        if y(4,1) < 0%simply make all red node located at the same position on the figure;
            y = circshift(y,[4,0]);
        end
        % % row 1
        subplot(6,12,spi{1}+(ncol-1)*3);
        plot(y([1:8 1],1),y([1:8 1],2),'.-.','linewidth',2,'color',[0.6 0.6 0.6]);
        hold on;
        for nface = 1:8
            plot(y(nface,1),y(nface,2),'.','color',colors(nface,:),'markersize',50,'markerface',colors(nface,:));
        end
        hold off;
        xlim([-1 1]);
        ylim([-1 1]);
        axis square;axis off
        title(titles{ncol},'fontweight','normal','horizontalalignment','center','verticalalignment','middle','fontsize',15);
        if ncol == 1
            text(-1.7,max(ylim)/2-.4,sprintf('MDS'),'fontsize',12,'rotation',90,'horizontalalignment','center');
        end
        
        % row 2
        if ncol < 4
            subplot(6,12,spi{2}+(ncol-1)*3);
            imagesc(-w(1,:)'*w(1,:),[d(ncol) u(ncol)]);
            axis off;axis square;
            try
                small_texth();small_textv();
            end
            if ncol == 1
                text(-6.5,max(ylim)/2,sprintf('Covariance\nComponents'),'fontsize',12,'rotation',90,'horizontalalignment','center');
            end
            
            subplot(6,12,spi{3}+(ncol-1)*3);
            imagesc(-w(2,:)'*w(2,:),[d(ncol) u(ncol)]);
            axis off;axis square;
            try
                small_texth();
            end
        else
            subplot(6,12,spi{2}+(ncol-1)*3);
            imagesc(-w(1:2,:)'*w(1:2,:),[d(ncol) u(ncol)]);
            axis off;axis square;
            try
                small_texth();small_textv();
            end
            subplot(6,12,spi{3}+(ncol-1)*3);
            imagesc(-w(3,:)'*w(3,:),[d(ncol) u(ncol)]);
            axis off;axis square;
            try
                small_texth();
            end
        end
        
        %last row
        subplot(6,12,spi{4}+(ncol-1)*3);
        %
        % row 3
        
        axis square
        imagesc(model,[.1 2]);
        if ncol == 1
            pos = get(gca,'position');
            hc  = colorbar('eastoutside');
            set(gca,'position',pos);
            try
                hc.Position(3:4) = hc.Position(3:4)./2;
                set(hc,'Ticks',[0.1 2],'TickLabels',[0 2],'box','off');
            end
        end
        %
        %         imagesc(model);colorbar
        axis off;
        axis square;
        h = text(.5,4,'CS+');set(h,'HorizontalAlignment','right','fontsize',10);
        h = text(.5,8,'CS-');set(h,'HorizontalAlignment','right','fontsize',10);
        h = text(.5,2,sprintf('90%c',char(176)));set(h,'HorizontalAlignment','right','fontsize',6);
        h = text(.5,6,sprintf('-90%c',char(176)));set(h,'HorizontalAlignment','right','fontsize',6);
        
        h = text(4,9,'CS+');set(h,'HorizontalAlignment','center','fontsize',10);
        h = text(8,9,'CS-');set(h,'HorizontalAlignment','center','fontsize',10);
        h = text(2,9,sprintf('90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',6);
        h = text(6,9,sprintf('-90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',6);
        
        if ncol == 1
            text(-2.5,max(ylim)/2,sprintf('Theoretical\nSimilarity\Matrices'),'fontsize',12,'rotation',90,'horizontalalignment','center');
        end
        
        %         set(gca,'xtick',[4 8],'xticklabel',{'CS+' 'CS-'},'yticklabel','')
    end
    
    %     SaveFigure(sprintf('~/Dropbox/feargen_lea/manuscript/figures/figure_01B.png'),'-r300');
    
elseif strcmp(varargin{1},'figure_02B')
    %% Presents evidence for learning manipulation based on explicit ratings as well as skin conductance responses.
    force_scr  = 0;
    force_rate = 0;
    
    p                 = Project;
    subs     = FPSA_FearGen('get_subjects');
    scrsubs  = subs(ismember(subs,p.subjects(p.subjects_scr)));
    scrpath           = sprintf('%sdata/midlevel/SCR_N%d.mat',path_project,length(scrsubs));
    %% SCR
    if ~exist(scrpath)||force_scr == 1
        g        = Group(scrsubs);
        out      = g.getSCR(2.5:5.5);
        save(scrpath,'out');
        clear g
    else
        load(scrpath)
    end
    av       = mean(out.y);
    sem      = std(out.y)./sqrt(length(scrsubs));
    %fit baseline to see if there's tuning
    data.y   = out.y(:,1:8);
    data.x   = repmat(-135:45:180,[length(scrsubs) 1])';
    data.ids = NaN;
    base     = Tuning(data);
    base.GroupFit(3);
    %same for test (cond not possible)
    data.y   = out.y(:,19:26);
    data.x   = repmat(-135:45:180,[length(scrsubs) 1]);
    data.ids = NaN;
    test     = Tuning(data);
    test.GroupFit(3);
    params   = test.groupfit.Est;
    %     params(3)= deg2rad(params(3)); %for VonMises Fit
    
    nulltrials = out.y(:,[9 18 27]);
    
    CI = 1.96*std(nulltrials)./sqrt(length(nulltrials)); %2times because in two directions (MEAN plusminus) % this is for plotting nulltrial CI later
    
    %are SCRS CS+ vs CS- sign. different?
    [h,pval,ci,teststat] = ttest(out.y(:,13),out.y(:,17))
    [h,pval,ci,teststat] = ttest(out.y(:,22),out.y(:,26))
    
    % single subject fits for parameter plot
    scr_ampl = nan(length(scrsubs),4);
    sc = 0;
    for sub = scrsubs(:)'
        sc = sc+1;
        for ph = [2 4]
            s = Subject(sub);
            s.get_fit('scr',ph)
            scr_ampl(sc,ph) = s.get_fit('scr',ph).params(1);
        end
    end
    %% plot SCR
    grayshade = [.8 .8 .8];
    figure(1022);clf;
    h = gcf;
    h.Position = [78         319        1372         775];
    subplot(2,4,1);
    pa = patch([-180 225 225 -180],[mean(nulltrials(:,1))-CI(1) mean(nulltrials(:,1))-CI(1) mean(nulltrials(:,1))+CI(1) mean(nulltrials(:,1))+CI(1)],'r','EdgeColor','none');
    set(pa,'FaceAlpha',.9,'FaceColor',grayshade,'EdgeColor','none')
    set(gca,'ytick',[-.5 0 .5],'fontsize',13)
    line([-180 225],repmat(mean(nulltrials(:,1)),[1 2]),'Color',[.5 .5 .5],'LineWidth',1.5)
    hold on;
    Project.plot_bar(-135:45:180,av(1:8));axis square;box off;hold on;
    errorbar(-135:45:180,av(1:8),sem(1:8),'k.','LineWidth',1.5);
    line([-150 195],repmat(mean(av(1:8)),[1 2]),'Color','k','LineWidth',2)
    %     ylim([-1 1]);
    YLIM = .75
    ylim([-YLIM YLIM])
    ylabel(sprintf('SCR\n(z-score)'))
    
    
    subplot(2,4,2);
    pa = patch([-180 225 225 -180],[mean(nulltrials(:,2))-CI(2) mean(nulltrials(:,2))-CI(2) mean(nulltrials(:,2))+CI(2) mean(nulltrials(:,2))+CI(2)],'r','EdgeColor','none');
    set(pa,'FaceAlpha',.9,'FaceColor',grayshade,'EdgeColor','none')
    line([-180 225],repmat(mean(nulltrials(:,2)),[1 2]),'Color',[.5 .5 .5],'LineWidth',1.5)
    hold on;
    Project.plot_bar(-135:45:180,av(10:17));axis square;box off;hold on;
    errorbar(-135:45:180,av(10:17),sem(10:17),'k.','LineWidth',1.5);
    set(gca,'YTick',0:2,'fontsize',13)    
    ylim([-.5 2.1])
%     line([0 180],[max(ylim) max(ylim)]+range(ylim)/10,'Color','k','LineWidth',1.5);
%     text(40,[max(ylim)]+range(ylim)/10+.1,'***','FontSize',20);
    
    
    
    subplot(2,4,3);
    pa = patch([-180 225 225 -180],[mean(nulltrials(:,3))-CI(3) mean(nulltrials(:,3))-CI(3) mean(nulltrials(:,3))+CI(3) mean(nulltrials(:,3))+CI(3)],'r','EdgeColor','none');
    set(pa,'FaceAlpha',.9,'FaceColor',grayshade,'EdgeColor','none')
    line([-180 225],repmat(mean(nulltrials(:,3)),[1 2]),'Color',[.5 .5 .5],'LineWidth',1.5)
    hold on;
    Project.plot_bar(-135:45:180,av(19:26));axis square;box off;hold on;
    errorbar(-135:45:180,av(19:26),sem(19:26),'k.','LineWidth',1.5);
    set(gca,'ytick',[-.5 0 .5],'fontsize',13)
    x = -150:0.1:195;
    plot(test.groupfit.x_HD,test.groupfit.fit_HD,'k-','LineWidth',2)
    %     ylim([-1 1])
    ylim([-YLIM YLIM])
%     line([0 180],[max(ylim) max(ylim)]+range(ylim)/10,'Color','k','LineWidth',1.5);
%     text(40,[max(ylim)]+range(ylim)/10+.1,'***','FontSize',20);
%     axis tight
    %     set(gca,'YTick',[0 1])
    
    %
    %
    differ = (out.y(:,13)-out.y(:,17))./out.y(:,17);
    mean(differ)
    std(differ)
    
    subplot(2,4,4)
    H=plot(linspace(1,2,63),scr_ampl(:,2),'ko',linspace(1,2,63)+2,scr_ampl(:,4),'ko',[1 2],[mean(scr_ampl(:,2)) mean(scr_ampl(:,2))],'k-',[3 4],[mean(scr_ampl(:,4)) mean(scr_ampl(:,4))],'k-');
    xlim([0 5]);ylim([-2.5 4]);
    box off;
    set(gca,'xtick',[1.5 3.5],'XTickLabel',{'Base.','Gen.'},'XTickLabelRotation',0,'fontsize',13,'ytick',[-4 -2 0 2 4]);    
    CCC=100;
    H(1).Color=[CCC CCC CCC]./255;
    H(2).Color=[CCC CCC CCC]./255;
    H(3).LineWidth = 4;
    H(4).LineWidth = 4;;    
%     line([1.5 3.5],[max(ylim) max(ylim)]+range(ylim)/10,'color','k','linewidth',1.5)      
%     text(2.5,[max(ylim)]+range(ylim)/10+.1,'***','fontsize',20,'horizontalalignment','center');
%     text(3.5,3.6,'***','fontsize',20,'horizontalalignment','center');
    axis tight
    xlim([0 5]);
    box off;
    set(gca,'color','none');
    box off;axis square;drawnow;alpha(.5);
    
    ylabel('Amplitude')
    ht = text(2.5,max(scr_ampl(:,4))+.4,'***','FontSize',18);
%     set(ht,'Rotation',0);
    ylim([-3 5])
    box off
    
    %% plot ratings
    subs              = FPSA_FearGen('get_subjects');
    ratepath           = sprintf('%sdata/midlevel/Ratings_N%d.mat',path_project,length(subs));
    %
    if ~exist(ratepath)||force_rate == 1
        g                 = Group(subs);
        ratings           = g.getRatings(2:4);
        save(ratepath,'ratings');
        clear g
    else
        load(ratepath)
    end
    % check all three phases for tuning
    for ph = 1:3
        data.y   = ratings(:,:,ph);
        data.x   = repmat(-135:45:180,[length(subs) 1]);
        data.ids = subs;
        t     = Tuning(data);
        t.GroupFit(3);
        fit(ph) = t.groupfit;
    end
    % get single sub Gaussian ampl
    for ph = 1:3
        sc = 0;
        for sub = subs(:)'
            sc = sc+1;
            s = Subject(sub);
            s.fit_method = 3;
            rate_ampl(sc,ph) = s.get_fit('rating',ph+1).params(1);
        end
    end
    
    %%
    for n = 1:3
        sp = n+4;
        subplot(2,4,sp)
        if n > 1
            l = line([-150 195],repmat(mean(mean(ratings(:,:,1))),[1 2]),'Color','k','LineWidth',2);
            set(l,'LineStyle',':')
        end
        hold on
        Project.plot_bar(-135:45:180,mean(ratings(:,:,n)));
        hold on;
        e        = errorbar(-135:45:180,mean(ratings(:,:,n)),std(ratings(:,:,n))./sqrt(size(ratings,1)),'k.');
        set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',[0 5 10],'FontSize',12)
        set(e,'LineWidth',2,'Color','k')
        hold on;
        if n>1
            
        end
        ylim([0 10])
        xlim([-180 225])
        axis square
        box off
    end
    %
    subplot(2,4,5);ylabel(sprintf('Shock\nExpectancy'))
    hold on;
    line([-150 195],repmat(mean(mean(ratings(:,:,1))),[1 2]),'Color','k','LineWidth',2); %null model in baseline
    subplot(2,4,6);
    plot(fit(2).x_HD,fit(2).fit_HD,'k','LineWidth',2)%add Groupfit line Cond
%     line([0 180],[8 8],'Color','k','LineWidth',1.5)
%     text(30,8.5,'***','FontSize',20)
    subplot(2,4,7);
    plot(fit(3).x_HD,fit(3).fit_HD,'k','LineWidth',2)%add Groupfit line Test
%     line([0 180],[8 8],'Color','k','LineWidth',1.5)
%     text(30,8.5,'***','FontSize',20)
   %
    subplot(2,4,8)    
    H=plot(linspace(1,2,74)+4,rate_ampl(:,3),'ko',linspace(1,2,74),rate_ampl(:,1),'ko',linspace(1,2,74)+2,rate_ampl(:,2),'ko',[1 2],[mean(rate_ampl(:,1)) mean(rate_ampl(:,1))],'k-',[3 4],[mean(rate_ampl(:,2)) mean(rate_ampl(:,2))],'k-',[5 6],[mean(rate_ampl(:,3)) mean(rate_ampl(:,3))],'k-');
    xlim([0 7]);
    %
    ylim([-12 12]);
    box off;
    set(gca,'xtick',[1.5 3.5 5.5],'XTickLabel',{'Base.','Cond.' 'Gen.'},'XTickLabelRotation',0,'fontsize',13,'ytick',[-10 -5 0 5 10],'xgrid','on');    
    CCC=100;
    H(1).Color=[CCC CCC CCC]./255;
    H(2).Color=[CCC CCC CCC]./255;
    H(3).Color=[CCC CCC CCC]./255;
    H(4).LineWidth = 4;
    H(5).LineWidth = 4;;    
    H(6).LineWidth = 4;;    
    %
%     line([1.5 3.5],[max(ylim) max(ylim)]+range(ylim)/10,'color','k','linewidth',1.5)      
%     text(2.5,[max(ylim)]+range(ylim)/10+.1,'***','fontsize',20,'horizontalalignment','center');
%     text(3.5,3.6,'***','fontsize',20,'horizontalalignment','center');        
    box off;
    set(gca,'color','none');
    box off;axis square;drawnow;alpha(.5);
    
    ylabel('Amplitude')
    box off
    
    %%
    keyboard
elseif strcmp(varargin{1},'figure_02B_get_params')
    % get single sub params to plot them in fig 2
    
    method = 3;
    force  = 0;
    
    pathparams = sprintf('%sdata/midlevel/params_fit_method_%d.mat',path_project,method);
    
    if ~exist(pathparams)|| force == 1
        % get single subject fit's parameter - SCR
        subs     = FPSA_FearGen('get_subjects');
        scrsubs = subs(ismember(subs,Project.subjects(Project.subjects_scr)));
        
        scr = [];
        scr.params =nan(length(scrsubs),4,2);
        scr.sub    =nan(length(scrsubs),4);
        scr.pval   =nan(length(scrsubs),4);
        sc = 0;
        for n =  scrsubs(:)'
            sc = sc+1;
            s    = Subject(n);
            s.fit_method = 3;
            for ph = [2 4] % Cond no fit, just 2 datapoints
                fit = s.get_fit('scr',ph);
                scr.params(sc,ph,:)    = fit.params(1);
                scr.pval(sc,ph)        = fit.pval;
                scr.sub(sc,ph)         = n;
            end
        end
        scr.valid    = scr.pval > -log10(.05);%selection criteria
        
        % get single subject fit's parameter - Ratings
        
        subs     = FPSA_FearGen('get_subjects');
        rate = [];
        rate.params =nan(length(subs),4,2);
        rate.sub    =nan(length(subs),4);
        rate.pval   =nan(length(subs),4);
        sc = 0;
        for n =  subs(:)'
            sc = sc+1;
            s    = Subject(n);
            s.fit_method = 3;
            for ph = 2:4
                fit = s.get_fit('rating',ph);
                rate.params(sc,ph,:)    = fit.params(1:2);
                rate.pval(sc,ph)        = fit.pval;
                rate.sub(sc,ph)         = n;
            end
        end
        rate.valid    = rate.pval > -log10(.05);%selection criteria
        save(pathparams,'rate','scr')
    else
        fprintf('Found saved file, loading it.\n')
        load(pathparams)
    end
    
    varargout{1} = scr;
    varargout{2} = rate;
    
    
    %%
    %     plot
    clf
    subplot(2,1,1)
    yyaxis left;boxplot(scr.params(:,:,1),'positions',1:4,'Width',.4,'Color','b','boxstyle','filled');
    set(gca,'FontSize',13,'YColor','b')
    yyaxis right;boxplot(scr.params(:,:,2),'positions',6:9,'Width',.4,'Color','k','boxstyle','filled');
    set(gca,'FontSize',13,'YColor','k','YTick',[0:50:150])
    xlim([0 10])
    line([0 5],[0 0],'Color','k')
    set(gca,'XTick',[2 4 7 9],'XTickLabel',{'B','T','B','T'})
    box off
    
    subplot(2,1,2)
    yyaxis left;boxplot(rate.params(:,:,1),'positions',1:4,'Width',.4,'Color','b','boxstyle','filled');
    set(gca,'FontSize',13,'YColor','b')
    yyaxis right;boxplot(rate.params(:,:,2),'positions',6:9,'Width',.4,'Color','k','boxstyle','filled');
    set(gca,'FontSize',13,'YColor','k','YTick',[0:50:150])
    xlim([0 10])
    set(gca,'XTick',[2 3 4 7 8 9],'XTickLabel',{'B','C','T','B','C','T'})
    line([5.5 10],[0 0],'Color','k')
    box off
    
    
    
elseif strcmp(varargin{1},'SFig_02_tuneduntuned')
    
    figure(1);
    g                 = Group(FPSA_FearGen('get_subjects'));
    ratings           = g.getRatings(2:4);
    g.tunings.rate{2} = Tuning(g.Ratings(2));
    g.tunings.rate{3} = Tuning(g.Ratings(3));
    g.tunings.rate{4} = Tuning(g.Ratings(4));
    
    g.tunings.rate{2}.GroupFit(8);
    g.tunings.rate{3}.GroupFit(8);
    g.tunings.rate{4}.GroupFit(8);
    %%
    f = figure(1022);
    set(f,'position',[0        0        794         481])
    clf
    for n = 2:4
        sn = n-1;
        subpl(n) =  subplot(2,3,sn);
        if n > 2
            l = line([-150 195],repmat(mean(mean(ratings(:,:,2))),[1 2]),'Color','k','LineWidth',2);
            set(l,'LineStyle',':')
        end
        hold on
        Project.plot_bar(-135:45:180,mean(ratings(:,:,sn)));
        %         Project.plot_bar(mean(ratings(:,:,sn)));
        hold on;
        e        = errorbar(-135:45:180,mean(ratings(:,:,sn)),std(ratings(:,:,sn))./sqrt(size(ratings,1)),'k.');
        set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',[0 5 10],'FontSize',12)
        %         SetFearGenBarColors(b)
        set(e,'LineWidth',2,'Color','k')
        ylim([0 10])
        xlim([-180 225])
        axis square
        box off
    end
    %
    subplot(2,3,1);ylabel('Rating of p(shock)','FontSize',12)
    hold on;
    %add Groupfit line
    params = [g.tunings.rate{3}.groupfit.Est; g.tunings.rate{4}.groupfit.Est];
    params = [params(:,1) params(:,2) deg2rad(params(:,3)) params(:,4)];
    x = linspace(-150,195,10000);
    
    subplot(2,3,1);
    line([-150 195],repmat(mean(mean(ratings(:,:,2))),[1 2]),'Color','k','LineWidth',2)
    
    subplot(2,3,2);
    plot(x,VonMises(deg2rad(x),params(1,1),params(1,2),params(1,3),params(1,4)),'k-','LineWidth',2)
    line([0 180],[8 8],'Color','k','LineWidth',1.5)
    text(30,8.5,'***','FontSize',20)
    
    subplot(2,3,3);
    plot(x,VonMises(deg2rad(x),params(2,1),params(2,2),params(2,3),params(2,4)),'k-','LineWidth',2)
    line([0 180],[8 8],'Color','k','LineWidth',1.5)
    text(30,8.5,'***','FontSize',20)
    
    clear g
    [~,exclsubs] = FPSA_FearGen('get_subjects');
    g                 = Group(exclsubs);
    ratings           = g.getRatings(2:4);
    g.tunings.rate{2} = Tuning(g.Ratings(2));
    g.tunings.rate{3} = Tuning(g.Ratings(3));
    g.tunings.rate{4} = Tuning(g.Ratings(4));
    
    g.tunings.rate{2}.GroupFit(8);
    g.tunings.rate{3}.GroupFit(8);
    g.tunings.rate{4}.GroupFit(8);
    
    for n = 2:4
        sn = n-1;
        subpl(n) =  subplot(2,3,sn+3);
        if n > 2
            l = line([-150 195],repmat(mean(mean(ratings(:,:,2))),[1 2]),'Color','k','LineWidth',2);
            set(l,'LineStyle',':')
        end
        hold on
        Project.plot_bar(-135:45:180,mean(ratings(:,:,sn)));
        %         Project.plot_bar(mean(ratings(:,:,sn)));
        hold on;
        e        = errorbar(-135:45:180,mean(ratings(:,:,sn)),std(ratings(:,:,sn))./sqrt(size(ratings,1)),'k.');
        set(gca,'XTick',-135:45:180,'XTickLabel',{'' '' '' 'CS+' '' '' '' 'CS-'},'YTick',[0 5 10],'FontSize',12)
        %         SetFearGenBarColors(b)
        set(e,'LineWidth',2,'Color','k')
        ylim([0 10])
        xlim([-180 225])
        axis square
        box off
    end
    %
    subplot(2,3,4);ylabel('Rating of p(shock)','FontSize',12)
    hold on;
    %add Groupfit line
    params = [g.tunings.rate{3}.groupfit.Est; g.tunings.rate{4}.groupfit.Est];
    params = [params(:,1) params(:,2) deg2rad(params(:,3)) params(:,4)];
    x = linspace(-150,195,10000);
    
    subplot(2,3,4);
    line([-150 195],repmat(mean(mean(ratings(:,:,2))),[1 2]),'Color','k','LineWidth',2)
    
    subplot(2,3,5);
    plot(x,VonMises(deg2rad(x),params(1,1),params(1,2),params(1,3),params(1,4)),'k-','LineWidth',2)
    line([0 180],[8 8],'Color','k','LineWidth',1.5)
    text(30,8.5,'***','FontSize',20)
    
    subplot(2,3,6);
    plot(x,VonMises(deg2rad(x),params(2,1),params(2,2),params(2,3),params(2,4)),'k-','LineWidth',2)
    line([0 180],[8 8],'Color','k','LineWidth',1.5)
    text(30,8.5,'***','FontSize',20)
    
elseif strcmp(varargin{1},'figure_03A')
    %% selected subjects are 44 and 47, 31.
    roi_contours            = 0;
    fixmat                  = FPSA_FearGen('get_fixmat');
    fixmat.kernel_fwhm      = 30;
    c                       = 0;
    nphase = 4;%we plot the generalization phase only.
    maps = [];
    for sub                     = varargin{2};
        %
        c  = 0;
        for ncond = circshift([0 45 90 135 180 -135 -90 -45],[1 3]);
            c    = c+1;
            v{c} = {'subject' sub 'deltacsp' ncond 'phase' nphase};
        end
        %cocktail blank correction
        fixmat.getmaps(v{:});
        dummy = fixmat.maps;
        dummy = dummy - repmat(mean(dummy,3),[1 1 8]);
        maps  = cat(3,maps,dummy);
    end
    
    FPSA_FearGen('plot_fdm',maps,roi_contours);
    %     SaveFigure(sprintf('~/Dropbox/feargen_lea/manuscript/figures/figure_03A_SingleSubjects_%02d_phase_%02d_contour_%02d.png',sub,nphase,roi_contours),'-r300');
    
elseif strcmp(varargin{1},'figure_03B')
    %% Produces the figure with fixation counts on 8 faces at different ROIs.
    %Draw the winning model on these count profiles (e.g. Gaussian or null
    %model).
    force    = 0;
    method   = 3;
    fs       = 12; %fontsize
    counts   = FPSA_FearGen('get_fixation_counts');
    counts   = diff(counts,1,4);
    counts   = counts*100;
    m_counts = nanmean(counts,3);
    s_counts = std(counts,1,3)./sqrt(size(counts,3));
    
    
    path_write = sprintf('%s/data/midlevel/ROI_fixcount_singlesub_fit_%d_N%d.mat',path_project,method,size(counts,3));
    %% fit single subs
    X_fit = [];
    Y_fit = [];
    pval = [];
    params = [];
    if ~exist(path_write) || force==1
        for sub = 1:size(counts,3)
            for nroi = 1:4
                data.y   = squeeze(counts(:,nroi,sub))';
                data.x   = repmat(-135:45:180,1);
                data.ids = sub';
                t        = [];
                t        = Tuning(data);
                t.SingleSubjectFit(method);
                pval(nroi,sub) = 10.^-t.fit_results.pval;
                params(nroi,sub,:) = t.fit_results.params;
                X_fit(nroi,sub,:) = t.fit_results.x_HD;
                if t.fit_results.pval > -log10(.05)
                    Y_fit(nroi,sub,:) = t.fit_results.y_fitted_HD;
                else
                    Y_fit(nroi,sub,:) = repmat(mean(t.y(:)),[1 length(t.fit_results.y_fitted_HD)]);
                end
            end
        end
        save(path_write,'X_fit','Y_fit','params','pval')
    else
        load(path_write)
    end
    %print summary
    fitvalid   = pval <.05;
    nvalidsubs = sum(fitvalid,2);
    realtuned  = (pval<.05).*(abs(params(:,:,3))<45);
    nrealtuned = sum(realtuned,2);
    fprintf('Significant tunings (and no peakshift) found in the following number of subjects:\n')
    fprintf('Left eye:  %02d (%02d) of %02d subjects.\n',nvalidsubs(1),nrealtuned(1),length(pval))
    fprintf('Right eye: %02d (%02d) of %02d subjects.\n',nvalidsubs(2),nrealtuned(2),length(pval))
    fprintf('Nose:      %02d (%02d) of %02d subjects.\n',nvalidsubs(3),nrealtuned(3),length(pval))
    fprintf('Mouth:     %02d (%02d) of %02d subjects.\n',nvalidsubs(4),nrealtuned(4),length(pval))
    
    % find tuned people
    for nroi = 1:4
        tunedsubs{nroi} = find(pval(nroi,:)<.05);
    end
    %     %% plot single fits.
    %     for n = 1:length(pval)
    %         subplot(8,8,n);
    %         Project.plot_bar(-135:45:180,counts(:,2,n));hold on;
    %         plot(squeeze(X_fit(2,n,:)),squeeze(Y_fit(2,n,:)),'k','LineWidth',1.5);
    %         set(gca,'XTick',[],'YTick',ylim);
    %     end
    %     %% bar plots of parameter alpha
    %     figure;
    %     for nroi = 1:4
    %         m_all_alpha     = mean(params(nroi,:,1),2);
    %         sem_all_alpha   = std(params(nroi,:,1),0,2)./sqrt(length(params));
    %         m_pval_alpha   = mean(params(nroi,fitvalid(nroi,:),1),2);
    %         sem_pval_alpha = std(params(nroi,fitvalid(nroi,:),1),0,2)./sqrt(nvalidsubs(nroi));
    %         m_tuned_alpha   = mean(params(nroi,logical(realtuned(nroi,:)),1),2);
    %         sem_tuned_alpha = std(params(nroi,logical(realtuned(nroi,:)),1),0,2)./sqrt(nrealtuned(nroi));
    %         subplot(1,3,1)
    %         bar(nroi,m_all_alpha,'c'); hold on;
    %         errorbar(nroi,m_all_alpha,sem_all_alpha,'k.','LineWidth',1.5)
    %         subplot(1,3,2)
    %         bar(nroi,m_pval_alpha,'b'); hold on;
    %         errorbar(nroi,m_pval_alpha,sem_pval_alpha,'k.','LineWidth',1.5)
    %         subplot(1,3,3)
    %         bar(nroi,m_tuned_alpha,'r'); hold on;
    %         errorbar(nroi,m_tuned_alpha,sem_tuned_alpha,'k.','LineWidth',1.5)
    %     end
    %
    % t-tests
    
    for nroi = 1:4
        [hypo(nroi,1) pval(nroi,1)] = ttest(params(nroi,:,1));
        [hypo(nroi,2) pval(nroi,2)] = ttest(params(nroi,fitvalid(nroi,:),1));
    end
    %
    %     subplot(1,3,1);title('all subs');
    %     ylabel('Alpha VonMises \Delta %(test-base)')
    %     text(4, max(ylim)-.0075, pval2asterix(pval(end,1)),'HorizontalAlignment','center','fontsize',16);
    %     subplot(1,3,2);title('subs fit p<.05');
    %     subplot(1,3,3);title('subs fit p<.05 AND shift < 45');
    %     for n = 1:3
    %         subplot(1,3,n);set(gca,'XTick',1:4,'XTickLabel',{'LeftEye','RightEye','Nose','Mouth'});
    %         box off
    %         axis square
    %     end
    %     EqualizeSubPlotYlim(gcf);
    
    figure;
    bp = boxplot(params(:,:,1)'); %third dimension is params, 1 is Ampl
    set(bp,'LineWidth',1.5)
    set(gca,'XTick',1:4,'XTickLabel',{'LeftEye','RightEye','Nose','Mouth'},'YTick',-20:20:40,'FontSize',14);
    ylabel('Amplitude \Delta% (Test-Base)')
    box off
    %% fit group
    X_fit = [];
    Y_fit = [];
    for nroi = 1:4
        data.y   = squeeze(counts(:,nroi,:))';
        data.x   = repmat(-135:45:180,size(counts,3),1);
        data.ids = [1:size(counts,3)]';
        t        = [];
        t        = Tuning(data);
        t.GroupFit(8);
        X_fit = [X_fit ;t.groupfit.x_HD];
        if t.groupfit.pval > -log10(.05)
            Y_fit = [Y_fit ;t.groupfit.fit_HD];
        else
            Y_fit = [Y_fit ;repmat(mean(t.y(:)),[1 length(t.groupfit.fit_HD)])];
        end
    end
    
    %% plot the 4 count-profiles for whole group
    figure;set(gcf,'position',[2176         705        1099         294]);
    t={'Right Eye', 'Left Eye' 'Nose' 'Mouth'};
    for n = 1:4
        H(n) = subplot(1,4,n);
        Project.plot_bar(-135:45:180,m_counts(:,n,1,1));
        set(gca,'XTickLabel','');
        hh=title(sprintf('%s',t{n}),'fontsize',fs,'fontweight','normal');
        if n == 1
            ylabel(sprintf('\\Delta%%\n(after-before)'));
        end
        hold on;
        errorbar(-135:45:180,m_counts(:,n,1,1),s_counts(:,n,1,1),'ko');
        try
            plot(X_fit(n,:),Y_fit(n,:),'linewidth',4,'color',[0 0 0 .5]);
        catch
            plot(X_fit(n,:),Y_fit(n,:),'linewidth',4,'color',[0 0 0]);
        end
        hold off;
        set(gca,'linewidth',1.2,'YTick',[-8 -4 0 4 8],'fontsize',fs,'YTickLabel',{'-8' '' '0' '' '8'});
        if n ~= 1;
            set(gca,'xticklabel','');
        end
        if n == 2 | n == 4 | n == 3
            set(gca,'yticklabel',[]);
        end
        ylim([-8 8]);
        if n < 3
            h = text(0,1.25,'CS+');set(h,'HorizontalAlignment','center','fontsize',12);
            h = text(180,1.25,'CS-');set(h,'HorizontalAlignment','center','fontsize',12);
            h = text(90,1.25,sprintf('90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',9);
            h = text(-90,1.25,sprintf('-90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',9);
        else
            h = text(0,-1.25,'CS+');set(h,'HorizontalAlignment','center','fontsize',12);
            h = text(180,-1.25,'CS-');set(h,'HorizontalAlignment','center','fontsize',12);
            h = text(90,-1.25,sprintf('90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',9);
            h = text(-90,-1.25,sprintf('-90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',9);
        end
        try
            set(gca,'XAxisLocation','origin')
        end
        set(gca,'XGrid','on','YGrid','off')
    end
    subplotChangeSize(H,.025,.025);
    %
    %     supertitle(sprintf('Subject: %d, Runs: %d - %d', current_subject_pool, runs(1) , runs(end)) ,1);
    %     SaveFigure(sprintf('~/Dropbox/feargen_lea/manuscript/figures/figure_3B_CountTuning_SubjectPool_%d.png',current_subject_pool));
elseif strcmp(varargin{1},'figure_04A')
    %% observed similarity matrices
    clf
    sim     = varargin{2};
    %
    cormatz = squareform(nanmean(sim.correlation));
    %     cormatz = CancelDiagonals(cormatz,NaN);
    [d u]   = GetColorMapLimits(cormatz,2.5);
    labels  = {sprintf('-135%c',char(176)) sprintf('-90%c',char(176)) sprintf('-45%c',char(176)) 'CS+' sprintf('+45%c',char(176)) sprintf('+90%c',char(176)) sprintf('+135%c',char(176)) 'CS-' };
    labels  = {'' sprintf('-90%c',char(176)) '' 'CS+' '' sprintf('+90%c',char(176)) '' 'CS-' };
    d       = .62;
    u       = .85;
    fs      = 12;
    %
    set(gcf,'position',[0 0         995         426]);
    %subplot(9,6,[1 2 3 7 8 9 13 14 15]);
    H(1) = subplot(1,3,1);
    h = imagesc(cormatz(1:8,1:8),[d u]);
    %     set(h,'alphaData',~diag(diag(true(8))));
    axis square;axis off;
    h = text(0,4,'CS+');set(h,'HorizontalAlignment','center','fontsize',12,'rotation',45);
    h = text(0,8,'CS-');set(h,'HorizontalAlignment','center','fontsize',12,'rotation',45);
    h = text(0,1,sprintf('90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',8,'rotation',45);
    h = text(0,6,sprintf('-90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',8,'rotation',45);
    h = text(4,9,'CS+');set(h,'HorizontalAlignment','center','fontsize',12,'rotation',45);
    h = text(8,9,'CS-');set(h,'HorizontalAlignment','center','fontsize',12,'rotation',45);
    h = text(2,9,sprintf('90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',8,'rotation',45);
    h = text(6,9,sprintf('-90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',8,'rotation',45);
    title('Baseline','fontweight','normal','fontsize',15);
    %     subplot(9,6,[4 5 6 10 11 12 16 17 18]);
    H(2) = subplot(1,3,2);
    h=imagesc(cormatz(9:16,9:16),[d u]);
    %     set(h,'alphaData',~diag(diag(true(8))));
    axis square;axis off;
    h = text(4,9,'CS+');set(h,'HorizontalAlignment','center','fontsize',12,'rotation',45);
    h = text(8,9,'CS-');set(h,'HorizontalAlignment','center','fontsize',12,'rotation',45);
    h = text(2,9,sprintf('90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',8,'rotation',45);
    h = text(6,9,sprintf('-90%c',char(176)));set(h,'HorizontalAlignment','center','fontsize',8,'rotation',45);
    title('Generalization','fontweight','normal','fontsize',15);
    %
    [indices] = FPSA_FearGen('CompareB2T_RSA');
    [Y X]     = ind2sub([8 8],indices(:,1));
    Y         = Y - .25;
    X         = X - .45;
    hold on;
    fs        = 12;
    for N = 1:length(X);
        if X(N) > Y(N)
            if indices(N,2) < .05 & indices(N,2) > .01;
                if Y(N) == 3.75;
                    text(X(N),Y(N),'*','fontsize',fs,'color','k');
                else
                    text(X(N),Y(N),'*','fontsize',fs,'color','w');
                end
            elseif indices(N,2) < .01 & indices(N,2) > .005;
                text(X(N),Y(N),'**','fontsize',fs,'color','w');
            elseif indices(N,2) < .005 & indices(N,2) > .001;
                text(X(N),Y(N),'***','fontsize',fs,'color','w');
            end
        end
    end
    pos = get(gca,'position');
    %% colorbar
    h2              = colorbar;
    set(h2,'location','east');
    h2.AxisLocation ='out';
    h2.Box          = 'off';
    h2.TickLength   = 0;
    h2.Ticks        = [.62 .85];
    h2.TickLabels   = {'.6' '.8'};
    %     pos             = [pos(1)+pos(3)-.1 .11 .1 .01];
    pos = get(gca,'position')
    try
        set(h2,'Position',[pos(1)+pos(3)+.004 .268 .01 .25])
    end
    % 	set(h2,'Position',pos)
    
    % plot the similarity to cs+
    %     subplot(9,6,[25:27 19:21])
    H(3) = subplot(1,3,3);
    Y         = FPSA_FearGen('get_mdscale',squareform(mean(sim.correlation)),2);
    y         = reshape(Y,length(Y)/16,16)';
    a =      0;
    rm     = [cos(deg2rad(a)) -sin(deg2rad(a)) ;sin(deg2rad(a)) cos(deg2rad(a)) ];
    y      = (rm*y')';
    colors    = GetFearGenColors;
    colors    = [colors(1:8,:);colors(1:8,:)];
    %
    try
        plot(y([1:8 1],1),y([1:8 1],2),'--','linewidth',3,'color',[0 0 0 .5]);
    catch
        plot(y([1:8 1],1),y([1:8 1],2),'--','linewidth',3,'color',[0 0 0]);
    end
    hold on;
    for nface = 1:8
        try
            scatter(y(nface,1),y(nface,2),500,'markerfacecolor',colors(nface,:),'markeredgecolor',colors(nface,:),'markerfacealpha',.0,'markeredgealpha',.75,'linewidth',2);
            plot(y(nface,1),y(nface,2),'o','markerfacecolor',colors(nface,:),'markeredgecolor',colors(nface,:),'linewidth',2);
        catch
            scatter(y(nface,1),y(nface,2),500,'markerfacecolor',colors(nface,:),'markeredgecolor',colors(nface,:),'linewidth',2);
            plot(y(nface,1),y(nface,2),'o','markerfacecolor',colors(nface,:),'markeredgecolor',colors(nface,:),'linewidth',2);
        end
    end
    box off;axis square;axis tight;axis off
    %
    try
        plot(y([1:8 1]+8,1),y([1:8 1]+8,2),'-','linewidth',3,'color',[0 0 0 1]);
    catch
        plot(y([1:8 1]+8,1),y([1:8 1]+8,2),'-','linewidth',3,'color',[0 0 0]);
    end
    hold on;
    for nface = 1:8
        try
            scatter(y(nface+8,1),y(nface+8,2),500,'markerfacecolor',colors(nface,:),'markeredgecolor',colors(nface,:),'markerfacealpha',1,'markeredgealpha',1);
        catch
            scatter(y(nface+8,1),y(nface+8,2),500,'markerfacecolor',colors(nface,:),'markeredgecolor',colors(nface,:));
        end
    end
    %     text(y(8+4,1)-.02,y(8+4,2)+.065,'CS+','FontWeight','normal','fontsize',12);
    %     text(y(8+8,1)-.08,y(8+8,2)+.07,['CS-'],'FontWeight','normal','fontsize',12);
    xlim([-.4 .4]);
    ylim([-.4 .4]);
    box off;axis square;axis off
    subplotChangeSize(H(3),.025,.025);
    subplotChangeSize(H(3),.025,.025);
    subplotChangeSize(H(3),.025,.025);
    %
    % legend
    plot(-.4+.06,.4,'ko','markersize',12)
    text(-.37+.06,.4,'Baseline','fontsize',12);
    hold on;
    plot(-.4+.06,.34,'ko','markersize',12,'markerfacecolor','k');
    text(-.37+.06,.34,'Generaliz.','fontsize',12)
    hold off;
    %%
    %     subplotChangeSize(H,.025,.025);
    %subplot(9,6,[22 23 24 28 29 30])
    %FPSA_FearGen('model_rsa_singlesubject_plot',1:100);
    %     SaveFigure('~/Dropbox/feargen_lea/manuscript/figures/figure_04A.png','-transparent');
    
    
elseif strcmp(varargin{1},'get_fixation_counts')
    %% Collects fixation counts and reports how they change with conditions on 4 different ROIs before and after learning.
    % these numbers are reported in the manuscript.
    filename       = sprintf('counttuning_runs_%02d_%02d.mat',runs(1),runs(end));
    fixmat         = FPSA_FearGen('get_fixmat');
    fixmat.unitize = 0;
    subjects       = unique(fixmat.subject);
    force          = 0;
    c = 0;
    for ns = subjects(:)'
        fprintf('Counting fixations in subject: %03d.\n',ns);
        c = c+1;
        p = 0;
        for phase = [2 4]
            p          = p + 1;
            path_write = sprintf('%s/data/sub%03d/p%02d/midlevel/%s.mat',path_project,ns,phase,filename);
            %
            if exist(filename) ==0 | force;
                %
                cc = 0;
                v  = [];
                for ncond = fixmat.realcond
                    cc   = cc+1;
                    v{cc} = {'subject' ns 'deltacsp' ncond 'phase' phase};
                end
                %
                fixmat.getmaps(v{:});
                %
                for iii = 1:size(fixmat.maps,3)
                    dummy_counts(iii,:) = fixmat.EyeNoseMouth(fixmat.maps(:,:,iii),1);
                end
                save(path_write,'dummy_counts');
            else
                load(path_write);
            end
            count(:,:,c,p) = dummy_counts;%[faces organs subjects phases]
        end
    end
    varargout{1} = count;
    groups.g1 = repmat([1:8]',[1 5 61 2]);
    groups.g2 = repmat(1:5,[8 1 61 2]);
    groups.g3 = repmat(reshape(1:61,[1 1 61]),[8 5 1 2]);
    groups.g4 = repmat(reshape(1:2,[1 1 1 2]),[8 5 61 1]);
    
    varargout{2} = groups;
    %% Compute fixation density and their change.
    P = mean(mean(mean(count(:,:,:,1:2),4),3));
    fprintf('Fixation density in percentage in Baseline + Generalization:\n')
    fprintf('%25s %3.5g\n','Left Eye:', P(1))
    fprintf('%25s %3.5g\n','Right Eye:',P(2))
    fprintf('%25s %3.5g\n','Nose:',P(3))
    fprintf('%25s %3.5g\n','Eyes+Nose:',sum(P(1:3)))
    fprintf('%25s %3.5g\n','Mouth:',P(4))
    fprintf('%25s %3.5g\n','Other:',P(5))
    %
    P = mean(mean(count(:,:,:,2),3))-mean(mean(count(:,:,:,1),3));
    fprintf('Change in Fixation Density in percentage (Generalization - Baseline):\n');
    fprintf('%25s %3.5g\n','Delta Left Eye:', P(1))
    fprintf('%25s %3.5g\n','Delta Right Eye:',P(2))
    fprintf('%25s %3.5g\n','Delta Nose:',P(3))
    fprintf('%25s %3.5g\n','Delta Eyes+Nose:',sum(P(1:3)))
    fprintf('%25s %3.5g\n','Delta Mouth:',P(4))
    fprintf('%25s %3.5g\n','Delta Other:',P(5))
elseif strcmp(varargin{1},'anova_count_tuning');
    %% will fit a model to the density changes.
    [count groups] = FPSA_FearGen('get_fixation_counts');
    %remove the last ROI
    count(:,5,:,:) = [];
    groups.g1(:,5,:,:) = [];
    groups.g2(:,5,:,:) = [];
    groups.g3(:,5,:,:) = [];
    groups.g4(:,5,:,:) = [];
    
    
    Y = Vectorize(count(:,:,:,1)-count(:,:,:,2));
    
    t = table(Y(:),abs(groups.g1(1:1952)-4)',categorical( groups.g2(1:1952)'),categorical( groups.g3(1:1952)'),'variablenames',{'count' 'faces' 'roi' 'subjects'});
    a = fitlm(t,'count ~ 1 + faces + roi + faces*roi')
    
    
elseif strcmp(varargin{1},'behavior_correlation');
    %% Computes correlation with behavior
    
    b      = FPSA_FearGen('get_behavior');%amplitude of the scr response
    fixmat = FPSA_FearGen('get_fixmat');
    %
    a2 = FPSA_FearGen('fix_counts',fixmat,1,15);
    a  = FPSA_FearGen('beta_counts',fixmat,1,15);
    %%
    try
        b.rating_03_center  = abs(b.rating_03_center);
        b.rating_04_center  = abs(b.rating_04_center);
        b.scr_04_center     = abs(b.scr_04_center);
        %corrected improvement:
        %
        %         b.rating_04_center_improvement  = abs(b.rating_03_center-b.rating_04_center)*-1
        %         b.rating_04_center              = ;
        %         b.scr_04_center                 = abs(b.scr_04_center);
    end
    b.rating_04_sigma_y     = [];c
    b.rating_03_sigma_y     = [];
    b.rating_03_offset      = [];
    b.rating_04_offset      = [];
    b.subject               = [];
    b.scr_04_offset         = [];
    b.scr_04_sigma_y        = [];
    
    %     b.scr_04_logkappa       =  log(b.scr_04_kappa);
    %     b.rating_03_logkappa    =  log(b.rating_03_kappa);
    %     b.rating_04_logkappa    =  log(b.rating_04_kappa);
    
    b.si                    =  b.rating_03_kappa    - b.rating_04_kappa;
    %     b.silog                 =  b.rating_03_logkappa - b.rating_04_logkappa;
    %
    addpath('/home/onat/Documents/Code/Matlab/CircStat/');
    dummy = deg2rad([b.rating_04_center  b.rating_03_center]);
    %     b.sip                    =  ((pi - circ_mean(dummy'))./circ_std(dummy'))';
    b.sip                    =  diff(dummy')';
    %
    vnames                  = sort(b.Properties.VariableNames);
    bb = table();
    for n = vnames
        bb.(n{1})=double(b.(n{1}));
    end
    data                    = [a2(:,:,1) a2(:,:,2) a(:,:,1) a(:,:,2) table2array(bb)];
    %%
    figure;
    imagesc(corrcov(nancov(data)).^2,[-.4 .4])
    hold on;
    plot([6 6]-.5+0,ylim,'r')
    plot([6 6]-.5+5,ylim,'r','linewidth',4)
    plot([6 6]-.5+10,ylim,'r')
    plot([6 6]-.5+15,ylim,'r','linewidth',4)
    plot([6 6]-.5+18,ylim,'r')
    plot([6 6]-.5+21,ylim,'r')
    %     plot([6 6]-.5+19,ylim,'r','linewidth',4)
    %     plot([6 6]-.5+16,ylim,'r')
    %     plot([6 6]-.5+22,ylim,'r')
    %     plot([6 6]-.5+25,ylim,'r','linewidth',4)
    %     plot(xlim,[6 6]-.5+25,'r','linewidth',4)
    %     plot(xlim,[6 6]-.5+22,'r')
    %     plot(xlim,[6 6]-.5+19,'r','linewidth',4)
    %     plot(xlim,[6 6]-.5+16,'r')
    plot(xlim,[6 6]-.5+21,'r');
    plot(xlim,[6 6]-.5+18,'r');
    plot(xlim,[6 6]-.5+15,'r','linewidth',4)
    plot(xlim,[6 6]-.5+10,'r')
    plot(xlim,[6 6]-.5+5,'r','linewidth',4)
    plot(xlim,[6 6]-.5+0,'r')
    hold off
    colorbar;axis square;
    set(gca,'ytick',1:size(data,2),'yticklabel',['eyel' 'eyer' 'nose' 'mouth' 'all' 'eyel' 'eyer' 'nose' 'mouth' 'all' 'eyel' 'eyer' 'nose' 'mouth' 'all' 'eyel' 'eyer' 'nose' 'mouth' 'all' vnames],'ticklabelinterpreter','none');axis square;
    set(gca,'xtick',1:size(data,2),'xticklabel',['eyel' 'eyer' 'nose' 'mouth' 'all' 'eyel' 'eyer' 'nose' 'mouth' 'all' 'eyel' 'eyer' 'nose' 'mouth' 'all' 'eyel' 'eyer' 'nose' 'mouth' 'all' vnames],'ticklabelinterpreter','none','XTickLabelRotation',90);
    %     SaveFigure('~/Dropbox/selim/Office/RSAofFDM/FigureCache/BehavioralCorrelation.png','-transparent')
    elsel
    model = corrcov(getcorrmat([2 2],5,0,1));% - corrcov(getcorrmat([1 1],1,0,1))
    %
elseif strcmp(varargin{1},'get_path_project');
    varargout{1} = path_project;
else
    fprintf('No action with this name %s is present...\n',varargin{1});
    varargout ={};
    return;
end

%% References:
% (1) https://github.com/selimonat/fancycarp.git
% (2) https://github.com/selimonat/globalfunctions.git
% (3) An extensive dataset of eye movements during viewing of complex images
% Nature Scientific Data, 2017
% Niklas Wilming, Selim Onat, Jos? P. Ossand?n, Alper A??k, Tim C.
% Kietzmann, Kai Kaspar, Ricardo R. Gameiro, Alexandra Vormberg & Peter K?nig
% (4) Comparing the similarity and spatial structure of neural
% representations: A pattern-component model
% Diedrichsen Jm et al.
% NeuroImage, 2011
