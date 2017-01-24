function [varargout]=FearCloud_RSA(varargin);

%% GET THE FIXATION DATA
path_project      = Project.path_project;
correct           = 1;
condition_borders = {1:8 9:16};
tbootstrap        = 1000;
method            = 'correlation';
block_extract     = @(mat,y,x,z) mat((1:8)+(8*(y-1)),(1:8)+(8*(x-1)),z);
current_subject_pool =1;


if strcmp(varargin{1},'get_subjects');
    %% returns the subjects based on the CURRENT_SUBJECT_POOL variable;
    filename = sprintf('%s/midlevel/subjectpool_%03d.mat',path_project,current_subject_pool);
    force = 0;
    if exist(filename) == 0 | force
        if current_subject_pool == 0;
            subjects = Project.subjects_bdnf(Project.subjects_ET);
        elseif current_subject_pool == 1%find tuned people;
            
            fprintf('finding tuned subjects first...\n');
            p=[];sub=[];pval=[];;
            for n = Project.subjects_bdnf(Project.subjects_ET);
                s    = Subject(n);
                p    = [p ; s.get_fit('rating',4).params];
                pval = [pval ; s.get_fit('rating',4).pval];
                sub  = [sub;n];
            end
            valid    = (abs(p(:,3)) < 45) & pval > -log10(.05);
            fprintf('Found %03d valid subjects...\n',sum(valid));
            subjects = sub(valid);
            save(filename,'subjects');
        end
    else
        load(filename);
    end
    varargout{1} = subjects;
elseif strcmp(varargin{1},'get_trialcount')    
    %% goes subject by subject and counts the number of trials in phase VARARGIN{2}
    phase = varargin{2};
    fixmat = FearCloud_RSA('get_fixmat');
    s = 0;
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
    %% load the fixation data from the baseline and test phases
    filename = sprintf('%s/midlevel/fixmat_subjectpool_%03d.mat',path_project,current_subject_pool);
    fix      = [];
    force = 0;
    if exist(filename) == 0 | force
        subjects = FearCloud_RSA('get_subjects',current_subject_pool);
        fix      = Fixmat(subjects,[2 4]);
        save(filename,'fix');
    else
        load(filename)
    end
    varargout{1} = fix;
elseif strcmp(varargin{1},'get_behavior')
    %%
    force    = 1;
    filename = sprintf('%s/midlevel/get_behavior_subjectpool_%02d.mat',path_project,current_subject_pool);
    if exist(filename) == 0 | force
        fixmat   = FearCloud_RSA('get_fixmat');%will return fixmat depending on the current_subject_pool
        subjects = unique(fixmat.subject)';
        % get the SCR from phase 4 and 3. THe phase 3 is just the difference
        % between CS+ and CS?
        p = [];p2 = [];scr_amp_03=[];        
        for ns = subjects(:)'
            fprintf('subject:%03d...\n',ns);
            dummy = Subject(ns).get_fit('scr',4).param_table;
            if ~isempty(dummy)
                p     = [p ; dummy];
            else
                p     = [p ; num2cell(nan(1,size(p,2)))];
            end
            %
            dummy     = Subject(ns).get_scr(3);
            if ~isempty(dummy.y)
                scr_amp_03 = [scr_amp_03;dummy.y_mean(4)-dummy.y_mean(8)];
            else
                scr_amp_03 = [scr_amp_03;NaN];
            end
            %
            p2 = [p2;Subject(ns).get_fit('rating',3).param_table Subject(ns).get_fit('rating',4).param_table];
        end        
        p = [p p2 table(subjects(:),'VariableName',{'subject'}) table(scr_amp_03,'VariableName',{'scr_03_amp'})];
        save(filename,'p');
    else
        load(filename)
    end
    varargout{1} = p;    
elseif  strcmp(varargin{1},'get_fixmap')
    %% load fixation map for subject recorded at both phases for fixations FIX.
    % maps are mean corrected for each phase separately.
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
        fixmat.getmaps(v{:});
        maps = cat(2,maps,demean(fixmat.vectorize_maps')');
    end
    varargout{1} = maps;
elseif  strcmp(varargin{1},'get_fixmap_oddeven')
    %% load fixation map for subject recorded at both phases for fixations FIX.
    % maps are mean corrected for each phase separately.
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
elseif strcmp(varargin{1},'get_rsa')
    %% COMPUTE THE SIMILARITY MATRIX
    %sim = FearCloud_RSA('get_rsa',1:100)
    fixations = varargin{2};
    filename  = sprintf('%s/midlevel/rsa_all_firstfix_%03d_lastfix_%03d_subjectpool_%03d.mat',path_project,fixations(1),fixations(end),current_subject_pool);
    force = 0;
    %
    if exist(filename) ==0 | force;
        fixmat   = FearCloud_RSA('get_fixmat',1);
        subc     = 0;
        for subject = unique(fixmat.subject);
            subc                    = subc + 1;
            maps                    = FearCloud_RSA('get_fixmap',fixmat,subject,fixations);
            fprintf('Subject: %03d, Method: %s\n',subject,method);
            sim.(method)(subc,:)    = pdist(maps',method);%
        end
        save(filename,'sim');
    else
        load(filename);
    end
    varargout{1} = sim;
    
elseif strcmp(varargin{1},'get_rsa2')
    %% same as get_rsa but works using a fixmat as input.    
    force      = 0;
    %%
    window_size    = varargin{2};
    window_overlap = varargin{3};    
    t              = 0:1:(window_size-1);
    start_times    = 0:window_overlap:1500-window_size+1;
    time           = repmat(start_times',1,length(t)) + repmat(t,length(start_times),1);
    %%
    filename  = sprintf('%s/midlevel/rsa2_all_windowsize_%03d_window_overlap_%03d_subjectpool_%03d.mat',path_project,window_size,window_overlap,current_subject_pool);
    %%
    if exist(filename) ==0 | force;
        tc = 0;
        for t = 1:size(time,1)-1
            tc = tc+1            
            fixmat   = FearCloud_RSA('get_fixmat');
            fixmat.UpdateSelection('start',time(t,:),'stop',time(t+1,:));
            fixmat.ApplySelection;
            subc     = 0;
            for subject = unique(fixmat.subject);
                subc                    = subc + 1;
                maps                    = FearCloud_RSA('get_fixmap',fixmat,subject,1:100);
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
    
elseif strcmp(varargin{1},'get_rsa_oddeven')
    %% COMPUTE THE SIMILARITY MATRIX using odd and even trials.
    %sim = FearCloud_RSA('get_rsa',1:100)
    fixations = varargin{2};
    filename  = sprintf('%s/midlevel/rsa_all_oddeven_firstfix_%03d_lastfix_%03d_subjectpool_%03d.mat',path_project,fixations(1),fixations(end),current_subject_pool);
    force     = 0;
    %
    if exist(filename) ==0 | force
        fixmat   = FearCloud_RSA('get_fixmat');
        subc     = 0;
        for subject = setdiff(unique(fixmat.subject),58);
            subc                    = subc + 1;
            maps                    = FearCloud_RSA('get_fixmap_oddeven',fixmat,subject,fixations);
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
        
elseif strcmp(varargin{1},'plot_rsa');
    %% plot correlation matrices without fisher
    figure;
    sim     = varargin{2};
    cormatz = 1-squareform(nanmean(sim.correlation));
    cormatz = CancelDiagonals(cormatz,NaN);
    [d u]   = GetColorMapLimits(cormatz,2.5);
    imagesc(cormatz,[d u]);
    axis square;colorbar
    set(gca,'fontsize',15)
    axis off
    %
elseif strcmp(varargin{1},'get_block')
    %% will get the Yth, Xth block from the RSA.
    sim = varargin{2};
    y   = varargin{3};
    x   = varargin{4};
    r   = [];
    for ns = 1:size(sim.correlation,1)
        dummy = squareform(sim.correlation(ns,:));
        r     = cat(3,r,block_extract(dummy,y,x,1));
    end
    varargout{1} = r;
elseif strcmp(varargin{1},'get_design_matrix');
    %% Linear Model on that with constant, physical similarity, aversive generalization components
    x             = [pi/4:pi/4:2*pi];
    const         = .144*ones(8);
    const         = squareform(CancelDiagonals(const,0));
    % phys = Scale([cos(x') sin(x')]*[cos(x') sin(x')]');
    phys          = [cos(x') sin(x')]*[cos(x') sin(x')]';
    phys          = squareform(CancelDiagonals(phys,0));
    gen           = make_gaussian2D(8,8,2,2,4,4); %90 degrees bc this is approx. the mean fwhm in testphase
    gen           = squareform(CancelDiagonals(gen,0));
    
    X             = [const(:) phys(:) gen(:)];
    X(:,2:3)      = OrthogonolizeTwoVectors(X(:,2:3));
    X(:,2:3)      = [zscore(X(:,2:3))];
    
    varargout{1}  = X;
elseif strcmp(varargin{1},'get_betas')
    %% compute loadings on these
    sim    = varargin{2};
    tsub   = size(sim.correlation,1);
    tblock = length(squareform(sim.correlation(1,:)))/8;
    fprintf('Found %02d blocks\n',tblock);
    betas  = [];
    X      = FearCloud_RSA('get_design_matrix');
    for nblock = 1:tblock
        n      = 0;
        data   = FearCloud_RSA('get_block',sim,nblock,nblock);
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
    
    sim     = FearCloud_RSA('get_rsa',1:100);
    [a b c] = FearCloud_RSA('get_betas_singlesubject',sim);
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
    X      = FearCloud_RSA('get_design_matrix');
    for nblock = 1:tblock
        data   = FearCloud_RSA('get_block',sim,nblock,nblock);
        for n = 1:tsub
            Y                  = ( 1-squareform(mean(data(:,:,n),3)) );
            betas(n,:,nblock)  = X\Y';
        end
    end
    % get errorbars for that
    ci           = std(betas)./sqrt(tsub);
    varargout{1} = squeeze(mean(betas));
    varargout{2} = [mean(betas)-ci/2 ;mean(betas)+ci/2];
    varargout{3} = betas;
    %
elseif strcmp(varargin{1},'plot_betas')
    %%
    if nargin == 2
        fix        = varargin{2};
        sim        = FearCloud_RSA('get_rsa',fix);
        [betas ci] = FearCloud_RSA('get_betas',sim);
    elseif nargin==3
        betas = varargin{2};
        ci    = varargin{3};
    end
    %
    figure
    color = {[1 0 0] [.5 0 0];[0 0 1] [0 0 .5];[.8 .8 .8] [.4 .4 .4]}';
    c= -1;
    xticks =[];
    for n = 1:size(betas,1);%betas
        c=c+1.2;
        for m = 1:size(betas,2)%phases
            c = c+1;
            h=bar(c,betas(n,m),1,'facecolor',color{m,n},'edgecolor',color{m,n});
            hold on;
            errorbar(c,betas(n,m),betas(n,m)-ci(1,n,m),betas(n,m)-ci(2,n,m),'k')
            xticks = [xticks c];
        end
    end
    ylim([-.15 .12]);
    hold off;
    box off
    set(gca,'xtick',xticks,'xticklabel','','color','none','xticklabel',{'before' 'after' 'before' 'after' 'before' 'after' },'XTickLabelRotation',45)
    ylabel('\beta weights')
    xlabel('regressors')
    
elseif strcmp(varargin{1},'searchlight')
%%    
    fixmat   = varargin{2};
    b1       = varargin{3};
    b2       = varargin{4};
    filename = DataHash({fixmat.kernel_fwhm,b1,b2});
    %
    tsub     = length(unique(fixmat.subject));
    fun      = @(block_data) FearCloud_RSA('fun_handle',block_data.data);%what we will do in every block
    phc = 0;
    for phase = [2 4];
        subc  = 0;
        phc   = phc + 1;
        conds = condition_borders{phc};
        for subject = unique(fixmat.subject);
            subc             = subc + 1;
            path_write = sprintf('%ssub%03d/p%02d/midlevel/%s.mat',path_project,subject,phase,filename);
            cprintf([1 0 0],'Processing subject %03d\ncache name: %s\n',subject,path_write);
            if exist(fileparts(path_write)) == 0;mkdir(fileparts(path_write));end;%create midlevel folder if not there.
            if exist(path_write) == 0
                % create the query cell
                maps             = FearCloud_RSA('get_fixmap',fixmat,subject,1:100);
                maps             = reshape(maps(:,conds),[500 500 length(conds)]);
                out              = blockproc(maps,[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true,'UseParallel',true);
                save(path_write,'out');
            else
                cprintf([0 1 0],'Already cached...\n');
                load(path_write);
            end	    
            B1(:,:,:,subc,phc)   = out;
            %
            %                 c = 0;
            %                 for m = 1:2
            %                     for n = 1:3
            %                         c = c +1;
            %                         subplot(2,3,c)
            %                         imagesc(nanmean(B1(:,:,n,:,m),4));colorbar;
            %                         drawnow;
            %                     end
            %                 end
        end
    end
    varargout{1} = B1;
                  
elseif strcmp(varargin{1},'searchlight_stimulus')
    %% applies the search light analysis to the V1 representations.
        
    b1          = varargin{2};
    b2          = varargin{3};
    noise_level = varargin{4};
    filename    = 'stimulus_searchlight';
    path_write  = sprintf('%smidlevel/%s_noiselevel_%02d.mat',path_project,filename,noise_level);
    fun         = @(block_data) FearCloud_RSA('fun_handle',block_data.data);%what we will do in every block    
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
    b1       = varargin{3};
    b2       = varargin{4};
    %
    tsub     = length(unique(fixmat.subject));
    fun      = @(block_data) FearCloud_RSA('fun_handle',block_data.data);%what we will do in every block
    bs       = 0;
    while bs < 1000
        bs                = bs+1;
        fprintf('Processing bs %03d\n',bs);
        % craete the query cell
        subject          = randsample(1:tsub,tsub,1);
        maps             = FearCloud_RSA('get_fixmap',fixmat,subject,1:100);
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
elseif strcmp(varargin{1},'fun_handle')
    %%
    maps = varargin{2};
    maps = reshape(maps,[size(maps,1)*size(maps,2) size(maps,3)]);
    if all(sum(abs(maps)))
        Y            = 1-pdist(maps','correlation');
        X            = FearCloud_RSA('get_design_matrix');
        betas(1,1,:) = X\Y';
    else
        betas(1,1,:)= [NaN NaN NaN];
    end
    varargout{1} = betas;
elseif strcmp(varargin{1},'fix_counts')
%%    
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
    
elseif strcmp(varargin{1},'plot_counts')
    %%
    c = varargin{2};
    t      = {'Before' 'After'};
    figure;
    for np = 1:2;
        subplot(1,2,np);
        violin(c(:,1:2,np),[]);
        set(gca,'xtick',1:2,'xticklabel',{'left' 'right'})
        title(t{np});
        ylabel('Fixation Counts');
    end    
elseif strcmp(varargin{1},'beta_counts')
    %%
    fixmat      = FearCloud_RSA('get_fixmat');
    b1          = 1;
    b2          = 15;
    out         = FearCloud_RSA('searchlight',fixmat,b1,b2);
    for np = 1:2
        for ns = 1:size(out,4);
            for beta = 2%1:size(out,3)
                map            = out(:,:,beta,ns,np);
                count(ns,:,np) = fixmat.EyeNoseMouth(map,0);
            end
        end
    end
    varargout{1} = count;
    cr = count;
    mean((cr(:,1,1)-cr(:,2,1))-(cr(:,1,2)-cr(:,2,2)))
elseif strcmp(varargin{1},'get_mdscale')
    %%
    sim                         = varargin{2};%sim is a valid similarity matrix;    
    ndimen                      = varargin{3};
    Criterion                   ='strain' ;    
    viz                         = 0;    
    [dummy stress disparities]  = mdscale(sim,ndimen,'Criterion',Criterion,'start','cmdscale','options',statset('display','final','tolfun',10^-12,'tolx',10^-12));
    dummy                       = dummy';
    Y                           = dummy(:);                       
    varargout{1} = Y;
    if viz
        FearCloud_RSA('plot_mdscale',Y);
    end
elseif strcmp(varargin{1},'plot_mdscale')
    %%plots the results of the mds
    Y      = varargin{2};
    ndimen = length(Y)./16;    
    y      = reshape(Y,length(Y)/16,16)';%to make it easy plotting put coordinates to different columns;
    if ndimen == 2
        plot(y([1:8 1],1),y([1:8 1],2),'o-','linewidth',3);
        hold on;
        plot(y([1:8 1]+8,1),y([1:8 1]+8,2),'ro-','linewidth',3);
        hold off;
        for n = 1:16;text(y(n,1),y(n,2),mat2str(mod(n-1,8)+1),'fontsize',25);end
    elseif ndimen == 3
        plot3(y([1:8 1],1),y([1:8 1],2),y([1:8 1],3),'o-','linewidth',3);
        hold on;
        plot3(y([1:8 1]+8,1),y([1:8 1]+8,2),y([1:8 1]+8,3),'ro-','linewidth',3);
        hold off;
        for n = 1:16;text(y(n,1),y(n,2),y(n,3),mat2str(mod(n-1,8)+1),'fontsize',25);end
    end
    
elseif strcmp(varargin{1},'get_mdscale_bootstrap')
    %sim = FearCloud_RSA('get_rsa',1:100);
    %FearCloud_RSA('get_mdscale_bootstrap',sim.correlation);    
    %%    
    sim      = varargin{2};%this sim.correlation, not yet squareformed.
    dimen    = varargin{3};
    tsubject = size(sim,1);
    subjects = 1:tsubject;    
    tbs      = 100;
    nbs      = 0;
    y        = nan(16*dimen,tbs);
    while nbs < tbs
        fprintf('Bootstrap: %03d of %03d...\n',nbs,tbs);
        sub          = randsample(subjects,tsubject,1);
        simmat       = squareform(mean(sim(sub,:)));
        y(:,nbs+1)   = FearCloud_RSA('get_mdscale',simmat,dimen);
        nbs          = nbs +1;
    end
    y = zscore(y);
    %% align to the mean separately for each phase
    ya      = [];
    tpoints = size(y,1);
    yc       = permute(reshape(y,[dimen tpoints/dimen tbs]),[2 1 3]);%for alignment we need columns again
    for ph = 1%:2
        i      = [1:16]+8*(ph-1);%indices of the faces for before and after
        E_mean = mean(yc(i,:,:),3);%average of this phase
        for ns = 1:tbs;
            [d z transform] = procrustes(E_mean' , yc(i,:,ns)' , 'Reflection',false);
            ya(i,:,ns)      = z';
        end
    end    
    %%
    FearCloud_RSA('plot_mdscale',Vectorize(mean(ya,3)'));
    %%    
%     yam = mean(ya,3);
%     plot(yam([1:8 1],1),yam([1:8 1],2),'o-','linewidth',3);
%     hold on;
%     plot(yam([1:8 1]+8,1),yam([1:8 1]+8,2),'ro-','linewidth',3);    
%     for node = 1:16;        
%         hold on
%         text(yam(node,1),yam(node,2),mat2str(mod(node-1,8)+1),'fontsize',25);        
%         error_ellipse(squeeze([ya(node,1,:);ya(node,2,:)])','color','k','linewidth',1);        
%     end        
%     axis square  
%     varargout{1} = y;
    %%
elseif strcmp(varargin{1},'anova')
    
    fixmat   = FearCloud_RSA('get_fixmat');    
    cr       = FearCloud_RSA('beta_counts',fixmat,1,15);
    tsubject = length(unique(fixmat.subject))
    
    Y        = [cr(:,1,1) cr(:,2,1) cr(:,1,2) cr(:,2,2)];    
    figure;
    errorbar(mean(Y),std(Y)./sqrt(size(Y,1)));
    y        = [cr(:,1,1);cr(:,2,1);cr(:,1,2);cr(:,2,2)];
    side     = [ones(tsubject,1);ones(tsubject,1)*2;ones(tsubject,1);ones(tsubject,1)*2];
    phase    = [ones(tsubject,1);ones(tsubject,1);ones(tsubject,1)*2;ones(tsubject,1)*2];
    anovan(y,{side(:) phase(:)},'model','full')

elseif strcmp(varargin{1},'figure02')
    p = Project;
    figure(1);
    g = Group(FearCloud_RSA('get_subjects'));
    ratings = g.getRatings(2:4);
    g.tunings.rate{2} = Tuning(g.Ratings(2));
    g.tunings.rate{3} = Tuning(g.Ratings(3));
    g.tunings.rate{4} = Tuning(g.Ratings(4));
    g.tunings.rate{2}.GroupFit(8);
    g.tunings.rate{3}.GroupFit(8);
    g.tunings.rate{4}.GroupFit(8);
    
    for n = 2:4
        sn = n-1;
        subpl(n) =  subplot(2,3,sn);
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
    subplot(2,3,1);ylabel('Rating of p(shock)','FontSize',12)
    hold on;
    %add Groupfit line
    params = [g.tunings.rate{3}.groupfit.Est; g.tunings.rate{4}.groupfit.Est];
    params = [params(:,1) params(:,2) deg2rad(params(:,3)) params(:,4)];
    x = -150:0.1:195;
    subplot(2,3,1);
    line([-150 195],repmat(mean(mean(ratings(:,:,2))),[1 2]),'Color','k','LineWidth',2)
    subplot(2,3,2);
    plot(x,VonMises(deg2rad(x),params(1,1),params(1,2),params(1,3),params(1,4)),'k-','LineWidth',2)
    line([0 180],[9 9],'Color','k','LineWidth',2)
    text(45,9,'***','FontSize',20)
    subplot(2,3,3);
    plot(x,VonMises(deg2rad(x),params(2,1),params(2,2),params(2,3),params(2,4)),'k-','LineWidth',2)
    line([0 180],[9 9],'Color','k','LineWidth',2)
    text(45,9,'***','FontSize',20)
    subplot(2,3,1);title('Baseline','FontSize',14);
    subplot(2,3,2);title('Conditioning','FontSize',14);
    subplot(2,3,3);title('Generalization','FontSize',14);
    
    % SCR
    g = Group(p.subjects_bdnf(p.subjects_scr));
    out = g.getSCR(2.5:5.5);
    av = mean(out.y);
    sem = std(out.y)./sqrt(length(g.ids));
    %fit baseline to see if there's tuning
    data.y = out.y(:,1:8);
    data.x = repmat(-135:45:180,[68 1])';
    data.ids = NaN;
    base = Tuning(data);
    base.GroupFit(8);
    %same for test (cond not possible)
    data.y = out.y(:,17:24);
    data.x = repmat(-135:45:180,[68 1]);
    data.ids = NaN;
    test = Tuning(data);
    test.GroupFit(8);
    params = test.groupfit.Est;
    params(3) = deg2rad(params(3));
    figure(1);
    % plot SCR
    subplot(2,3,4);
    ylabel('SCR (z-score)')
    b = bar(-135:45:180,av(1:8));SetFearGenBarColors(b);axis square;box off;hold on;
    errorbar(-135:45:180,av(1:8),sem(1:8),'k.','LineWidth',1.5);
    line([-150 195],repmat(mean(av(1:8)),[1 2]),'Color','k','LineWidth',2)
    ylim([-.55 .55])
    subplot(2,3,5);
    b = bar(-135:45:180,av(9:16));SetFearGenBarColors(b);axis square;box off;hold on;
    errorbar(-135:45:180,av(9:16),sem(9:16),'k.','LineWidth',1.5);
    subplot(2,3,6);
    b = bar(-135:45:180,av(17:24));SetFearGenBarColors(b);axis square;box off;hold on;
    errorbar(-135:45:180,av(17:24),sem(17:24),'k.','LineWidth',1.5);
    x = -150:0.1:195;
    plot(x,VonMises(deg2rad(x),params(1),params(2),params(3),params(4)),'k-','LineWidth',2)
    ylim([-.55 .55])
    for n = 4:6
        subplot(2,3,n)
        set(gca,'XTick',[0 180],'XTickLabel',{'CS+' 'CS-'},'FontSize',12)
        xlim([-180 225])
    end
    
elseif strcmp(varargin{1},'figure03')
    
    sim     = varargin{2};
    cormatz = 1-squareform(nanmean(sim.correlation));
    cormatz = CancelDiagonals(cormatz,NaN);
    [d u]   = GetColorMapLimits(cormatz,2.5);
    labels = {sprintf('-135%c',char(176)) sprintf('-90%c',char(176)) sprintf('-45%c',char(176)) 'CS+' sprintf('+45%c',char(176)) sprintf('+90%c',char(176)) sprintf('+135%c',char(176)) 'CS-' };
    labels = {'' sprintf('-90%c',char(176)) '' 'CS+' '' sprintf('+90%c',char(176)) '' 'CS-' };
    d = -.3;u = .15;
    fs = 12;
    figure;
    set(gcf,'position',[2132          23         600        1048]);
    subplot(9,6,[1 2 3 7 8 9 13 14 15]);
    h = imagesc(cormatz(1:8,1:8),[d u]);
    %     contourf(CancelDiagonals(cormatz(1:8,1:8),mean(diag(cormatz(1:8,1:8),-1))),4);
    axis square;
    if ~ismac
        set(gca,'fontsize',fs,'xtick',1:8,'ytick',1:8,'XTickLabelRotation',45,'xticklabels',labels,'fontsize',fs,'YTickLabelRotation',45,'yticklabels',labels)
    else
        set(gca,'fontsize',fs,'xtick',1:8,'ytick',1:8,'xticklabels',labels,'fontsize',fs,'yticklabels',labels)
    end
    %     set(h,'alphaData',~diag(ones(1,8)));
    title('Before');
    
    subplot(9,6,[4 5 6 10 11 12 16 17 18]);
    h=imagesc(cormatz(9:16,9:16),[d u]);
    %     contourf(CancelDiagonals(cormatz(9:16,9:16),mean(diag(cormatz(9:16,9:16),-1))),4);
    axis square;h2 = colorbar;set(h2,'location','east');h2.Position = [.91 .65 0.02 .1];h2.AxisLocation='out'
    if ~ismac
        set(gca,'fontsize',fs,'xtick',1:8,'XTickLabelRotation',45,'xticklabels',labels,'fontsize',fs,'YTickLabelRotation',45,'yticklabels',{''})
    else
        set(gca,'fontsize',fs,'xtick',1:8,'xticklabels',labels,'fontsize',fs,'yticklabels',{''})
    end
    title('After')
    if ~ismac
        set(h2,'box','off','ticklength',0,'ticks',[d 0 u],'fontsize',fs)
    end
    %axis off
    %     set(h,'alphaData',~diag(ones(1,8)));
    %% plot the similarity to cs+
    subplot(9,6,[25:27 19:21])    
    Y       = FearCloud_RSA('get_mdscale',squareform(mean(sim.correlation)),2);
    y      = reshape(Y,length(Y)/16,16)';
    plot(y([1:8 1],1),y([1:8 1],2),'o-','linewidth',3);
    label={'' '' '' 'cs+' '' '' '' 'cs-'};
    for n = 1:8;text(y(n,1),y(n,2),label{n},'fontsize',15);end    
    box off;axis square;xlim([-.6 .6]);ylim([-.6 .6]);
    set(gca,'XTickLabel',[],'YTickLabel',[],'xtick',0,'YTick',0);grid on;
    hold off
    %
    subplot(9,6,[28:30 22:24])    
    plot(y([1:8 1]+8,1),y([1:8 1]+8,2),'ro-','linewidth',3);
    label={'' '' '' 'cs+' '' '' '' 'cs-'};
    for n = [1:8]+8;text(y(n,1),y(n,2),label{n-8},'fontsize',15);end    
    box off;axis square;xlim([-.6 .6]);ylim([-.6 .6]);    
    set(gca,'XTickLabel',[],'YTickLabel',[],'xtick',0,'YTick',0);grid on;
    hold on
    plot(y([1:8 1],1),y([1:8 1],2),'o-.','linewidth',2,'color',[114 189 255]/255);
    hold off
    
    %%
    subplot(9,6,[31:32 37:38])
    X = FearCloud_RSA('get_design_matrix');
    imagesc(squareform(X(:,1)),[-1 1]);axis square;axis off
    title(sprintf('Constant\nSimilarity'),'fontsize',8)
    subplot(9,6,[33 34 39:40])
    X = FearCloud_RSA('get_design_matrix');
    imagesc(squareform(X(:,2)),[-1 1]);axis square;axis off
    title(sprintf('Perceptual\nSimilarity'),'fontsize',8)
    subplot(9,6,[35:36 41:42])
    X = FearCloud_RSA('get_design_matrix');
    imagesc(squareform(X(:,3)),[-1 1]);axis square;axis off
    title(sprintf('CS+\nSimilarity'),'fontsize',8)
    %%
    [betas ci] = FearCloud_RSA('get_betas',sim);
    
    location = {[43 44 49 50] [45 46 51 52] [47 48 53 54]};
    color = {[1 0 0] [.5 0 0];[0 0 1] [0 0 .5];[.8 .8 .8] [.4 .4 .4]}';
    c= -1;
    xticks =[];
    for n = 1:size(betas,1);%betas
        subplot(9,6,location{n});
        for m = 1:size(betas,2)%phases
            h=bar(m,betas(n,m),1,'facecolor',color{m,n},'edgecolor',color{m,n});
            hold on;
            errorbar(m,betas(n,m),betas(n,m)-ci(1,n,m),betas(n,m)-ci(2,n,m),'k')
            box off;
            if n ==2
                plot([1 2],[.14 .14],'k-');
                plot([1.5],[.15],'k*');
            elseif n == 3
                plot([1 2],[.035 .035],'k-');
                plot([1.5],[.04],'k*');
            end
        end
        xlim([0 3])
        hold off;
        if ~ismac
            set(gca,'xtick',[1 2],'xticklabel','','color','none','xticklabel',{'before' 'after' 'before' 'after' 'before' 'after' },'XTickLabelRotation',45);
        else
            set(gca,'xtick',[1 2],'xticklabel','','color','none','xticklabel',{'before' 'after' 'before' 'after' 'before' 'after' });
        end
        SetTickNumber(gca,3,'y');
        axis square
        if n == 1
            ylabel('\beta');
        end
    end
    SaveFigure(sprintf('~/Dropbox/feargen_lea/manuscript/figures/figure03_subjectpool_%02d.png',current_subject_pool));
elseif strcmp(varargin{1},'figure04');
        
    fixmat      = FearCloud_RSA('get_fixmat');
    M           = FearCloud_RSA('searchlight',fixmat,1,15);
    Mori        = squeeze(nanmean(M,4));
    %%
    M           = Mori;
    crop_amount = [0 0];
    M           = M( 1+crop_amount(1):end-crop_amount(1), 1+crop_amount(2):end-crop_amount(2),:,:,:,:,:);
%     M           = padarray(M,[crop_amount],'replicate','both');
    M           = reshape(M,[size(M,1) size(M,2) 6]);
    fs          = 15;%font size
    stim        = fixmat.stimulus;
    stim        = stim( 1+crop_amount(1):end-crop_amount(1), 1+crop_amount(2):end-crop_amount(2),:,:,:,:,:);
    %
    G           = make_gaussian2D(51,51,45,45,26,26);
    G           = G./sum(G(:));
    mask        = conv2(M(:,:,4),G,'same');
    
    k = .22;
    M(:,:,4) = conv2(M(:,:,4),G,'same');
    M(:,:,1) = conv2(M(:,:,1),G,'same');
    %    
    figure(12122);clf
    set(gcf,'position',[1952         361        1743         714]);
    spi = [1 2 3 5 6 7];
    cmap = repmat([0 .5 ; 0 .15; 0 0.1],2,1);
    [X Y]   = meshgrid(1:size(M,1),1:size(M,2));    
    for n = 1:6
        h(n)       = subplot(2,4,spi(n));
        imagesc(stim,cmap(n,:));
        hold on;
%         imagesc(M(:,:,n).*r);
        if n == 1 | n == 4
           [dummy,h2] = contourf(M(:,:,n),[-.3 -.3+k -.3+k*2 -.3+k*3 -.3+k*4 -.3+k*5],'color','none');axis ij;                            
        else
           [~,h2]  = contourf(M(:,:,n).*(mask <= k),20,'color','none');axis ij;                                          
        end        
        
        hold on           
        [~,h3]     = contourf(mask,[k k ]);axis ij;                    
        h3.Fill    = 'off';
        hold off           
        %    
        h4=colorbar;set(gca,'xticklabel','','yticklabel','')
        set(h4,'box','off','ticklength',0,'ticks',cmap(n,:),'fontsize',fs);      
        hold off;
        drawnow;        
        contourf_transparency(h2,0.5);;           
        colormap jet 
%         axis off;
                
        if spi(n) == 1
            ylabel('BEFORE','fontsize',15);
        elseif spi(n) == 4;
            ylabel('AFTER','fontsize',15);
        end
        %
        if spi(n) == 1
            title(sprintf('Constant\nSimilarity'));
        elseif spi(n) == 2
            title(sprintf('Perceptual\nSimilarity'));            
        elseif spi(n) == 3
            title(sprintf('CS+\nSimilarity'));
%         elseif n == 5
%             title(sprintf('Fixation\n probability'))
        end
        axis square;        
    end
    %%   
    %% 4th column
    h(7)       = subplot(2,4,4);    
    v = [];
    c = 0;
    %     fixmat.unitize = 0;
    for sub = unique(fixmat.subject)
        c    = c+1;
        v{c} = {'subject' sub 'deltacsp' fixmat.realcond 'phase' 2};
    end
    fixmat.getmaps(v{:});
    map     = nanmean(fixmat.maps,3);
    [d u]   = GetColorMapLimits(map,7);
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);
    hold on;
    hi       = imagesc(map,[0 u]);
    hold off
    set(hi,'alphaData',.5);
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','');
    set(h3,'box','off','ticklength',0,'ticks',[0 u],'fontsize',fs)    
    %================================================================
    h(8)       = subplot(2,4,8);
    v = [];
    c = 0;
    for sub = unique(fixmat.subject)
        c    = c+1;
        v{c} = {'subject' sub 'deltacsp' fixmat.realcond 'phase' 4};
    end
    fixmat.getmaps(v{:});
    map     = mean(fixmat.maps,3);
    [d u] = GetColorMapLimits(map,7);
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);
    hold on;
    hi       = imagesc(map,[0 u]);
    set(hi,'alphaData',.5);
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','');
    end
    set(h3,'box','off','ticklength',0,'ticks',[0 u],'fontsize',fs)
    hold off    
    %%    
    for n = 1:8
        subplotChangeSize(h(n),.05,.05)
        axis square;
    end
    colormap jet
    SaveFigure(sprintf('~/Dropbox/feargen_lea/manuscript/figures/figure04_ver2_subjectpool_%02d.png',current_subject_pool));
elseif strcmp(varargin{1},'figure05');
    %%
    figure;
    fixmat = FearCloud_RSA('get_fixmat');
    c = FearCloud_RSA('fix_counts',fixmat);
    subplot(1,2,1);
    c = cat(2,c(:,1:2,1),c(:,1:2,2));
    c = c(:,[1 3 2 4]);
    bar(mean(c),1,'k');
    hold on;
    errorbar(mean(c),std(c)./sqrt(65),'ro');
    hold off;
    title(sprintf('Fixation\nCount'));
    lab = @() set(gca,'xticklabel',{'left-before' 'left-after' 'right-before' 'right-after'  },'XTickLabelRotation',45,'box','off');
    lab();
    SetTickNumber(gca,3,'y');
    %=========================
    subplot(1,2,2);
    fixmat = FearCloud_RSA('get_fixmat');
    c = FearCloud_RSA('beta_counts',fixmat,1,15);
    c = cat(2,c(:,1:2,1),c(:,1:2,2));
    c = c(:,[1 3 2 4]);
    bar(mean(c),1,'k');
    hold on;
    errorbar(mean(c),std(c)./sqrt(65),'ro');
    title(sprintf('Physical\nSimilarity'));
    lab();
    SetTickNumber(gca,3,'y');
    plot([3 4],[.0112 .0112],'k-')
    plot([3.5],[.0114],'k*')
    hold off;
    SaveFigure(sprintf('~/Dropbox/feargen_lea/manuscript/figures/figure05_subjectpool_%02d.png',current_subject_pool));
elseif strcmp(varargin{1},'figure_single_subject')
    fs = 15;
    fixmat                  = FearCloud_RSA('get_fixmat');
    fixmat.kernel_fwhm      = 25;
    sub                     = varargin{2};    
    %
    c           = 0;    
    for nphase = [4]
        for ncond = [-45 0 45 90 135 180 -135 -90];
            c    = c+1;
            v{c} = {'subject' sub 'deltacsp' ncond 'phase' nphase};
        end
    end
    fixmat.getmaps(v{:});
    %fixmat.maps = cat(3,fixmat.maps(:,:,1:8)-repmat(mean(fixmat.maps(:,:,1:8),3),[1 1 8]),fixmat.maps(:,:,9:end)-repmat(mean(fixmat.maps(:,:,9:end),3),[1 1 8]));    
    fixmat.maps = fixmat.maps(:,:,1:8)-repmat(mean(fixmat.maps(:,:,1:8),3),[1 1 8]);
    % 
    %[d u] = GetColorMapLimits(fixmat.maps(:),1);
    grids  = linspace(min(fixmat.maps(:)),max(fixmat.maps(:)),15);    
    t      = {'-45' 'CS+' '+45' '90' '135' 'CS-' '-135' '-90'};
    figure;set(gcf,'position',[1952         361        1743         714]);
    for n = 1:8
        hhhh(n)=subplot(2,4,n);
        imagesc(fixmat.stimulus);
        hold on
        [~,h2] = contourf(fixmat.maps(:,:,n),grids,'color','none');
        caxis([grids(4) grids(end-6)]);
        drawnow;
        contourf_transparency(h2,0.4);;                   
        %if n == 8
            %h4 = colorbar;        
            %set(h4,'box','off','ticklength',0,'ticks',[[grids(4) grids(end-4)]],'fontsize',fs);      
        %end
        hold off
        axis image;        
        axis off;
        title(t{n});        
    end    
%     colormap jet
    for n = 1:8;
        subplotChangeSize(hhhh(n),.1,.1);
    end
    SaveFigure(sprintf('~/Dropbox/feargen_lea/manuscript/figures/SingleSubjects_%02d.png',sub));
elseif strcmp(varargin{1},'behavior_correlation');
        
    
    b      = FearCloud_RSA('get_behavior');%amplitude of the scr response
    fixmat = FearCloud_RSA('get_fixmat');
    %
    a2 = FearCloud_RSA('fix_counts',fixmat,1,15);
    a  = FearCloud_RSA('beta_counts',fixmat,1,15);    
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
    b.rating_04_sigma_y     = [];
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
    SaveFigure('~/Dropbox/selim/Office/RSAofFDM/FigureCache/BehavioralCorrelation.png','-transparent')
end
