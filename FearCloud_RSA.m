function [varargout]=FearCloud_RSA(varargin);

%% GET THE FIXATION DATA
path_project      = Project.path_project;
correct           = 1;
condition_borders = {1:8 9:16};
tbootstrap        = 1000;
method            = 'correlation';
block_extract     = @(mat,y,x,z) mat((1:8)+(8*(y-1)),(1:8)+(8*(x-1)),z);
current_subject_pool =0;


if strcmp(varargin{1},'get_subjects');
        
    filename = sprintf('%s/midlevel/subjectpool_%03d.mat',path_project,current_subject_pool);    
    if exist(filename) == 0
        if current_subject_pool == 0;
            subjects = Project.subjects_bdnf;
        elseif current_subject_pool == 1%find tuned people;
            
            fprintf('finding tuned subjects first...\n');
            p=[];sub=[];pval=[];;
            for n = 1:82;
                s = Subject(n);
                p = [p ; s.feargen_rating(4).params];
                pval = [pval ; s.feargen_rating(4).pval];
                sub = [sub;n];
            end
            valid    = (abs(p(:,3)) < 45) & pval > -log10(.05);
            fprintf('Found %03d valid subjects...\n',sum(valid));
            subjects = sub(valid);            
            save(filename,'subjects');
        elseif current_subject_pool == 2;%males
            subjects = Project.subjects_bdnf(Project.gender == 1);
        elseif current_subject_pool == 3;%females
            subjects = Project.subjects_bdnf(Project.gender == 2);
        elseif current_subject_pool == 4
            subjects = Project.subjects_bdnf(Project.BDNF == 1);
        elseif current_subject_pool == 5
            subjects = Project.subjects_bdnf(Project.BDNF == 2);
        end
    else
        load(filename);
    end
    subjects = setdiff(subjects,[13 38]);
    %subject 13's phase02 has no valid eye data, exluding that too.
    %subject 30's correlation matrix is full of nans, will not investigate
    %it further but just exclude.
    varargout{1} = subjects;
elseif strcmp(varargin{1},'get_fixmat');
    %% load the fixation data from the baseline and test phases
    filename = sprintf('%s/midlevel/fixmat_subjectpool_%03d.mat',path_project,current_subject_pool);
    fix = [];
    if exist(filename) == 0  
        subjects = FearCloud_RSA('get_subjects',current_subject_pool)
        fix      = Fixmat(subjects,[2 4]);
        save(filename,'fix');
    else
        load(filename)
    end    
    varargout{1} = fix;
elseif strcmp(varargin{1},'get_behavior')
    fixmat = FearCloud_RSA('get_fixmat');
    p = [];
    for ns = unique(fixmat.subject)
        fprintf('subject:%03d...\n',ns);
        s = Subject(ns);
        dummy = s.feargen(4).params;
        p = [p ; dummy(1) vM2FWHM(dummy(2)) abs(dummy(3))];
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
elseif strcmp(varargin{1},'get_rsa')
    %% COMPUTE THE SIMILARITY MATRIX
    %sim = FearCloud_RSA('get_rsa',1:100)
    fixations = varargin{2};
    filename  = sprintf('%s/midlevel/rsa_all_firstfix_%03d_lastfix_%03d_subjectpool_%03d.mat',path_project,fixations(1),fixations(end),current_subject_pool);    
    %
    if exist(filename) ==0 ;
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
    const         = ones(8);
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
            Y                  = ( 1-squareform(mean(data(:,:,i),3)) );
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
    if nargin == 1
        sim        = FearCloud_RSA('get_rsa',1:100);
        [betas ci] = FearCloud_RSA('get_betas',sim);
    else
        betas = varargin{2};
        ci    = varargin{3};
    end
    %
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
            fprintf('Processing subject %03d\ncache name: %s\n',subject,filename);
            path_write = sprintf('%ssub%03d/p%02d/midlevel/%s.mat',path_project,subject,phase,filename);            
            if exist(fileparts(path_write)) == 0;mkdir(fileparts(path_write));end;%create midlevel folder if not there.
            if exist(path_write) == 0
                % create the query cell
                maps             = FearCloud_RSA('get_fixmap',fixmat,subject,1:100);
                maps             = reshape(maps(:,conds),[500 500 length(conds)]);
                out              = blockproc(maps,[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true,'UseParallel',true);
                save(path_write,'out');
            else
                fprintf('already cached...\n');
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
elseif strcmp(varargin{1},'plot_searchlight')'
    %
   
    %%
    
elseif strcmp(varargin{1},'searchlight_bs')
    
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
    
    fixmat      = varargin{2};
    b1          = varargin{3};
    b2          = varargin{4};
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
    
elseif strcmp(varargin{1},'get_mdscale')    
    %%    
    sim    = varargin{2};
%     simmat = squareform(mean(sim));

Criterion='metricstress' ;
ndimen = 2;
if ndimen == 2
    init   = [cosd(-135:45:180)' sind(-135:45:180)'];
    y      = mdscale(sim,ndimen,'Criterion',Criterion,'start',repmat(init,2,1)+rand(16,2).*.1);
    plot(y([1:8 1],1),y([1:8 1],2),'o-','linewidth',3);
    hold on;
    plot(y([1:8 1]+8,1),y([1:8 1]+8,2),'ro-','linewidth',3);
    hold off;
    for n = 1:16;text(y(n,1),y(n,2),mat2str(mod(n-1,8)+1),'fontsize',25);end
elseif ndimen == 3
    init      = [[cosd(-135:45:180)' sind(-135:45:180)' zeros(8,1)];[cosd(-135:45:180)' sind(-135:45:180)' zeros(8,1)]];
    y      = mdscale(sim,ndimen,'Criterion',Criterion,'start',init+rand(16,3).*.01);
    plot3(y([1:8 1],1),y([1:8 1],2),y([1:8 1],3),'o-','linewidth',3);
    hold on;
    plot3(y([1:8 1]+8,1),y([1:8 1]+8,2),y([1:8 1]+8,3),'ro-','linewidth',3);
    hold off;
    for n = 1:16;text(y(n,1),y(n,2),mat2str(mod(n-1,8)+1),'fontsize',25);end
end

elseif strcmp(varargin{1},'anova')    
    
    y = [cr(:,1,1);cr(:,2,1);cr(:,1,2);cr(:,2,2)];
side = [ones(65,1);ones(65,1)*2;ones(65,1);ones(65,1)*2];
phase =[ones(65,1);ones(65,1);ones(65,1)*2;ones(65,1)*2];

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
    subplot(6,6,[1 2 3 7 8 9 13 14 15]);
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
    
    subplot(6,6,[4 5 6 10 11 12 16 17 18]);
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
%%
    subplot(6,6,19:20)
    X = FearCloud_RSA('get_design_matrix');
    imagesc(squareform(X(:,1)),[-1 1]);axis square;axis off
    title(sprintf('Constant\nSimilarity'))
    subplot(6,6,21:22)
    X = FearCloud_RSA('get_design_matrix');
    imagesc(squareform(X(:,2)),[-1 1]);axis square;axis off
    title(sprintf('Perceptual\nSimilarity'))
    subplot(6,6,23:24)
    X = FearCloud_RSA('get_design_matrix');
    imagesc(squareform(X(:,3)),[-1 1]);axis square;axis off
    title(sprintf('CS+\nSimilarity'))
    %%   
    [betas ci] = FearCloud_RSA('get_betas',sim);
    
    location = {[25 26 ] [27 28 ] [29 30 ]};
    color = {[1 0 0] [.5 0 0];[0 0 1] [0 0 .5];[.8 .8 .8] [.4 .4 .4]}';
    c= -1;
    xticks =[];
    for n = 1:size(betas,1);%betas        
        subplot(6,6,location{n});           
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
    SaveFigure('~/Dropbox/feargen_lea/manuscript/figures/figure03.png','-transparent');
elseif strcmp(varargin{1},'figure04');
    fixmat  = FearCloud_RSA('get_fixmat');    
    M      = FearCloud_RSA('searchlight',fixmat,1,15);    
    M      = squeeze(nanmean(M,4));
    M      = reshape(M,[500 500 6]);
    fs     = 15;
    %%    1st column
    figure;
    set(gcf,'position',[ 2132         528        1579         543]);
    d       = -.1;
    u       = .6;
    G      = make_gaussian2D(51,51,32,32,26,26);
    G      = G./sum(G(:));
    h       = subplot(2,4,1);
    map     = M(:,:,1);
    map     = inpaint_nans(map);
    mapc    = conv2(map,G,'valid');
    mapc    = padarray(mapc,[25 25],NaN);
    mapc    = inpaint_nans(mapc);    
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));    
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);    
    hold on;
    h       = imagesc(mapc,[d u]);
    set(h,'alphaData',Scale(abs(map))*.5+.5);    
    %
    [~,h2]  = contourf(X,Y,mapc,3);
    h2.Fill = 'off';    
    hold off
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','')
    set(h3,'box','off','ticklength',0,'ticks',[d u],'fontsize',fs);
    ylabel('BEFORE','fontsize',15)
    title(sprintf('Constant\nSimilarity'));
    %
    h       = subplot(2,4,5);
    map     = M(:,:,4);
    map     = inpaint_nans(map);
    mapc    = conv2(map,G,'valid');
    mapc    = padarray(mapc,[25 25],NaN);
    mapc    = inpaint_nans(mapc);    
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));    
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);    
    hold on;
    h       = imagesc(mapc,[d u]);
    set(h,'alphaData',Scale(abs(map))*.5+.5);    
    %
    [~,h2]  = contourf(X,Y,mapc,3);
    h2.Fill = 'off';    
    hold off
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','')
    set(h3,'box','off','ticklength',0,'ticks',[d u],'fontsize',fs)
    ylabel('AFTER','fontsize',15)
    %% 2nd column
    G      = make_gaussian2D(51,51,4.5,4.5,26,26);
    G      = G./sum(G(:));
    d       = 0;
    u       = .17;
    tcont   = 6;
    h       = subplot(2,4,2);
%     mask      = conv2(M(:,:,1),G,'same')>0.1;
    map       = M(:,:,2);
%     map(mask) = NaN;    
    map     = inpaint_nans(map);    
    mapc    = conv2(map,G,'valid');
    mapc    = padarray(mapc,[25 25],NaN);
    mapc    = inpaint_nans(mapc);    
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));    
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);    
    hold on;
    h       = imagesc(map,[d u]);
    set(h,'alphaData',Scale(abs(map))*.8+.2);    
    %
    [~,h2]  = contourf(X,Y,mapc,tcont);
    h2.Fill = 'off';    
    hold off
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','')
    set(h3,'box','off','ticklength',0,'ticks',[d u],'fontsize',fs)
    title(sprintf('Perceptual\nSimilarity'));
    %
    h       = subplot(2,4,6);
    map     = M(:,:,5);
    map     = inpaint_nans(map);
    mapc    = conv2(map,G,'valid');
    mapc    = padarray(mapc,[25 25],NaN);
    mapc    = inpaint_nans(mapc);    
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));    
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);    
    hold on;
    h       = imagesc(map,[d u]);
    set(h,'alphaData',Scale(abs(map))*.8+.2);    
    %
    [~,h2]  = contourf(X,Y,mapc,tcont);
    h2.Fill = 'off';    
    hold off
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','')
    set(h3,'box','off','ticklength',0,'ticks',[d u],'fontsize',fs)
     %% 3rd column
    G      = make_gaussian2D(51,51,4.5,4.5,26,26);
    G      = G./sum(G(:));
    d       = 0;
    u       = .17;
    tcont   = 6;
    h       = subplot(2,4,3);
%     mask      = conv2(M(:,:,1),G,'same')>0.1;
    map       = M(:,:,3);
%     map(mask) = NaN;    
    map     = inpaint_nans(map);    
    mapc    = conv2(map,G,'valid');
    mapc    = padarray(mapc,[25 25],NaN);
    mapc    = inpaint_nans(mapc);    
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));    
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);    
    hold on;
    h       = imagesc(map,[d u]);
    set(h,'alphaData',Scale(abs(map))*.8+.2);    
    %
    [~,h2]  = contourf(X,Y,mapc,tcont);
    h2.Fill = 'off';    
    hold off
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','')
    set(h3,'box','off','ticklength',0,'ticks',[0 u],'fontsize',fs)
    title(sprintf('CS+\nSimilarity'));
    %
    h       = subplot(2,4,7);
    map     = M(:,:,6);
    map     = inpaint_nans(map);
    mapc    = conv2(map,G,'valid');
    mapc    = padarray(mapc,[25 25],NaN);
    mapc    = inpaint_nans(mapc);    
    [X Y]   = meshgrid(1:size(map,1),1:size(map,2));    
    %plot the image;
    imagesc(X(1,:),Y(:,1)',fixmat.stimulus);    
    hold on;
    h       = imagesc(map,[d u]);
    set(h,'alphaData',Scale(abs(map))*.8+.2);    
    %
    [~,h2]  = contourf(X,Y,mapc,tcont);
    h2.Fill = 'off';        
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','');
    set(h3,'box','off','ticklength',0,'ticks',[0 u],'fontsize',fs)        
    hold off
    %% 4th column
    h       = subplot(2,4,4);
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
    h       = imagesc(map,[d u]);
    hold off
    set(h,'alphaData',.8);        
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','');
    set(h3,'box','off','ticklength',0,'ticks',[0 u],'fontsize',fs)    
    title(sprintf('Fixation\n probability'))
    %================================================================
    h       = subplot(2,4,8);
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
    h       = imagesc(map,[d u]);    
    set(h,'alphaData',.8);        
    h3=colorbar;axis image;set(gca,'xticklabel','','yticklabel','');
    set(h3,'box','off','ticklength',0,'ticks',[0 u],'fontsize',fs)    
    hold off            
    %%    
    colormap jet
    SaveFigure('~/Dropbox/feargen_lea/manuscript/figures/figure04.png','-transparent');
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
        SaveFigure('~/Dropbox/feargen_lea/manuscript/figures/figure05.png','-transparent');
end
end
