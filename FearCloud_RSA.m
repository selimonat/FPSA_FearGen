function [varargout]=FearCloud_RSA(varargin);

%% GET THE FIXATION DATA
path_project      = Project.path_project;
correct           = 1;
condition_borders = {1:8 9:16};
tbootstrap        = 1000;
method            = 'correlation';
block_extract     = @(mat,y,x,z) mat((1:8)+(8*(y-1)),(1:8)+(8*(x-1)),z);

if strcmp(varargin{1},'get_fixmat');
    %% load the fixation data from the baseline and test phases
    filename = sprintf('%s/midlevel/fixmat.mat',path_project);
    fix = [];
    if exist(filename) == 0
        
        fprintf('finding tuned subjects first...\n');
        p=[];sub=[];pval=[];;
        for n = 1:82;
            s = Subject(n);
            p = [p ; s.feargen(4).params];
            pval = [pval ; s.feargen(4).pval];
            sub = [sub;n];
        end
        valid    = (abs(p(:,3)) < 45) & pval > -log10(.05);
        fprintf('Found %03d valid subjects...\n',sum(valid));
        subjects = setdiff(sub(valid),[13 38]);
        %subject 13's phase02 has no valid eye data, exluding that too.
        %subject 30's correlation matrix is full of nans, will not investigate
        %it further but just exclude.
        fix      = Fixmat(subjects,[2 4]);
        save(filename,'fix');
    else
        load(filename)
    end
    varargout{1} = fix;
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
    fixations = varargin{2};
    filename  = sprintf('%s/midlevel/rsa_all_firstfix_%03d_lastfix_%03d.mat',path_project,fixations(1),fixations(end));    
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
        load(filename)
    end
    varargout{1} = sim;
elseif strcmp(varargin{1},'get_rsa_single')
    %% COMPUTE THE SIMILARITY MATRIX
    if nargin > 1
        correct =varargin{2};
    end
    filename = sprintf('%s/midlevel/rsa_all_single_correction_%d.mat',path_project,correct);
    tfix = 5;
    %
    if exist(filename) ==0 ;
        fix      = FearCloud_RSA('get_fixmat');
        tsub     = length(unique(fix.subject));
        subc     = 0;
        for subject = unique(fix.subject);
            subc = subc + 1;
            %creaete the query cell
            v    = [];
            c    = 0;
            for ph = [2 4]
                for nfix = 1:tfix
                    for cond = -135:45:180
                        c    = c+1;
                        v{c} = {'phase', ph, 'deltacsp' cond 'subject' subject 'fix' nfix};
                    end
                end
            end
            % plot and save fixation maps
            fix.getmaps(v{:});
            maps          = fix.vectorize_maps;
            if correct
                i         = 1:(8*tfix);
                maps(:,i) = demean(maps(:,i)')';
                i         = (8*tfix)+1:size(maps,2);
                maps(:,i) = demean(maps(:,i)')';
            end
            for method = {'euclidean' 'cosine' 'correlation'}
                fprintf('Subject: %03d, Method: %s\n',subject,method{1});
                sim.(method{1})(subc,:)       = pdist(maps',method{1});%both phases simultaneously           
            end
        end
        save(filename,'sim');
    else
        load(filename)
    end
    varargout{1} = sim;        
elseif strcmp(varargin{1},'plot_rsa');
    %% plot correlation matrices without fisher
    figure(1000);
    sim     = varargin{2};
    cormatz = 1-squareform(nanmean(sim.correlation));
    cormatz = CancelDiagonals(cormatz,NaN);        
    [d u] = GetColorMapLimits(cormatz,2.5);
    imagesc(cormatz,[d u]);    
    axis square;colorbar
    set(gca,'fontsize',15)
    axis off
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
            n               = n +1;
            i               = randsample(1:tsub,tsub,1);
            Y               = 1-squareform(nanmean(data(:,:,i),3));
            betas(n,nblock,:)  = X\Y';
        end
    end
    % get errorbars for that
    ci           = prctile(betas,[2.5 97.5]);
    varargout{1} = squeeze(mean(betas));
    varargout{2} = ci;
    %
elseif strcmp(varargin{1},'plot_betas')    
%     betas  = varargin{2};
%     ci     = varargin{3};
%     tbetas = size(betas,2);
    %% plot these betas
%     keyboard
% %     subplot(1,size(betas,2),1);
% %     bar(betas');hold on
% %     errorbar(1:tbetas,betas',betas-ci(1,1,1),ci(2,1,1)-mean(betas(:,1,1)),'k.');
% %     errorbar(2,mean(betas(:,2,1)),mean(betas(:,2,1))-ci(1,2,1),ci(2,2,1)-mean(betas(:,2,1)),'k.');
% %     axis square    
elseif strcmp(varargin{1},'searchlight')
    
    fixmat   = varargin{2}
    b1       = varargin{3};
    b2       = varargin{4};
    %
    tsub     = length(unique(fixmat.subject));    
    fun      = @(block_data) FearCloud_RSA('fun_handle',block_data.data);%what we will do in every block
    subc     = 0;
    for subject = unique(fixmat.subject);
        fprintf('Processing subject %03d\n',subject);
        subc = subc + 1;
        % craete the query cell
        maps      = FearCloud_RSA('get_fixmap',fixmat,subject,1:100);
        maps      = reshape(maps,[500 500 16]);        
        B1(:,:,:,subc,1) = blockproc(maps(:,:,1:8),[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true);
        B1(:,:,:,subc,2) = blockproc(maps(:,:,9:16),[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true);
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
        B1(:,:,:,bs,1) = blockproc(maps(:,:,1:8),[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true);
        B1(:,:,:,bs,2) = blockproc(maps(:,:,9:16),[b1 b1],fun,'BorderSize',[b2 b2],'TrimBorder', false, 'PadPartialBlocks', true);
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
      X            =  FearCloud_RSA('get_design_matrix');      
      betas(1,1,:) = X\Y';      
   else
       betas(1,1,:)= [NaN NaN NaN];
   end
   varargout{1} = betas;
elseif strcmp(varargin{1},'fix_counts')
    
    fixmat         = varargin{2};
    fixmat.unitize = 0;
    subjects = unique(fixmat.subject);
    c = 0;
    for ns = subjects(:)'
        fprintf('Counting fixations in subject: %03d.\n',ns)
        c = c+1;
        p = 0;
        for phase = [2 4]        
            p = p +1;
            fixmat.getmaps({'phase' phase 'subject' ns});
            dummy        = fixmat.maps;
            count(c,:,p) = fixmat.EyeNoseMouth(dummy);            
        end           
    end
    varargout{1} = count;
end
% %     %% loadings for single diagonal elements
% %     n = 0;
% %     tsub = size(corrmat,3);
% %     betasfix = [];
% %     while n < 1000
% %         i       = randsample(1:tsub,tsub,1);
% %         n       = n +1
% %         for db = 1:8%diagonal blocks
% %             Y               = block_extract(fisherz_inverse(nanmean(fisherz(corrmat(:,:,i)),3)),db,db,1);
% %             Y(invalid(:))   = [];
% %             Y               = Y';
% %             betasfix(n,db,:) = X\Y;
% %         end
% %     end
% %     betadiff = betasfix(:,5:8,:);%-betasfix(:,1:4,:);
% %     % get the ci
% %     ci2 = prctile(betadiff,[2.5 97.5]);
% %     
% %     %% test - base as barplots
% %     nfix = 4;
% %     for b = 1:3;
% %         subplot(sp(1),sp(2),9+b);
% %         bar(1:nfix,nanmean(betadiff(:,:,b)),'FaceColor',[.03 .1 .4],'EdgeColor','none')
% %         hold on
% %         errorbar(1,nanmean(betadiff(:,1,b)),nanmean(betadiff(:,1,b))-ci2(1,1,b),ci2(2,1,b)-nanmean(betadiff(:,1,b)),'k.','LineWidth',2);
% %         errorbar(2,nanmean(betadiff(:,2,b)),nanmean(betadiff(:,2,b))-ci2(1,2,b),ci2(2,2,b)-nanmean(betadiff(:,2,b)),'k.','LineWidth',2);
% %         errorbar(3,nanmean(betadiff(:,3,b)),nanmean(betadiff(:,3,b))-ci2(1,3,b),ci2(2,3,b)-nanmean(betadiff(:,3,b)),'k.','LineWidth',2);
% %         errorbar(4,nanmean(betadiff(:,4,b)),nanmean(betadiff(:,4,b))-ci2(1,4,b),ci2(2,4,b)-nanmean(betadiff(:,4,b)),'k.','LineWidth',2);
% %         xlabel('Fix','FontSize',12)
% %         xlim([0 5])
% %         %     ylim([-.2 .2])
% %         box off
% %         axis square
% %         set(gca,'FontSize',12)
% %         
% %     end
% %     subplot(sp(1),sp(2),10);
% %     ylabel('beta [a.u.]');
% %     ylim([-.05 .2]);set(gca,'YTick',[0 .2])
% %     subplot(sp(1),sp(2),11);
% %     ylim([-.05 .15]);set(gca,'YTick',[0 .15])
% %     subplot(sp(1),sp(2),12);
% %     ylim([-.05 .15]);set(gca,'YTick',[0 .15])