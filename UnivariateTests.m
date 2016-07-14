%% 
% this script provides results that could be used to complement the
% super-complex hyperplan analysis by doing basically down to earth
% comparisons of amount of fixation counts. The question it tries to answer
% is whether and where there are more fixations on the face in CS+ vs. CS?
% conditions. The simple comparison for 0 vs. 180 gives nice results which
% are in line with the hyperplane analysis. However one needs to have a
% control and here also the best control is to check whether 90 vs. -90
% gives similar results. The idea would be that if these differences are
% induces purely by differences of the images, it should also be present in
% any combination of opposite faces.

fix             = Fixmat(setdiff(Project.subjects_1500,7),[2 4]);
p = Project;
subjects        = intersect(find(p.getMask('ET_feargen')),p.subjects_1500)';%22 has to be kicked for svm
%% similarity of fixation maps
cmat  =[];
pval  = [];
sub_c = 0;
v = [];
correct  = 1;
for ns = setdiff(Project.subjects_1500,[20 22 7]);
    ns
    sub_c = sub_c + 1;cond_c = 0;
    for phase = [2 4]
        for ncond = -135:45:180
            cond_c = cond_c + 1;
            v{cond_c} = {'subject' ns 'phase' phase 'deltacsp' [ncond]};
        end
    end
    fix.getmaps(v{:});
    if correct == 1
        %correct phase specific cocktail blank
        fix.maps(:,:,1:8) = fix.maps(:,:,1:8) - repmat(mean(fix.maps(:,:,1:8),3),[1 1 8]);
        fix.maps(:,:,9:end) = fix.maps(:,:,9:end) - repmat(mean(fix.maps(:,:,9:end),3),[1 1 8]);
    end
    cmat(:,:,sub_c) = CancelDiagonals(  fix.corr,NaN);
    [~,pval(:,:,sub_c)] = corr(fix.vectorize_maps);
end
medmat = reshape(ifisherz(median(reshape(fisherz( cmat),[16 16 26]),3)),[16 16]);
imagesc(medmat)
box off
axis square
set(gca,'XTick',[4 8 12 16],'XTickLabel',{'CS+' 'CS-' 'CS+' 'CS-'},'YTick',[4 8 12 16],'YTickLabel',{'CS+' 'CS-' 'CS+' 'CS-'},'FontSize',12)
%this shows that similar faces generate similar fixation patterns i.e.
%fixation patterns are also circularly organized. Furthermore there is also
%circular similarity between baseline and test phases, which suggests that
%the behavior during the test phase is similar to the baseline. When taking
%the median across subjects, the most dissimilar pair of faces become the
%CS+ and CS?. This contradicts the pure circularity predictions and
%therefore suggests that these two faces are even more dissimilar than what
%is predicted by pure dissimilarity based on circularity. So let's look at
%what is different between CS+ and CS? faces.
%
%% same for single fixations
cmat  =[];
pval  = [];
sub_c = 0;
v = [];
correct  = 1;
tfix     = 5;

for ns = setdiff(Project.subjects_1500,[20 22 7]);
    ns
    sub_c = sub_c + 1;cond_c = 0;
    for phase = [2 4]
        for nfix = 1:tfix
            for ncond = -135:45:180
                cond_c = cond_c + 1;
                v{cond_c} = {'subject' ns 'phase' phase 'deltacsp' [ncond] 'fix' nfix};
            end            
        end
    end
    fix.getmaps(v{:})
            if correct == 1
                %correct phase specific cocktail blank
                fix.maps(:,:,1:40) = fix.maps(:,:,1:40) - repmat(mean(fix.maps(:,:,1:40),3),[1 1 40]);
                fix.maps(:,:,41:end) = fix.maps(:,:,41:end) - repmat(mean(fix.maps(:,:,41:end),3),[1 1 40]);
            end
     cmat(:,:,sub_c) = CancelDiagonals(  fix.corr,NaN);
     [~,pval(:,:,sub_c)] = corr(fix.vectorize_maps);
end
medmat = fisherz_inverse(nanmedian(fisherz( cmat),3));
imagesc(medmat)
%% differences of fixation maps
cmat     = [];
v        = [];
conds    = [-135:45:180;circshift(-135:45:180,[1 4])]
% conds    = [1:8;circshift(1:8,[1 4])];
M        = [];
cond_c   = 0;
for phase = 4
    for ncond = conds
        cond_c    = cond_c + 1;
        v{1} = {'phase' phase 'deltacsp' ncond(1) 'subject' subjects};
        v{2} = {'phase' phase 'deltacsp' ncond(2) 'subject' subjects};
        fix.getmaps(v{:});
        M    = cat(3,M,fix.maps(:,:,1)-fix.maps(:,:,2));
    end
end
fix.maps = M;
fix.plot
%this plot shows differences between fixation maps of opposite faces.
%  -135   -90   -45     0    45    90   135   180
%    45    90   135   180  -135   -90   -45     0
%in this plot only the first 4 faces carry information as same opposites
%are related to each other with a sign reversal. What I see here is that no
%matter the opposite pair I look at the right eye location is always "red" (except -135 vs. 45).
%Actually if fixation maps were to follow the facial features (which are
%changing circularly) they would be expected to change in a circular manner
%as well. This can be shown to be the case when using fileid instead of
%deltacsp. For example if the eyes are moving in space slightly due to
%different morphing levels, the difference maps taken from opposite pairs
%should basically display a rotating effect. Red and blue blobs should
%rotate across different opposite pairs (this is indeed what is happening).
%So the fact that this is not happening in the deltacsp aligned differences
%can result from the fact that different people received different faces,
%and therefore the constant influece of facial features is not anymore
%aligned with each other across faces. This is good, and the main reason
%why we switched faces for different participants. Still the mystery is
%that the right eye effect is also there in the -90 vs. 90 comparison. This
%is bad. However this could be understood if one assumes that the effect of
%conditioning is not peaking at the deltacsp 0 but somewhere closer to 45
%degrees. Then the orthogonal direction would be actually -135 vs. 45. And
%this is where we don't observe the eye effect anymore. In the next cell I
%compute some parameteric tests to show that fixation counts are
%significantly different.
%%
for kernel_fwhm = 40;
    fix.kernel_fwhm = kernel_fwhm;
    fix.unitize     = 1;
    %
    tfix = 1;
    clear D;
    for nfix = 1:tfix
        sub_c = 0;
        v = [];
        for ns = subjects
            sub_c  = sub_c + 1
            v{1}   = {'subject' ns 'phase' 4 'deltacsp' [-45]  'fix' 1:10};
            v{2}   = {'subject' ns 'phase' 4 'deltacsp' [135]  'fix' 1:10 };
            v{3}   = {'subject' ns 'phase' 4 'deltacsp' [45]   'fix' 1:10};
            v{4}   = {'subject' ns 'phase' 4 'deltacsp' [-135] 'fix' 1:10 };
            fix.getmaps(v{:});
            D(:,:,sub_c) = (fix.maps(:,:,1)-fix.maps(:,:,2)) - (fix.maps(:,:,3)-fix.maps(:,:,4));;
        end
    end
    D        = reshape(D,[size(D,1)*size(D,2) size(D,3) ]);
    [h p]    = ttest(D(:,:,1)');    
    fix.maps = reshape(-log10(p),[500 500 1]);
    fix.plot;   
end
%the effect on the right eye is positive and on the fronthead is negative.
%The effect is not there for Baseline. Also it disappears for the
%comparison 45 vs. -135. However the effect gets weaker if
%differences between 45 vs. 135 is compared to 45 vs. -135.

%% collect covariance matrices for each fixation and each subject. in the next cell we will test 
% different exploration strategies on these matrices.

for kernel_fwhm = 29
    M = [];
    C               = [];
    Cr              = [];
    fix.kernel_fwhm = kernel_fwhm;
    c_sub           = 0;
    for ns = subjects(:)'
        c_sub   = c_sub +1
        counter = 0;
        v       = [];
        for phase = [2 4]
            for nfix = 1:4
                for ncond = [-135:45:180];
                    counter      = counter + 1;
                    %                     trials = unique(fix.trialid((fix.deltacsp == ncond).*(fix.phase == phase).*(fix.subject == ns)==1));
                    v{counter}   = {'phase' phase 'deltacsp' ncond 'subject' ns 'fix' nfix(:)'};%, 'trialid', trials(1:end)};
                end
            end
        end
        fix.getmaps(v{:});
        fix.maps(:,:,1:32)      = fix.maps(:,:,1:32)   - repmat(mean(fix.maps(:,:,1:32),3),[ 1 1 32]);
        fix.maps(:,:,33:end)    = fix.maps(:,:,33:end) - repmat(mean(fix.maps(:,:,33:end),3),[ 1 1 32]);        
        M(:,:,:,c_sub)          = fix.maps;
        C(:,:,c_sub)            = fix.cov;
        Cr(:,:,c_sub)           = fix.corr;
    end   
end
%% correlation matrix decomposition. We will now decompose the covariance matrices into 4 different components: 
%_constant component: this models the correlation (or covariance) that is
%common to all faces. Basically this is the one that is mostly affected
%when doing the cocktail blank corrections. 
%_gaussian component centered on the CS+. This models the generalization of
%explorations strategies towards neighboring faces that are not accounted
%by the circular component.
%_diagonal component: This is just the variance of fixation maps. When
%fixations are all over the place there will be low variance, where as when
%all the fixations are located at the same location this will generate high
%probabilitiy values and thus high variance. For the fit to work nicely it needs to be included.
%_circular component: This is our interesting finding of circular
%similarity of fixation maps. I am still not sure what it really means. It
%would be desirable to find a name for this covariance hypothesis (besides
%circularity). It could be highly interesting to see whether this is
%related to discrimination or feargen performance (we have these weights on
%a subject by subject basis).
[X Y]         = meshgrid(linspace(0,2*pi-2*pi/8,8));
xdata         = [X(:) Y(:)]-2.3562;
block_extract = @(mat,y,x,z) mat((1:8)+(8*(y-1)),(1:8)+(8*(x-1)),z);
% the code above creates a function handle to extract blocks from the
% covariance matrices
clear amp_*
for db = 1:8%diagonal blocks
    Y                  = block_extract(C,db,db,1:size(C,3));
    [bla, param]       = fitcorrmat(Y,'gradient');
    amp_circ(db,:)     = param(1,:);
    amp_gau(db,:)      = param(2,:);
    amp_const(db,:)    = param(3,:);
    amp_diag(db,:)     = param(4,:);
end
%%
%results saved in Dropbox/CovarianceAnalysis.
figure(2);
subplot(4,1,1);bar(nanmean(amp_const,2));title('constant component')
subplot(4,1,2);bar(nanmean(amp_circ,2));title('circular component')
subplot(4,1,3);bar(nanmean(amp_gau,2));title('gaussian component');
subplot(4,1,4);bar(nanmean(amp_diag,2));title('diagonal component');

















