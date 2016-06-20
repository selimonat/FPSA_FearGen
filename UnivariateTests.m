%% 
%this script provides results that could be used to complement the
%super-complex hyperplan analysis by doing basically down to earth
%comparisons of amount of fixation counts. The question it tries to answer
%is whether and where there are more fixations on the face in CS+ vs. CS?
%conditions. The simple comparison for 0 vs. 180 gives nice results which
%are in line with the hyperplane analysis. However one needs to have a
%control and here also the best control is to check whether 90 vs. -90
%gives similar results. The idea would be that if these differences are
%induces purely by differences of the images, it should also be present in
%any combination of opposite faces.
fix             = Fixmat(setdiff(Project.subjects_1500,7),[2 4]);
p = Project;
subjects        = intersect(find(p.getMask('ET_feargen')),p.subjects_1500)';%22 has to be kicked for svm
%% similarity of fixation maps
cmat  =[];
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
end
imagesc(fisherz_inverse(   mean(fisherz( cmat),3)))
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

%% correlation matrix decomposition
[X Y]         = meshgrid(linspace(0,2*pi-2*pi/8,8));
xdata         = [X(:) Y(:)]-2.3562;

block_extract = @(mat,y,x,z) mat((1:8)+(8*(y-1)),(1:8)+(8*(x-1)),z);
Y             = block_extract(CC,6,6,1);

[model]=FitCorrMat(Y)



















