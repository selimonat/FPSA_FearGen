
%get fixations
clear all
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g = Group(subjects);
g.getSI(3)

v = [];
c = 0;
sc=0;
for sub = subjects'
    fix             = Fixmat(sub,[4]);
    sc=sc+1;
    cc=0;
    for cond = -135:45:180
        cc=cc+1;
        v = {'deltacsp' cond 'subject' sub};
        fix.getmaps(v);
        %compute differences
        i = unique(fix.file(fix.deltacsp==0));
        j = unique(fix.file(fix.deltacsp==cond)); %csp+delta
        v1diff = abs(im(:,:,i) - im(:,:,j));
        %crop the inner 500px
        v1diff = v1diff((size(v1diff,1)/2-fix.window):(size(v1diff,1)/2+fix.window)-1,(size(v1diff,2)/2-fix.window):(size(v1diff,2)/2+fix.window)-1);
        %and now reshape
        v1diff = reshape(v1diff,size(v1diff,1)*size(v1diff,2),size(v1diff,3));
        [r(sc,cc)] = corr2(fix.vectorize_maps,v1diff);
    end
    
end


%% prepare different V1 models
%get fixations
clear all
r=[];rr=[];
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);
fix4            = Fixmat(subjects,4);
fix2            = Fixmat(subjects,2);
g.getSI(3);%get the sharpening index
%
f0              = 'C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\V1\modV1\';%f0 needs the path of the images
%f is the output of CalibrateFace_V1Model and will lead to V1.mats
path2v1         = CalibrateFace_V1Model(f0);%generate V1 activity maps
dummy           = dir([path2v1 '*.mat']);
v1files         = [repmat([fileparts(path2v1) filesep],length(dummy),1) vertcat(dummy(:).name)];
tfiles          = size(v1files,1);
im              = [];
for i = 1:tfiles
    dummy       = load(v1files(i,:)); 
    im(:,:,i)   = dummy.v1;
end
%crop the inner 500px
im      = im((size(im,1)/2-fix4.window):(size(im,1)/2+fix4.window)-1,(size(im,2)/2-fix4.window):(size(im,2)/2+fix4.window)-1,:);
%%
smoothcount     = 0;
rr              = [];rrc = [];r =[];
for smooth = 15;%this loops through different smoothings (fwhm_kernel for fixmat)
    fix4.kernel_fwhm = smooth;%comes from for-loop
    fix2.kernel_fwhm = smooth;
    smoothcount      = smoothcount+1
    % get the corresponding fixmats.
    subc             = 0;
    for sub = subjects(:)'
        subc            = subc+1;
        cc              = 0;
        for cond = -135:45:180
            cc          = cc+1;
            v           = {'deltacsp' cond 'subject' sub 'fix' [1]};
            fix4.getmaps(v);            
            %compute differences
            i           = unique(fix4.file((fix4.deltacsp == 0)    & (fix4.subject == sub)));
            j           = unique(fix4.file((fix4.deltacsp == cond) & (fix4.subject == sub))); %csp+delta
            if length(i) ~= 1 && length(j) ~= 1
                keyboard
            end
            v1diff                   = abs(zscore(im(:,:,i) - im(:,:,j)));
            rr(subc,cc,smoothcount)  = v1diff(:)'*fix4.vectorize_maps;%testphase correlation
            %% get control map from phase 2 and compute control correlation
            fix2.getmaps(v);                        
            rrc(subc,cc,smoothcount) = v1diff(:)'*fix2.vectorize_maps;  %baseline correlation          
        end
    end
    r(1,smoothcount) = corr2(nanmean(rr(:,:,smoothcount),2),g.SI);%for one V1/smoothing combination
end
bar(nanmean(rr-rrc))

%% where is the most information in V1?
path2v1 = 'C:\Users\Lea\Documents\Experiments\FearCloud_Eyelab\V1\';
dummy           = dir([path2v1 '*.mat']);
v1files         = [repmat([fileparts(path2v1) filesep],length(dummy),1) vertcat(dummy(:).name)];
v1files         = v1files([1,5,9,13,17,21,25,29],:);

tfiles          = size(v1files,1);
im              = [];
for i = 1:tfiles
    dummy       = load(v1files(i,:)); 
    im(:,:,i)   = dummy.v1;
end

%correlation matrix
imagesc(corr(reshape(im,[400*400,8])));
axis square
box off
colormap(flipud(gray(256)));
cbh=colorbar('eastoutside');
set(cbh,'YTick',.96:.01:1)
set(gca,'FontSize',12)
ylabel('Face Number','FontSize',12)
xlabel('Face Number','FontSize',12)
 

 modulohelper = [8 2; 1 3;2 4; 3 5; 4 6; 5 7; 6 8; 7 1];
for a = 1:size(im,3)
    diff(:,:,a) = abs(im(:,:,a)-mean(im(:,:,modulohelper(a,:)),3));
end
