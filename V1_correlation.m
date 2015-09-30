

%get fixations
clear all
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_600);
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
p               = Project;
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g = Group(subjects);
g.getSI(3)

%f0 needs the path of the images
%f0 = ~/Google Drive/EthnoMaster/data/midlevel/V1/modV1/;
%f is the output of CalibrateFace_V1Model and will lead to V1.mats

smoothcount = 0;

for smooth = 15:5:50%this loops through different smoothings (fwhm_kernel for fixmat)
    smoothcount = smoothcount+1;
    scalecount = 0;
    for nScale = 1:10 %loops through different frequencies in V1-Model
        scalecount = scalecount+1;
        
        f = CalibrateFace_V1Model(f0,nScale);
        
        %load V1 images (computed before in this for-loop)
        path2v1 = f;
        dummy = dir([path2v1 '*.mat']);
        v1files    = [repmat([fileparts(path2v1) filesep],length(dummy),1) vertcat(dummy(:).name)];
        tfiles = size(v1files,1);
        im=[];
        c=0;
        for i = [1:tfiles]
            c=c+1;
            dummy = load(v1files(i,:)); dummy=dummy.v1;
            im(:,:,c) = dummy;
        end
        
        % get the corresponding fixmats.
        subc = 0;
        for sub = subjects'
            subc=subc+1;
            cc=0;
            
            fix             = Fixmat(sub,4); % needs to go subject by subject
            fix.maptype     = 'conv';
            fix.kernel_fwhm = smooth;%comes from for-loop
            
            for cond = -135:45:180
                cc=cc+1;
                v = {'deltacsp' cond 'subject' sub};
                fix.getmaps(v);
                %compute differences
                i        = unique(fix.file(fix.deltacsp==0));
                j        = unique(fix.file(fix.deltacsp==cond)); %csp+delta
                v1diff   = abs(im(:,:,i) - im(:,:,j));
                %crop the inner 500px
                v1diff   = v1diff((size(v1diff,1)/2-fix.window):(size(v1diff,1)/2+fix.window)-1,(size(v1diff,2)/2-fix.window):(size(v1diff,2)/2+fix.window)-1);
                %and now reshape
                v1diff   = reshape(v1diff,size(v1diff,1)*size(v1diff,2),size(v1diff,3));
                rr(subc,cc) = corr2(fix.vectorize_maps,v1diff);
            end
        end
        rr(smoothcount,scalecount) = mean(nanmean(r,2),1);%for one V1/smoothing combination
    end
end

imagesc(rr)

