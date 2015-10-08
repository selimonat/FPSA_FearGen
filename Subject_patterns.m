%% 
clear all
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);

phase=[1 2 3 4 5];
fix = Fixmat(subjects,phase);
v=[];
c=0;

for ph = phase
    for sub = subjects'
        c=c+1;
        v{c} = {'subject' sub 'phase' ph};
    end
end



%%
fix.maptype = 'conv';
% fix.binsize = 15;
fix.getmaps(v{:});
fix.maps = imresize(fix.maps,[50 50]);   %if conv, make it less data
%% fixation types
corrmat = corrcoef(fix.vectorize_maps).^2;
figure;imagesc(corrmat);
corrmat = corrmat(1:28,:);
kth = [0 ([1:4]*length(subjects))];

errorbar(1:5,mean(spdiags(corrmat(1:28,:),kth),1),std(spdiags(corrmat(1:28,:),kth),1)./sqrt(length(subjects)),'bo-',...
    'LineWidth',2.5)
ylabel('mean diag')
xlabel('nth diagonal')
set(gca,'XTick',[1 2 3 4 5])
ylim([0.5 1.2]);
%% correlation matrix of fixation similarities across phases
clear all
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);

phase=[1 2 3 4 5];
fix = Fixmat(subjects,phase);
%% mean for each phase
v = [];
c = 0;
for ph = phase
        c = c+1;
        v{c} = {'phase' ph};
end
fix.getmaps(v{:});
%if wanted, cocktailblank it
%fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 size(fix.maps,3)]);
corrmat = corr(fix.vectorize_maps);
imagesc(corrmat)
%% mean for each subject
v = [];
c = 0;
for sub = subjects'
        c = c+1;
        v{c} = {'subject' sub};
end
fix.getmaps(v{:});
corrmat = corr(fix.vectorize_maps);
imagesc(corrmat)
%% single subjects @ phases
v=[];
c=0;
for sub = subjects'
    for ph = phase
    c=c+1;
    v{c} = {'subject' sub 'phase' ph};
    end
end
fix.getmaps(v{:});
fix.maps = imresize(fix.maps,[50 50]);
corrmat = corr(fix.vectorize_maps);
imagesc(corrmat);colorbar;
set(gca,'XTick',5:5:140,'XTickLabel',[],'YTick',5:5:140,'YTickLabel',[])

% mean the 5x5 matrices along the diag
start = -4;
cc = 0;
for i = 1:length(subjects)
    start = start+5;
    cc = cc+1;
    subcor(cc) = mean(mean(corrmat(start:(start+4),start:(start+4))));
    start:(start+4)
end

%% mds scale different subjects
v=[];
c=0;
for sub = subjects'
    c=c+1;
    v{c} = {'subject' sub};
end
fix.getmaps(v{:});
fix.maps = imresize(fix.maps,[50 50]); % so the PC doesn't get out of memory
maps = fix.vectorize_maps;
ed = zeros(length(subjects));
for nf1 = 1:length(subjects)
    A = maps(:,nf1);
    for nf2 = 1:length(subjects)
        if nf2<nf1
            fprintf('Processing subjects %d-%d\n',nf1,nf2);            
            B = maps(:,nf2);
            ed(nf1,nf2) = norm(A(:)-B(:));
            ed(nf2,nf1) = ed(nf1,nf2);
        end
    end
end
imagesc(ed)
mds = mdscale(ed,2,'criterion','metricstress');

figure;plot(mds(:,1),mds(:,2),'o','MarkerFaceColor','b')
text(mds(:,1)+0.1e-4,mds(:,2),num2str(g.ids),'fontsize',10)

%% exemplary subjects vs group for phase 1 - 5
v=[];
c=0;
for sub = 10
    for ph = phase
    c=c+1;
    v{c} = {'subject' sub 'phase' ph};
    end
end
fix.getmaps(v{:});
close all

%% correlations within vs between subjects
clear all
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),subjects);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);

phase = [1 2 3 4 5];
fix   = Fixmat(subjects,phase);
v=[];
c=0;
for sub = subjects'
    for ph = phase
    c=c+1;
    v{c} = {'subject' sub 'phase' ph};
    end
end
fix.getmaps(v{:});
%cocktail blank
%fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 size(fix.maps,3)]);
fix.maps = imresize(fix.maps,[50 50]); 
subj_maps = fix.vectorize_maps;
imagesc(corr(subj_maps));
axis off
axis image
print -dbitmap;

corrmat = corr(subj_maps);
%extract a subject's consistence (across phase correlation)
subcorr=[];
start = -4;
cc = 0;
for i = 1:length(subjects)
    start = start+5;
    cc = cc+1;
    subcor(cc) = mean(squareform(corrmat(start:(start+4),start:(start+4))-1)+1);
    start:(start+4);
end
%% within phases, between subjects correlation
v=[];
c=0;
for ph = phase
    for sub = subjects'
        c=c+1;
        v{c} = {'subject' sub 'phase' ph};
    end
end
fix.getmaps(v{:});
%cocktail blank
%fix.maps = fix.maps - repmat(mean(fix.maps,3),[1 1 size(fix.maps,3)]);
fix.maps = imresize(fix.maps,[50 50]); 
subj_maps = fix.vectorize_maps;
corrmat = corr(subj_maps);
imagesc(corrmat);
axis off
axis image
print -dbitmap;

phasecor=[];
start = -length(subjects)+1;
cc=0;
for i = 1:length(phases)
    start = start+length(subjects);
    cc = cc+1;
    phasecor(cc) = mean(squareform(corrmat(start:start+(length(subjects)-1),start:start+(length(subjects)-1))-1)+1);
    start,start+(length(subjects)-1)
end


