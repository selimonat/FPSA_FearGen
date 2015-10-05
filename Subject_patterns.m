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

errorbar(1:5,mean(spdiags(corrmat(1:28,:),kth),1),std(spdiags(corrmat(1:28,:),kth),1)vm,,m.v.v,m. /sqrt(length(subjects)),'bo-',...
    'LineWidth',2.5)
ylabel('mean diag')
xlabel('nth diagonal')
set(gca,'XTick',[1 2 3 4 5])
ylim([0.5 1.2]);
