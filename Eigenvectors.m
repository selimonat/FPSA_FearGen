clear all
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('PMF');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);


phase=1;
fix = Fixmat(subjects,phase);
v=[];
c=0;
for sub = subjects'
    c=c+1;
  v{c} = {'subject' sub};
end

%%
fix.maptype = 'conv';
% fix.binsize = 15;
fix.getmaps(v{:});
fix.maps = imresize(fix.maps,[50 50]);   %if conv, make it less data
%show subjects' personal fixation maps
fix.plot
%% eigenvector computation
subject_maps = fix.vectorize_maps;
fprintf('starting eigenvector computation\n')
[e dv] = eig(cov(fix.vectorize_maps'));
fprintf('done\n')
close all
plot(diag(dv((400-size(fix.maps,1)):end,(400-size(fix.maps,1)):end)),'bo-')
%mostly the last 5 are most telling, so only take 5 best
eigen = fliplr(e);
%plot them on faces (best eigenvector first)
num = 10;
fix.maps = reshape(eigen(:,1:num),[size(fix.maps,1),size(fix.maps,1),num]);
fix.plot
%% compute loadings of subject on eigenvectors
ss=0;
subj_load=[];
for s = 1:length(subjects)
    ss=ss+1;
subj_load(ss,:) = subject_maps(:,s)'*eigen(:,1:num);
end
close all
%exemplary plotting some loadings
for i=1:5;subplot(1,5,i);bar(subj_load(i,:));xlim([0 num+1]);title(sprintf('subject %d',5+i));end

%take initial discrimination threshold, correlate it with these loadings
initial_alpha = mean([g.pmf.csp_before_alpha,g.pmf.csn_before_alpha],2);
[r,p]=corr(subj_load,initial_alpha);


%%

%% across all phases

clear all
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),subjects);
mask            = p.getMask('PMF');
subjects        = intersect(find(mask),subjects);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);

phase=[1 2 3 4 5];
fix = Fixmat(subjects,phase);
v=[];
c=0;
for sub = subjects'
    c=c+1;
    v{c} = {'subject' sub};
end
%%
fix.maptype = 'conv';
% fix.binsize = 15;
fix.getmaps(v{:});
fix.maps = imresize(fix.maps,[50 50]);   %if conv, make it less data
%show subjects' personal fixation maps
fix.plot
%% eigenvector computation
subject_maps = fix.vectorize_maps;
fprintf('starting eigenvector computation\n')
[e dv] = eig(cov(fix.vectorize_maps'));
fprintf('done\n')
close all
plot(diag(dv((480-size(fix.maps,1)):end,(480-size(fix.maps,1)):end)),'bo-')
dv = sort(diag(dv),'descend');plot(cumsum(dv)./sum(dv),'o-')
%mostly the last 5 are most telling, so only take 5 best
eigen = fliplr(e);
%plot them on faces (best eigenvector first)
num = 10;
fix.maps = reshape(eigen(:,1:num),[size(fix.maps,1),size(fix.maps,1),num]);
fix.plot
%% compute loadings of subject on eigenvectors
ss=0;
subj_load=[];
for s = 1:length(subjects)
    ss=ss+1;
subj_load(ss,:) = subject_maps(:,s)'*eigen(:,1:num);
end
close all
%exemplary plotting some loadings
for i=1:5;subplot(1,5,i);bar(subj_load(i,:));xlim([0 num+1]);title(sprintf('subject %d',5+i));end

%take initial discrimination threshold, correlate it with these loadings
initial_alpha = mean([g.pmf.csp_before_alpha,g.pmf.csn_before_alpha],2);
[r,p]=corr(subj_load,initial_alpha);
later_alpha   = mean([g.pmf.csp_after_alpha,g.pmf.csn_after_alpha],2);
[r,p]=corr(subj_load,later_alpha);

[mat labels]   = g.parameterMat;
mean_improvmnt = mean(mat(:,9:10),2);
[r,p]=corr(subj_load,mean_improvmnt);
incr_improvmnt = mat(:,11);
[r,p]=corr(subj_load,incr_improvmnt);

%SI
g.getSI(3)
[r,p]=corr(subj_load,g.SI);
[r,p]=corr(subj_load,g.sigma_cond);
[r,p]=corr(subj_load,g.sigma_test);

%% which eigenvector dominant in which phase?
%get phase maps
v=[];
c=0;
for ph = phase
    c=c+1;
    v{c} = {'phase' ph};
end
fix.maptype = 'conv';
fix.getmaps(v{:});
fix.maps = imresize(fix.maps,[50 50]); 
phase_maps = fix.vectorize_maps;

ps=0;
phase_load=[];
for s = 1:length(phase)
    ps=ps+1;
phase_load(ps,:) = subject_maps(:,s)'*eigen(:,1:num);
end
close all
%plot the loadings per phase
for i=1:5;subplot(1,5,i);bar(subj_load(i,:));xlim([0 num+1]);title(sprintf('ph %d',i));end

