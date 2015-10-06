%% discrimination only
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = intersect(find(mask),Project.subjects_600);
g               = Group(subjects);

phase=5;
fix = Fixmat(subjects,phase);

v = [];
c = 0;
for i = 1:length(subjects)
    correct = sort(g.subject{i}.load_paradigm(phase).psi.log.globaltrial(g.subject{i}.load_paradigm(phase).psi.log.response==1));
    correct = sort([correct*2;correct*2-1])';
    wrong   = sort(g.subject{i}.load_paradigm(phase).psi.log.globaltrial(g.subject{i}.load_paradigm(phase).psi.log.response==0));
    wrong   = sort([wrong*2;wrong*2-1])';
    for trialid = {correct wrong}
        c    = c+1;
        v{c} = {'trialid' trialid{1}};
    end
end

fix.getmaps(v{:});

dummy = fix.maps;
fix.maps = mean(dummy(:,:,1:2:end),3)- mean(dummy(:,:,2:2:end),3);
fix.plot

%% try statistical test on that
dummy2 =dummy(:,:,1:2:end)-dummy(:,:,2:2:end);
dummy2=reshape(dummy2,500*500,size(dummy2,3));

h = NaN(500*500,1);
pval = NaN(500*500,1);
for i = 1:length(dummy2)
    i
   distr = dummy2(i,:);
   if sum(dummy2)~= 0
   [h(i) pval(i)] = ttest(distr);
   end
end
    

%% correct vs wrong in interval around threshold
clear all
p               = Project;
mask            = p.getMask('ET_Discr');
subjects        = intersect(find(mask),Project.subjects_1500);
mask            = p.getMask('PMF');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);

phase=1;
fix = Fixmat(subjects,phase);

v = [];
c = 0;
for i = 1:length(subjects)
    sub = g.subject{i}.id
    alpha_csp = g.subject{i}.pmf.params1(3,1);
    alpha_csn = g.subject{i}.pmf.params1(4,1);
    %direct neighbors N
    N_csp = [11.25*floor(alpha_csp/11.25) 11.25*ceil(alpha_csp/11.25)];
    N_csn = [11.25*floor(alpha_csn/11.25) 11.25*ceil(alpha_csn/11.25)];
    %shorten these guys a bit
    x        = abs(g.subject{i}.load_paradigm(phase).psi.log.x);
    trial    = g.subject{i}.load_paradigm(phase).psi.log.globaltrial;
    response = g.subject{i}.load_paradigm(phase).psi.log.response;
    %read out relevant trials (within interval around threshold)
    csp_interval_trials = trial(1,ismember(x(1,:),N_csp));
    csn_interval_trials = trial(2,ismember(x(2,:),N_csn));
    csp_correct = trial(1,(response(1,:)==1));
    csp_wrong   = trial(1,(response(1,:)==0));
    csn_correct = trial(2,(response(2,:)==1));
    csn_wrong   = trial(2,(response(2,:)==0));
    %now only trials within interval
    csp_correct = intersect(csp_interval_trials,csp_correct);
    csp_wrong   = intersect(csp_interval_trials,csp_wrong);
    csn_correct = intersect(csn_interval_trials,csn_correct);
    csn_wrong   = intersect(csn_interval_trials,csn_wrong);
    
    
    %sanity check
    if  ~isempty(intersect(csp_correct,csp_wrong))
        fprintf('something wrong in the trial computation')
        pause
    end
    
    correct_csp = sort([csp_correct*2 csp_correct*2-1])';
    wrong_csp   = sort([csp_wrong*2 csp_wrong*2-1])';
    correct_csn = sort([csn_correct*2 csn_correct*2-1])';
    wrong_csn   = sort([csn_wrong*2 csn_wrong*2-1])';
    
    correct =[csp_correct csn_correct];
    wrong   =[csp_wrong csn_wrong];
    
    for trialid = {correct wrong}
        c    = c+1;
        v{c} = {'subject' sub 'trialid' trialid{1}};
    end
end

fix.getmaps(v{:});
% 
dummy = fix.maps;
fix.maps = mean(dummy(:,:,1:2:end),3)- mean(dummy(:,:,2:2:end),3);
fix.plot

%%
