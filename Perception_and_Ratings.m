%% general stats
p               = Project;
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),Project.subjects_600);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);
g.getSI(3)

[mat labels] = g.parameterMat;

SI_norm =mean(g.sigma_cond-g.sigma_test)./mean(g.sigma_cond);

%plot something for that
means  = [mean(mean(mat(:,1:4),2)),mean(mat(:,end-2)),mean(mat(:,end-1))];
errors = [std(mean(mat(:,1:4),2)),std(mat(:,end-2)),std(mat(:,end-1))]/sqrt(length(subjects));
[hbarr errbar]=barwitherr(errors,means);
ylim([0 65])
set(errbar,'LineWidth',2)
axis off

%% model Ratings as Gaussian

%%1500ms
p               = Project;
g1 = Group([6 7 9 10 11 12 13 14 15 16 17 18 19 20 21 22 24 25 26 28 30 33 34 35]);%people with acceptable fits for ratings and PMF
n1 = length(g1.ids);
CSP_subject = mean(g1.pmf.params1([1,3],1,:),1);
CSP_mean = mean(CSP_subject);
g1.getSI(3);
FWHM=[];
for i=1:length(g1.tunings{4}.singlesubject)
    FWHM(i,:)= [g1.tunings{3}.singlesubject{i}.Est(2) g1.tunings{4}.singlesubject{i}.Est(2)];
end
FWHM_mean = mean(FWHM,1);
AMPL=[];
for i=1:length(g1.tunings{4}.singlesubject)
    AMPL(i,:)= [g1.tunings{3}.singlesubject{i}.Est(1) g1.tunings{4}.singlesubject{i}.Est(1)];
end
AMPL_mean = mean(AMPL,1);

alpha_incr=[];
for i=1:length(g1.tunings{4}.singlesubject)
alpha_incr = [alpha_incr; ((g1.pmf.params1(1,1,i)-g1.pmf.params1(3,1,i))-(g1.pmf.params1(2,1,i)-g1.pmf.params1(4,1,i)))];
end


%%
%600ms
p               = Project;
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),Project.subjects_600);
g2              = Group(subjects);
n2 = length(g2.ids);

g2.getSI(3);
FWHM=[];
for i=1:length(g2.tunings{4}.singlesubject)
    FWHM(i,:)= [g2.tunings{3}.singlesubject{i}.Est(2) g2.tunings{4}.singlesubject{i}.Est(2)];
end
FWHM_mean = mean(FWHM,1);
AMPL=[];
for i=1:length(g2.tunings{4}.singlesubject)
    AMPL(i,:)= [g2.tunings{3}.singlesubject{i}.Est(1) g2.tunings{4}.singlesubject{i}.Est(1)];
end
AMPL_mean = mean(AMPL,1);

alpha_incr=[];
for i=1:length(g2.tunings{4}.singlesubject)
alpha_incr = [alpha_incr; ((g2.pmf.params1(1,1,i)-g2.pmf.params1(3,1,i))-(g2.pmf.params1(2,1,i)-g2.pmf.params1(4,1,i)))];
end


%% look at perception and ratings  - %do good discriminators generalize less?
clear all;
p         = Project;
rate_mask = find(p.getMask('RATE'));
pmf_mask  = find(sum(p.getMask('PMF'),2) == 4);
subs      = intersect(intersect(rate_mask,pmf_mask),Project.subjects_1500);
g         = Group(subs);
g.getSI(3);
[mat labels] = g.parameterMat;


med  = median(mean(mat(:,[1 2 3 4]),2));

good = g.ids(find(mean(mat(:,[1 2 3 4]),2) < med));
bad  = g.ids(find(mean(mat(:,[1 2 3 4]),2) >= med));
alpha_good = mean(mean(mat(ismember(g.ids,good),[1 2 3 4]),2));
alpha_bad  = mean(mean(mat(ismember(g.ids,bad),[1 2 3 4]),2));

sigma_cond_good = mat(ismember(g.ids,good),12);
sigma_cond_bad  = mat(ismember(g.ids,bad),12);
sigma_test_good = mat(ismember(g.ids,good),12);
sigma_test_bad  = mat(ismember(g.ids,bad),12);
SI_good         = mat(ismember(g.ids,good),end);
SI_bad          = mat(ismember(g.ids,bad),end);

[H,P,CI,STATS] = ttest(sigma_cond_good,sigma_cond_bad)
[H,P,CI,STATS] = ttest(sigma_test_good,sigma_test_bad)
[H,P,CI,STATS] = ttest(SI_good,SI_bad)
%nope.




