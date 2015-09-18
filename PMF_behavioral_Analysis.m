%analyze PMF data on the behavioral level

%%1500ms

g = Group(intersect(Project.subjects_1500,find(sum(PMFmask,2)==4))); %use only these people that have data for all 4 PMFs
n = length(g.ids);

%600ms
g = Group(intersect(Project.subjects_600,find(sum(PMFmask,2)==4))); 
n = length(g.ids);


%%compute mean/std
alpha.mean=mean(g.pmf.params1(:,1,:),3);
alpha.std=std(squeeze(g.pmf.params1(:,1,:)),0,2);
alpha.sem=alpha.std/sqrt(n);
beta.mean=mean(g.pmf.params1(:,2,:),3);
beta.std=std(squeeze(g.pmf.params1(:,2,:)),0,2);
beta.sem=beta.std/sqrt(n);

gamma.mean=mean(g.pmf.params1(:,3,:),3);
gamma.std=std(squeeze(g.pmf.params1(:,3,:)),0,2);
gamma.sem=alpha.std/sqrt(n);
lambda.mean=mean(g.pmf.params1(:,4,:),3);
lambda.std=std(squeeze(g.pmf.params1(:,4,:)),0,2);
lambda.sem=beta.std/sqrt(n);


%% Scatterplot and TTEST
fig1=figure;
subplot(2,2,1)%CS+ alpha
plot(squeeze(g.pmf.params1(1,1,:)),squeeze(g.pmf.params1(3,1,:)),'ro','MarkerFaceColor','r');
ylabel(sprintf('after'))
axis square
DrawIdentityLine(gca);
title('CS+ alpha')
box off
hold on;
plot(alpha.mean(1),alpha.mean(3),'k+','MarkerSize',10)
subplot(2,2,2)%CS+ beta
plot(log10(squeeze(g.pmf.params1(1,2,:))),log10(squeeze(g.pmf.params1(3,2,:))),'ro','MarkerFaceColor','r');
title('CS+ log10(beta)')
axis square
DrawIdentityLine(gca);
box off
hold on;
plot(log10(beta.mean(1)),log10(beta.mean(3)),'k+','MarkerSize',10)

subplot(2,2,3)%CS- alpha
plot(squeeze(g.pmf.params1(2,1,:)),squeeze(g.pmf.params1(4,1,:)),'bo','MarkerFaceColor','b');
xlabel('before')
ylabel(sprintf('after'))
axis square
DrawIdentityLine(gca);
title('CS- alpha')
box off
hold on;
plot(alpha.mean(2),alpha.mean(4),'k+','MarkerSize',10)

subplot(2,2,4)%CS- beta
plot(log10(squeeze(g.pmf.params1(2,2,:))),log10(squeeze(g.pmf.params1(4,2,:))),'bo','MarkerFaceColor','b');
xlabel(sprintf('before'))
axis square
DrawIdentityLine(gca);
title('CS- log10(beta)')
box off
hold on;
plot(log10(beta.mean(2)),log10(beta.mean(4)),'k+','MarkerSize',10)


t=supertitle('Alpha and Beta before and after aversive learning (600ms)',1);
set(t,'FontSize',14)

[h,p,ci]=ttest(squeeze(g.pmf.params1(1,1,:)),squeeze(g.pmf.params1(3,1,:)))%CS+ alpha before/after p<0.001***
[h,p,ci]=ttest(squeeze(g.pmf.params1(2,1,:)),squeeze(g.pmf.params1(4,1,:)))%CS- alpha before/after p =0.11
[h,p,ci]=ttest(squeeze(g.pmf.params1(1,2,:)),squeeze(g.pmf.params1(3,2,:)))%CS+ beta before/after  p =0.36
[h,p,ci]=ttest(squeeze(g.pmf.params1(2,2,:)),squeeze(g.pmf.params1(4,2,:)))%CS- beta before/after  p =0.53

%%Median Split


