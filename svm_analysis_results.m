%% 2 - subjects within phases
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\svm_analysis\20151121T163214\subjects_inphase_rand0\result.mat')
a = squeeze(mean(result,3));
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,27]);%scale by rowsums
mean(diag(mean(scaled,3)))
imagesc(mean(scaled,3)); axis square;
colorbar;caxis([0 1])
figure
for i = 1:5
    subplot(1,5,i);imagesc(scaled(:,:,i));
    title(['phase' num2str(i)]);colorbar;caxis([0 1]);axis off;axis square;colorbar off;
end
for i = 1:5;m(i) = mean(diag(scaled(:,:,i)));s(i) = std(diag(scaled(:,:,i)));end
[h e] = barwitherr(s,m);
set(e,'LineWidth',1.5)
%% 3 - phases (unbalanced)
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\svm_analysis\20151117T163214\phases_rand0\result.mat')
a = squeeze(mean(result,3)); %mean over bootstraps
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,5]); %scale by rowsums
mean(diag(mean(scaled,3)))
imagesc(mean(scaled,3)); axis square; set(gca,'XTick',1:5,'YTick',1:5)
colorbar;caxis([0 1])

%% 4 - phases within subjects
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\svm_analysis\20151121T163214\phases_insubject_rand0\result.mat')
a = squeeze(mean(result,3)); %mean over bootstraps
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,5]); %scale by rowsums
mean(diag(mean(scaled,3)))
imagesc(mean(scaled,3)); axis square; set(gca,'XTick',[1:5],'YTick',1:5)
colorbar;caxis([0 1])
for i = 1:5;
    subplot(1,5,i);
    imagesc(scaled(:,:,i));
    title(['phase' num2str(i)]);
    colorbar;caxis([0 0.92]);
    axis off;axis square;
    colorbar off;
end
colorbar
%classification for the 27 subjects, barplot
for n=1:27;acc(n) = mean(diag(scaled(:,:,n)));err(n) = std(diag(scaled(:,:,n)));end
[h,e]=barwitherr(err,acc)
set(e,'LineWidth',1.5)
% confusion business
conf45 = squeeze(scaled(4,5,:));
conf41 = squeeze(scaled(4,1,:));
p = Project;
mask = p.getMask('RATE');
sub_c = unique(labels.sub);
ind = ismember(sub_c,intersect(find(mask),sub_c));
conf45 = conf45(ind);
conf41 = conf41(ind);
g = Group(sub_c(ind));
g.getSI(3);
[mat tags] = g.parameterMat;
[r,p]=corr(conf45,mat(:,12))
[r,p]=corr(conf45,mat(:,13))
[r,p]=corr(conf45,mat(:,14))
[r,p]=corr(conf41,mat(:,12))
[r,p]=corr(conf41,mat(:,13))
[r,p]=corr(conf41,mat(:,14))
%% 6 - conditions
a = squeeze(mean(result(:,:,:,:,2),3));
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,8]);%scale by rowsums
imagesc(mean(scaled,3))
b = mean(scaled,3)
plot(b(4,:))
hold;plot(b(8,:),'r')

%% behavioral
w = mean(w,2);
%sigma_test
load('C:\Users\onat\Google Drive\EthnoMaster\data\midlevel\singletrialfixmaps\1500\SI_N24\phase4\alleigen.mat')
eigen = eigen(:,1:142);
Sigmatest_hp = reshape(eigen*w,[50 50]);
fix = Fixmat(6,4);
fix.getmaps({'subject' 6})
fix.maps = Sigmatest_hp;
fix.plot
%% 12 - csp csn classification
%within subs
for i = 1:2
a = squeeze(mean(result(:,:,:,:,i),3));
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,2]);%scale by rowsums
subplot(1,2,i);imagesc(mean(scaled,3));axis square;colorbar;caxis([0 1]);
subplot(1,2,i);set(gca,'XTick',1:2,'YTick',1:2,'XTickLabel',{'CSP','CSN'},'YTickLabel',{'CSP','CSN'})
end
subplot(1,2,1);title('baseline')
subplot(1,2,2);title('testphase')
t = supertitle('Classifying CSP vs CSN within subjects',1)
set(t,'FontSize',14)
[h,p]=ttest(scaled(1,1,:),0.5)
[h,p]=ttest(scaled(2,2,:),0.5)
%across subs
for i = 1:2
a = squeeze(mean(result(:,:,:,i),3));
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,2]);%scale by rowsums
subplot(1,2,i);imagesc(mean(scaled,3));axis square;colorbar;caxis([0 1]);
subplot(1,2,i);set(gca,'XTick',1:2,'YTick',1:2,'XTickLabel',{'CSP','CSN'},'YTickLabel',{'CSP','CSN'})
end
subplot(1,2,1);title('baseline')
subplot(1,2,2);title('testphase')
t = supertitle('Classifying CSP vs CSN grouping subjects',1)
set(t,'FontSize',14)
a = squeeze(result(:,:,:,1));
scaled = (a./sum(a(:)))./repmat(sum((a./sum(a(:))),2),[1,2]);%scale by rowsums
[h,p]=ttest(scaled(1,1,:),0.5)
[h,p]=ttest(scaled(2,2,:),0.5)