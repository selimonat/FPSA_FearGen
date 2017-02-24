%% run
SVM_simulations_BDNF(2);
SVM_simulations_BDNF(4,'1');
SVM_simulations_BDNF(4,'2');
SVM_simulations_BDNF(4,'3');
SVM_simulations_BDNF(4);


%% plot

clear all
clf
neig = 8;

load(['C:\Users\Lea\Documents\Experiments\project_bdnf\data\midlevel\svm_analysis\findingparams\SVM_simulations_phase2_NEV8_FWHM30_alltrials.mat'])
M = mean(mean(result,2),3);
SE = std(mean(result,2),[],3)./sqrt(size(result,3));
figure(100);
subplot(1,5,1);
b = bar(M);SetFearGenBarColors(b);
hold on;
errorbar(M,SE,'k.','LineWidth',2)
ylim([.3 .7])
set(gca,'YTick',.3:.1:.7,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'});
box off
axis square
ylabel('Classified as CS+')

for tt = 1:4
    if tt< 4
        load(['C:\Users\Lea\Documents\Experiments\project_bdnf\data\midlevel\svm_analysis\findingparams\SVM_simulations_phase4_NEV8_FWHM30_part' sprintf('%d.mat',tt)])
    else
        load('C:\Users\Lea\Documents\Experiments\project_bdnf\data\midlevel\svm_analysis\findingparams\SVM_simulations_phase4_NEV8_FWHM30_alltrials.mat')
    end
    RR(:,tt) = R;
    AVEmap(:,:,tt) = AVEHP(:,:);
    M(:,tt) = mean(mean(result,2),3);
    SE(:,tt) = std(mean(result,2),[],3)./sqrt(size(result,3));
end

for tt = 1:4
subplot(1,5,tt+1)
b = bar(RR(:,tt));SetFearGenBarColors(b);
hold on;
errorbar(M(:,tt),SE(:,tt),'k.','LineWidth',2)
ylim([.3 .7])
set(gca,'YTick',.3:.1:.7,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'});
box off
axis square
end

labels = {'baseline' 'test_1' 'test_2' 'test_3' 'test_{all}'};
for n = 1:5
    subplot(1,5,n)
    ylim([.25 .75])
    set(gca,'YTick',[.3 .5 .7])
    xlim([0 9])
    title(labels{n})
end
%% add chance level
hold on;
load(['C:\Users\Lea\Documents\Experiments\project_bdnf\data\midlevel\svm_analysis\findingparams\SVM_simulations_phase2_NEV8_FWHM30_r1_alltrials.mat'])
M = mean(mean(result,2),3);
SE = std(mean(result,2),[],3)./sqrt(size(result,3));
subplot(1,5,1);
l = line(xlim,[M(1) M(8)]);set(l,'Color','k','LineStyle',':','LineWidth',2)

for tt = 1:4
    if tt< 4
        load(['C:\Users\Lea\Documents\Experiments\project_bdnf\data\midlevel\svm_analysis\findingparams\SVM_simulations_phase4_NEV8_FWHM30_r1__part' sprintf('%d.mat',tt)])
    else
        load('C:\Users\Lea\Documents\Experiments\project_bdnf\data\midlevel\svm_analysis\findingparams\SVM_simulations_phase4_NEV8_FWHM30_r1_alltrials.mat')
    end
    RR(:,tt) = R;
    AVEmap(:,:,tt) = AVEHP(:,:);
    M(:,tt) = mean(mean(result,2),3);
    SE(:,tt) = std(mean(result,2),[],3)./sqrt(size(result,3));
end

for tt = 1:4
subplot(1,5,tt+1)
l = line(xlim,[M(1,tt) M(8,tt)]);set(l,'Color','k','LineStyle',':','LineWidth',2)
end

%% try to staple them
close all
figure(101);
for tt = [3 2 1]
b = bar((1:8)+tt*.2,RR(:,tt));SetFearGenBarColors(b);set(get(b,'Children'),'FaceAlpha',tt/3-.1)
hold on
ylim([.3 .7])
set(gca,'YTick',.3:.1:.7,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'});
end


figure(1);
fix = Fixmat(6,2);
fix.maps = AVEmap;
fix.plot('linear',[1 3])
