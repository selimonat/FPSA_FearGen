%% Discrimination
p               = Project;
allsubs = sort ([Project.subjects_1500 Project.subjects_600]);
mask            = p.getMask('ET_discr');
subjects        = intersect(find(mask),allsubs);
mask            = p.getMask('PMF');
subjects        = intersect(find(sum(mask,2)==4),subjects);
g               = Group(subjects);
g.getSI(3);
[mat,tag] = g.parameterMat;

%prepare fixmaps
phases =1;
fix = Fixmat(subjects,phases);
% collect single trials
scale            = .1;%use .1 to have the typical 2500 pixel maps
final_size       = prod(fix.rect(end-1:end).*scale);
trialnumber      = [400 120 124 240 400];
ttrial           = length(subjects)*sum(trialnumber(phases));
datamatrix       = NaN(final_size,ttrial);
labels.sub       = NaN(1,ttrial);
labels.g1500      = NaN(1,ttrial);
labels.phase     = NaN(1,ttrial);
labels.trial     = NaN(1,ttrial);
labels.cond      = NaN(1,ttrial);
labels.pos       = NaN(1,ttrial);
v = [];
c=0;
for sub = g.ids'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr};
            fprintf('subject %d phase %d trial %d\n',sub,ph,tr);
            fix.getmaps(v);
            if ~any(isnan(fix.maps(:)))                
                c                   = c+1;
                %scale it if necessary
                if scale ~= 1
                    fix.maps        = imresize(fix.maps,scale,'method','bilinear');
                end                                
                datamatrix(:,c)     = fix.vectorize_maps;
                labels.sub(c)       = sub;
                labels.g1500(c)     = ismember(sub,Project.subjects_1500);
                labels.phase(c)     = ph;
                labels.trial(c)     = tr;
                labels.cond(c)      = unique(fix.deltacsp(fix.selection));
                if ismember(ph,[1 5])
                    labels.pos(c)   = mod(tr,2);%1 is 1st, 0 is 2nd
                else
                    labels.pos(c)   = NaN;
                end
            end
        end
    end
end

%cut the nans
todelete = isnan(sum(datamatrix));
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end

 path = fullfile(homedir,'Google Drive','EthnoMaster','data','midlevel','singletrialfixmaps','phase1');
 save(fullfile(savepath,'datamatrix.mat'),'datamatrix')
 save(fullfile(savepath,'labels.mat'),'labels')
%% 
%% feargen 
clear all
p               = Project;
allsubs = sort ([Project.subjects_1500 Project.subjects_600]);
mask            = p.getMask('ET_feargen');
subjects        = intersect(find(mask),allsubs);
mask            = p.getMask('RATE');
subjects        = intersect(find(mask),subjects);
g               = Group(subjects);
g.getSI(3);
[mat,tag] = g.parameterMat;

%prepare fixmaps
phases = 4;
fix = Fixmat(subjects,phases);
% collect single trials
scale            = .1;%use .1 to have the typical 2500 pixel maps
final_size       = prod(fix.rect(end-1:end).*scale);
trialnumber      = [400 120 124 240 400];
ttrial           = length(subjects)*sum(trialnumber(phases));
datamatrix       = NaN(final_size,ttrial);
labels.sub       = NaN(1,ttrial);
labels.g1500      = NaN(1,ttrial);
labels.phase     = NaN(1,ttrial);
labels.trial     = NaN(1,ttrial);
labels.cond      = NaN(1,ttrial);
labels.pos       = NaN(1,ttrial);
v = [];
c=0;
for sub = g.ids'
    for ph = phases
        for tr = 1:max(fix.trialid(fix.phase == ph))
            v = {'subject' sub, 'phase' ph 'trialid' tr};
            fprintf('subject %d phase %d trial %d\n',sub,ph,tr);
            fix.getmaps(v);
            if ~any(isnan(fix.maps(:)))                
                c                   = c+1;
                %scale it if necessary
                if scale ~= 1
                    fix.maps        = imresize(fix.maps,scale,'method','bilinear');
                end                                
                datamatrix(:,c)     = fix.vectorize_maps;
                labels.sub(c)       = sub;
                labels.g1500(c)     = ismember(sub,Project.subjects_1500);
                labels.phase(c)     = ph;
                labels.trial(c)     = tr;
                labels.cond(c)      = unique(fix.deltacsp(fix.selection));
                if ismember(ph,[1 5])
                    labels.pos(c)   = mod(tr,2);%1 is 1st, 0 is 2nd
                else
                    labels.pos(c)   = NaN;
                end
            end
        end
    end
end

%cut the nans
todelete = isnan(sum(datamatrix));
fprintf('Will delete %g trials...\n',sum(todelete));
datamatrix(:,todelete)=[];
labels.sub(:,todelete)=[];
labels.phase(:,todelete)=[];
labels.trial(:,todelete)=[];
labels.cond(:,todelete)=[];
labels.pos(:,todelete)=[];

c = 0;
for l = unique(labels.sub)
    c = c + 1;
    labels.easy_sub(labels.sub == l) = c;
end
 path = fullfile(homedir,'Google Drive','EthnoMaster','data','midlevel','singletrialfixmaps','phase4');
 save(fullfile(savepath,'datamatrix.mat'),'datamatrix')
 save(fullfile(savepath,'labels.mat'),'labels')
