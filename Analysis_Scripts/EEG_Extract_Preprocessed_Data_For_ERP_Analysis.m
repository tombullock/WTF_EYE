function EEG_Extract_Preprocessed_Data_For_ERP_Analysis(sn)

%{
ERP_Extract_Preprocessed_Data_For_ERP_Analysis
Author: Tom Bullock
Date: 05.28.20

Grab the preprocessed EEG data ready for ERP analysis

%}

% clear
% close all

% set dirs
rDir = '/home/waldrop/Desktop/WTF_EYE';
eegDir = [rDir '/' 'EEG_Prepro2_Avg_Baseline'];
destDir = [rDir '/' 'EEG_Processed'];

% add analysis scripts to path
addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts'))

% % select subs
% subjects = [1:7,9:14,16:20,22:27,31];

% for iSub=1:length(subjects)
%     
%     sn=subjects(iSub)
    
% merge data across sessions and split by condition
for iSession=1:2
    
    load([eegDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sn,iSession)])
    
    % reject artifact trials from eeg and beh
    eegs = eeg.data;
    eegs = eegs(:,:,~eeg.arf.artIndCleaned);
    allBehStruct = allBehStruct(~eeg.arf.artIndCleaned);
    
    bothSessions(iSession).eegs=eegs;
    bothSessions(iSession).allBehStruct=allBehStruct;
    
    clear eegs allBehStruct
    
end

eegs = cat(3,bothSessions(1).eegs, bothSessions(2).eegs);
beh = cat(2,bothSessions(1).allBehStruct, bothSessions(2).allBehStruct);
clear bothSessions

% create condition vector
clear condBin
for i=1:length(beh)
    if      beh(i).memCondition==1 && beh(i).eyesCondition==2; condBin(i)=1;behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==2; condBin(i)=2;behBin(i)=beh(i);
    elseif  beh(i).memCondition==1 && beh(i).eyesCondition==6; condBin(i)=3;behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==6; condBin(i)=4;behBin(i)=beh(i);
    end
end

% split data by condition
for iCond=1:4
    allConds(iCond).eegs = eegs(:,:,find(condBin==iCond));
    allConds(iCond).beh = beh(find(condBin==iCond));
end

% get channel locs
chanlocs = eeg.chanlocsOriginal;
times = eeg.times;

% save data
save([destDir '/' sprintf('sj%02d_eeg.mat',sn)],'allConds','chanlocs','times')

clear allConds beh eegs bothSessions allBehStruct
    
%end
