% insert this into EEG script to get condition structures into the correct
% counterbalancing order to match with the continuous EEG file

load('sj02_allBeh.mat')

% one long struct for EEG (order according to cbOrder)
longStruct=[masterStruct(cbOrder(1)).allTrialData,...
            masterStruct(cbOrder(2)).allTrialData,...
            masterStruct(cbOrder(3)).allTrialData,...
            masterStruct(cbOrder(4)).allTrialData];

% load EEG
EEG=pop_loadset('sj01_wtfEye.set')

% resample data
EEG = pop_resample(EEG,256);

% save data in matlab format
save('sj01_wtfEye.mat','EEG')

% epoch data
load('sj01_wtfEye.mat')
EEG=pop_epoch(EEG,{102},[0 3]);

