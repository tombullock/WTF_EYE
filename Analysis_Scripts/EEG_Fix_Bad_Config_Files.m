%{
EEG_Fix_Bad_Config_Files
Author: Tom Bullock
Date: 01.27.20
%}

filenameOiginal = 'sj02_se02_wtfEye-Deci.bdf';
filenameSave = 'sj02_se02_wtfEye.mat';

% load file
EEG = pop_biosig(filenameOiginal);

% crop excess channels
EEG = pop_select(EEG,'channel',1:72);

% replace missing chan info with chan info from another file
load('D:\WTF_EYE\EEG_Chanlocs.mat')
EEG.chanlocs=chanlocs;

% save file
save(filenameSave,'EEG')










% chanlocs=tmp.EEG.chanlocs;
% save('EEG_Chanlocs.mat','chanlocs')