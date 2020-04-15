%{
EEG_Merge_Broken_Files
Author: Tom Bullock
Date: 11.21.19
%}

clear
close all

% add paths
cd('D:\\WTF_EYE\\Analysis_Scripts');
cd('D:\\WTF_EYE\\Dependancies\\eeglab14_1_2b');
eeglab
close all
cd('D:\\WTF_EYE\\Analysis_Scripts');

% set dir
sourceDir = 'D:\\WTF_EYE\\EEG_Raw';

sjNum=14;
session=1;

% set file parts
filePart1 = sprintf('sj%02d_se%02d_wtfEye.bdf',sjNum,session);
filePart2 = sprintf('sj%02d_se%02d_wtfEye_1.bdf',sjNum,session);

% load parts
EEG1=pop_biosig([sourceDir '\\' filePart1]);
EEG2=pop_biosig([sourceDir '\\' filePart2]);

% % special case three files merged
filePart3 = sprintf('sj%02d_se%02d_wtfEye_2.bdf',sjNum,session);
EEG3=pop_biosig([sourceDir '\\' filePart3]);
% EEG = pop_mergeset(EEG,EEG3);

% merge
EEG = pop_mergeset(EEG1,EEG2);
EEG = pop_mergeset(EEG,EEG3);

% save
save([sourceDir '\\' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum,session)],'EEG','-v7.3')