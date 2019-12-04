%{
EEG_Merge_Broken_Files
Author: Tom Bullock
Date: 11.21.19
%}

clear
close all

% add paths
cd('C:\\Users\\BOSS-EEG\\Desktop\\WTF_EYE\\Analysis_Scripts');
cd('C:\\Users\\BOSS-EEG\\Desktop\\WTF_EYE\\Dependancies\\eeglab14_1_2b');
eeglab
close all
cd('C:\\Users\\BOSS-EEG\\Desktop\\WTF_EYE\\Analysis_Scripts');

% set dir
sourceDir = 'C:\\Users\\BOSS-EEG\\Desktop\\WTF_EYE\\EEG_Raw';

sjNum=3;
session=2;

% set file parts
filePart1 = sprintf('sj%02d_se%02d_wtfEye.bdf',sjNum,session);
filePart2 = sprintf('sj%02d_se%02d_wtfEye_1.bdf',sjNum,session);

% load parts
EEG1=pop_biosig([sourceDir '\\' filePart1]);
EEG2=pop_biosig([sourceDir '\\' filePart2]);

% merge
EEG = pop_mergeset(EEG1,EEG2);

% save
save([sourceDir '\\' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum,session)],'EEG','-v7.3')