%{
IEM_Cross_TT_Compile_Slopes
Author: Tom Bullock
Date: 04.15.20

Loops through all the slope files and compiles them all into MASTER
matrices for plotting/stats analysis in R

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Cross';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% select subs
subjects = [1:7,9:14,16:20,22:27,31];

for iSub=1:length(subjects)
    sjNum = subjects(iSub);
    load([sourceDir '/' sprintf('sj%02d_IEM_Alpha_CrossTT.mat',sjNum)])
    all_rSl(iSub,:,:,:) = em.rSl.total;
    all_pSl(iSub,:,:,:) = em.pSl.total;
end

% save compiled data
save([destDir '/' 'Cross_TT_Slope_MASTER.mat'],'all_rSl','all_pSl','-v7.3')
