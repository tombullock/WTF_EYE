function SNR_Analysis_Compute(subjects)

%{
SNR_Analysis_Compute
Author:Tom Bullock
Date: 05.30.20
%}

sourceDirAlpha = '/home/waldrop/Desktop/WTF_EYE/EEG_Bandpassed';
sourceDirNonAlpha = '/home/waldrop/Desktop/WTF_EYE/EEG_Bandpassed_Non_Alpha';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% % select subjects
% subjects = [1:7,9:14,16:20,22:27,31];

% select electrodes
P1_elects = [25,26,27,29,30,62:64];

% loop through subs
for iSub=1:length(subjects)
    
   sjNum=subjects(iSub)
   
   for iCond=1:4
       
       load([sourceDirAlpha '/' sprintf('sj%02d_cond%02d_Alpha.mat',sjNum,iCond)])
       
       alphaMat = mean((abs(band.eeg).^2),3);
       
       clear band
       
       load([sourceDirNonAlpha '/'  sprintf('sj%02d_cond%02d_non_Alpha.mat',sjNum,iCond)])
       
       nonAlphaMat = std((abs(band.eeg).^2),0,3);
       
       clear band
       
       snrMat(iSub,iCond,:,:) = alphaMat./nonAlphaMat;
       
       
   end
end

save([destDir '/' 'SNR_Data.mat'],'snrMat')