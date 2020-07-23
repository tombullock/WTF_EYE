clear
close all

load('/home/waldrop/Desktop/WTF_EYE/Beh_Data_Processed/sj01_allBeh.mat')
load('/home/waldrop/Desktop/WTF_EYE/EYE/Processed_EYE/s1_112.mat')

for i=1:length(masterStruct(1,1).allTrialData)
    
   checkEvents(i,1:2) = [masterStruct(1,1).allTrialData(i).locTrigger, eyeEpoch.targCode(i)];
   checkEvents(i,3) = checkEvents(i,1)-checkEvents(i,2);
    
end