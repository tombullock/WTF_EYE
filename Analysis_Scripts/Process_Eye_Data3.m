%{
Process_Eye_Data_For_IEM_Stuff
Author:Tom Bullock
Date: 04.29.20

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Synchronized_EYE_UPDATED';
destDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Sync_Final';

% select subject number (just one)
sjNum=14;

% load data
load([sourceDir '/' sprintf('sj%02d_eye_beh_sync.mat',sjNum)])

% loop through sessions and conditions
for iCond=1:4
    for iSession=1:2
        
        % isolate eye data from struct
        eyeData = masterStruct(iSession,iCond).eyeData;
        
        % isolate beh data from struct
        trialData = masterStruct(iSession,iCond).allTrialData;
        
        % add subject exceptions
        if sjNum==10 && iCond==2 && iSession==2
           trialData(end) = []; 
        elseif sjNum==23 && iCond==1 && iSession==1
            trialData(end) = [];
        elseif sjNum==14 && iCond==2 && iSession==1
            trialData(end) = [];
        end
        
        % get vector of broken fixation trials
        brokenFixTrials = [trialData.brokenTrial];
        
        % remove broken fixation trials from trial data
        trialData(brokenFixTrials==1) = [];
        
        % remove broken fixation trials from eye data
        eyeData.allTimes = eyeData.allTimes(brokenFixTrials==0,:);
        eyeData.allXpos = eyeData.allXpos(brokenFixTrials==0,:);
        eyeData.allYpos = eyeData.allYpos(brokenFixTrials==0,:);
        eyeData.allPA = eyeData.allPA(brokenFixTrials==0,:);
        eyeData.targCode = eyeData.targCode(brokenFixTrials==0,:);
        eyeData.targTime = eyeData.targTime(brokenFixTrials==0,:);
        
        % loop through trials
        for iTrial=1:size(eyeData.allTimes,1)
            
            % get target time
            targTime = eyeData.targTime(iTrial);
            
            % get all times
            allTimes = eyeData.allTimes(iTrial,:);
            
            % get index of target time relative to all times (find closet value)
            [v,idx] = min(abs(targTime-allTimes));
            
            % generate time seg index for target onset > end retention period (1000 ms)
            tIdx = (idx-125):(idx+988); % final 11 samples lost?
            
            % create new eye movement and pupil size mats based on trial time segment
            try
                eyeProcessed(iSession,iCond).ex(iTrial,:) = eyeData.allXpos(iTrial,tIdx);
                eyeProcessed(iSession,iCond).ey(iTrial,:) = eyeData.allYpos(iTrial,tIdx);
                eyeProcessed(iSession,iCond).pa(iTrial,:) = eyeData.allPA(iTrial,tIdx);
            catch
                eyeProcessed(iSession,iCond).ex(iTrial,:) = nan(1,length(tIdx));
                eyeProcessed(iSession,iCond).ey(iTrial,:) = nan(1,length(tIdx));
                eyeProcessed(iSession,iCond).pa(iTrial,:) = nan(1,length(tIdx));
            end
            
            clear targTime allTimes v idx tIdx 
            
        end
        
        eyeProcessed(iSession,iCond).trialData = trialData;
        eyeProcessed(iSession,iCond).times = (-125:988)*2;
        
        clear eyeData trialData brokenFixTrials 
        
    end
    
end

save([destDir '/' sprintf('sj%02d_eye_beh_sync_pro.mat',sjNum)],'eyeProcessed')








% 
% % compute euclidian distance from baseline eye position and plot
% for i=1:size(eyeEpoch.allXpos,1)
%     
%     bx=mean(eyeEpoch.allXpos(i,1:250));
%     by=mean(eyeEpoch.allYpos(i,1:250));
%     
%     x=eyeEpoch.allXpos(i,:);
%     y=eyeEpoch.allYpos(i,:);
%     
%     
%     for j=1:length(x)
%         e(j) = sqrt(((bx-x(j))^2) + ((by-y(j))^2));  
%     end
%     
%     plot(e)
%     
%     pause(1)
% end