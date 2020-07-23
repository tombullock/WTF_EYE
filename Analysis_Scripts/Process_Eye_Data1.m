%{
Process_Eye_Data
Author: Tom Bullock
Date: 02.13.20

Notes: epoch around Trigger 102 and create matrix of all trials 

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Raw_EYE';
destDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Processed_EYE_UPDATED';

cd(sourceDir)
d=dir('s23*.mat');

for f=1:length(d)
    
    disp(num2str(f))
    
    % load data file
    load([d(f).folder '/' d(f).name])
    
    triggerInfo = eyeMat.event.Messages.info;
    triggerTimes = eyeMat.event.Messages.time;
    sampleTimes = double(eyeMat.RawEdf.FSAMPLE.time);
    xPos = double(eyeMat.RawEdf.FSAMPLE.gx(2,:));
    yPos = double(eyeMat.RawEdf.FSAMPLE.gy(2,:));
    pa = double(eyeMat.RawEdf.FSAMPLE.pa(2,:));
    
    cnt=0;
    for i=1:length(triggerInfo)
        
        if strcmp(triggerInfo(i),'TRIGGER 102')
            cnt=cnt+1;
            
            % get timestamp index (102) for each trial
            [~,startTimeIdx] = min(abs(triggerTimes(i)-sampleTimes));
            
            % get timestamps, x/y positions and pupil diameter
            eyeEpoch.allTimes(cnt,:) = sampleTimes(startTimeIdx:startTimeIdx+1730); % 1750 = 3500 ms post 102 (should be enough to capture full trial?)
            eyeEpoch.allXpos(cnt,:) = xPos(startTimeIdx:startTimeIdx+1730);
            eyeEpoch.allYpos(cnt,:) = yPos(startTimeIdx:startTimeIdx+1730);
            eyeEpoch.allPA(cnt,:) = pa(startTimeIdx:startTimeIdx+1730);
            
            % get code and timestamp for target (if target present)
            if ismember(str2double(triggerInfo{i+1}(end-2:end)),201:208)
                eyeEpoch.targCode(cnt,:) = str2double(triggerInfo{i+1}(end-2:end));
                eyeEpoch.targTime(cnt,:) = triggerTimes(i+1);
            else
                eyeEpoch.targCode(cnt,:) = 0;
                eyeEpoch.targTime(cnt,:) = 0;
            end
        end
        
    end
    
    save([destDir '/' d(f).name],'eyeEpoch')
    
    clear cnt eyeEpoch eyeMat pa sampleTimes triggerInfo triggerTimes xPos yPos
        
end