%{
EEG_Compile_Bandpassed_Data
Author: Tom Bullock
Date: 04.22.20

Grabs bandpassed data and compiles into a big grand matrix for alpha topo
plotting purposes
%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/EEG_Bandpassed';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% set subject numbers
subjects = [1:7,9:14,16:20,22:27,31];

for iSub=1:length(subjects)
    sjNum=subjects(iSub)
    for iCond=1:4 
        
        % load bandpassed data
        load([sourceDir '/' sprintf('sj%02d_cond%02d_Alpha.mat',sjNum,iCond) ])  
        
        % baseline correct
        thisBand = band.eeg;
        
        thisBand = thisBand - mean(thisBand(:,1:128,:),2);
        
        
        %allBand = allBand - mean(allBand(:,:,:,:,1:128),5);
        band.eeg = thisBand;
        
        clear thisBand

        % convert hilbert complex doubles to power (induced) before avg
        band.eeg = abs(band.eeg).^2;
        
        % get stimulus location list
        clear stimLocs
        for b=1:length(band.beh) 
            stimLocs(b) = band.beh(b).stimLoc;
        end
        
        % parse by stimulus location (average over location trials)
        for iLoc=1:8
            allBand(iSub,iCond,iLoc,:,:) = mean(band.eeg(:,:,find(stimLocs==iLoc)),3);       
        end
        
        % log info to save as struct (all band files contain the same info)
        dataInfo.chanlocs = band.chanlocs;
        dataInfo.freqs = band.freqs;
        dataInfo.times = band.times;
        dataInfo.srate = band.srate;
        
        clear band
        
    end
    
    
end



% save data
save([destDir '/' 'Bandpassed_Data_ALpha.mat'],'allBand','dataInfo','-v7.3')
