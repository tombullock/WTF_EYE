%{
Spatial_IEM
Original Author: Joshua Foster (heavily edited by Tom Bullock, UCSB Attention Lab)
Date:11.15.19

Load preprocessed EEG and Beh data.
Merge across sessions and split into condition
Get time-frequency data
Either 1) Run IEM or 2) just get bandpassed data ...

%}

function SNR_Analysis_Get_Non_Alpha(sn)

%
% if nargin==0
%     subjects = [1:7,9:14,16:20,22:27,31];
% end
%
% %clear
% %close all

% set dirs
rDir = '/home/waldrop/Desktop/WTF_EYE';
eegDir = [rDir '/' 'EEG_Prepro2_Avg_Baseline'];
bandpassedDir = [rDir '/' 'EEG_Bandpassed_Non_Alpha'];
%iemDir = [rDir '/' 'IEM_Results_TT_Within_Tmp' ];
% select subjects
%subjects = 4;

addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts'))

% loop through subs
%for s=1:length(subjects)
%    sn=subjects(s);

% import IEM settings
em = IEM_Settings;
% nChans = em.nChans;
% nBins = em.nBins;
% nIter = em.nIter;
% nPerms = em.nPerms;
% nBlocks = em.nBlocks;
freqs = em.frequencies;
bands = em.bands;
% times = em.time;
% nFreqs = size(em.frequencies,1);
% nSamps = length(em.time);
% Fs = em.Fs;
% basisSet = em.basisSet;

% merge data across sessions and split by condition
for iSession=1:2
    
    load([eegDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sn,iSession)])
    
    % reject artifact trials from eeg and beh
    eegs = eeg.data;
    eegs = eegs(:,:,~eeg.arf.artIndCleaned);
    allBehStruct = allBehStruct(~eeg.arf.artIndCleaned);
    
    bothSessions(iSession).eegs=eegs;
    bothSessions(iSession).allBehStruct=allBehStruct;
    
    clear eegs allBehStruct
    
end

eegs = cat(3,bothSessions(1).eegs, bothSessions(2).eegs);
beh = cat(2,bothSessions(1).allBehStruct, bothSessions(2).allBehStruct);
clear bothSessions

% create condition vector
clear condBin
for i=1:length(beh)
    if      beh(i).memCondition==1 && beh(i).eyesCondition==2; condBin(i)=1;behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==2; condBin(i)=2;behBin(i)=beh(i);
    elseif  beh(i).memCondition==1 && beh(i).eyesCondition==6; condBin(i)=3;behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==6; condBin(i)=4;behBin(i)=beh(i);
    end
end

% split data by condition
for iCond=1:4
    allConds(iCond).eegs = eegs(:,:,find(condBin==iCond));
    allConds(iCond).beh = beh(find(condBin==iCond));
end

% find location bin with min count per condition (to balance IEM)
for iCond=1:4
    
    clear beh posBin minCnt
    beh=allConds(iCond).beh;
    
    for i=1:length(beh)
        posBin(i)=beh(i).stimLoc;
    end
    
    binCnt = [];
    for bin = 1:8
        binCnt(bin) = sum(posBin == bin);
    end
    
    allConds(iCond).minCnt = min(binCnt);
    allPosBin(iCond).posBin=posBin;
    
    clear posBin
end

% loop through conditions
for iCond=1:4
    
    % set posBin and minCnt for this cond
    posBin = allPosBin(iCond).posBin;
    
    for iMin=1:4
        minCntMat(iMin)=allConds(iMin).minCnt;
    end
    
    minCnt = min(minCntMat)-1; % minus 1 to give all bins variability across iterations
    
    %%clear eegs
    
    
    % loop through frequency bands (typically alpha, theta)
    for f=1%:size(freqs,1)
        
        clear data bandEEG dataFilt1
        
        % grab eegs
        eegs = allConds(iCond).eegs;
        
        % apply butterworth filter (3rd order, bandpass)
        [z1,p1] = butter(3, [freqs(f,1), freqs(f,2)]./(eeg.sampRate/2),'stop');
        data = double(eegs);
        bandEEG = NaN(size(data,1),size(data,2),size(data,3));
        for x = 1:size(data,1)
            for y = 1:size(data,3)
                dataFilt1 = filtfilt(z1,p1,data(x,:,y));
                bandEEG(x,:,y) = dataFilt1;
            end
        end
        
        % apply hilbert tranform to bandpassed data
        eegs = [];
        for i=1:size(bandEEG,1) % chan loop
            i
            for j=1:size(bandEEG,3) % trial loop
                eegs(i,:,j) = hilbert(squeeze(bandEEG(i,:,j)));
            end
        end
        
%         % take std over trials
%         eegs = std(eegs,3);
%         
%         % average over trials to reduce matrix size
%         eegs = squeeze(mean(eegs,3))
        
        %average over all trials to reduce the size of the matrix
        %%rawHilbert = squeeze(std(fdata_evoked,1));
        
        % create bandpassed data structure and save for sub\cond\freq
        band.eeg = eegs;
        band.freqs = freqs(f,:);
        band.chanlocs = eeg.chanLabels;
        band.srate = eeg.sampRate;
        band.beh = allConds(iCond).beh;
        band.times = eeg.times;
        save([bandpassedDir '/' sprintf('sj%02d_cond%02d_non_%s',sn,iCond,bands{f})],'band','-v7.3');
        
    end
end


