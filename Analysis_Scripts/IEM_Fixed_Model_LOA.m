%{
IEM_Fixed_Model
Author: Tom Bullock, UCSB Attention Lab (using code borrowed from Joshua Foster)
Date created: 11.15.19
Date updated: 10.23.20

Purpose:

Load preprocessed EEG and Beh data.
Merge across sessions and split into condition
Get time-freq data
Run IEM, training on all conditions and testing each condition separately

Inputs: sn (subject number)

Outputs: 

%}

function IEM_Fixed_Model_LOA(sn)

% set dirs
rDir = '/home/waldrop/Desktop/WTF_EYE';
eegDir = [rDir '/' 'EEG_Prepro2_Avg_Baseline'];
iemDir = [rDir '/' 'IEM_Results_TT_Within_Fixed' ];

% add paths
addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts'))

% import IEM settings from structure
em = IEM_Fixed_Model_LOA_Settings;
nBins = em.nBins;
nIter = em.nIter;
freqs = em.frequencies;
bands = em.bands;
times = em.time;
nFreqs = size(em.frequencies,1);
nSamps = length(em.time);
Fs = em.Fs;
basisSet = em.basisSet;

% merge data across sessions and split by condition
for iSession=1:2
    
    % load data
    load([eegDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sn,iSession)])
    
    % reject artifact trials from eeg and beh
    eegs = eeg.data;
    eegs = eegs(:,:,~eeg.arf.artIndCleaned);
    allBehStruct = allBehStruct(~eeg.arf.artIndCleaned);
    
    % combine sessions
    bothSessions(iSession).eegs=eegs;
    bothSessions(iSession).allBehStruct=allBehStruct;
    
    clear eegs allBehStruct
    
end

% create merged session matrices
eegs = cat(3,bothSessions(1).eegs, bothSessions(2).eegs); % [elects x times x trials]
beh = cat(2,bothSessions(1).allBehStruct, bothSessions(2).allBehStruct); % [trials]
clear bothSessions

% create vector of condition labels
clear condBin
for i=1:length(beh)
    if      beh(i).memCondition==1 && beh(i).eyesCondition==2; condBin(i)=1;%behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==2; condBin(i)=2;%behBin(i)=beh(i);
    elseif  beh(i).memCondition==1 && beh(i).eyesCondition==6; condBin(i)=3;%behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==6; condBin(i)=4;%behBin(i)=beh(i);
    end
end

% split eeg and beh data by condition
for iCond=1:4
    allConds(iCond).eegs = eegs(:,:,find(condBin==iCond));
    allConds(iCond).beh = beh(find(condBin==iCond));
end

clear condBin

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
    
    allConds(iCond).minCnt = min(binCnt); % add minCnt to allConds structure
    allPosBin(iCond).posBin=posBin; % create structure with just position bins
    
    clear posBin
end

% create a continuous vector of position bins across all conditions from allPosBin
posBin = [...
    allPosBin(1).posBin,...
    allPosBin(2).posBin,...
    allPosBin(3).posBin,...
    allPosBin(4).posBin...
    ];

% find min bin cnt across conditions (-1 ensures variable
minCnt = min([allConds.minCnt]) - 1;

% create index to keep track 
condIdx = [...
    repmat(1,1,length(allPosBin(1).posBin)),...
    repmat(2,1,length(allPosBin(2).posBin)),...
    repmat(3,1,length(allPosBin(3).posBin)),...
    repmat(4,1,length(allPosBin(4).posBin))...
    ];

% loop through frequency bands
for f=1:size(freqs,1) 
        
    % apply butterworth filter (3rd order, bandpass)
    disp('Applying Butterworth Filter')
    [z1,p1] = butter(3, [freqs(f,1), freqs(f,2)]./(eeg.sampRate/2),'bandpass');
    data = double(eegs);
    bandEEG = NaN(size(data,1),size(data,2),size(data,3));
    for x = 1:size(data,1)
        for y = 1:size(data,3)
            dataFilt1 = filtfilt(z1,p1,data(x,:,y));
            bandEEG(x,:,y) = dataFilt1;
        end
    end
    
    clear data dataFilt1
    
    % apply hilbert tranform to bandpassed data
    disp('Applying Hilbert Tranformation')
    eegs = [];
    for i=1:size(bandEEG,1) % chan loop
        %i
        for j=1:size(bandEEG,3) % trial loop
            eegs(i,:,j) = hilbert(squeeze(bandEEG(i,:,j)));
        end
    end
    
    clear bandEEG
    
    % remove bad channels for IEM analyses
    badChanIdx = [];
    cnt=0;
    for i=1:length(eeg.badChannels)
        thisBadChan=eeg.badChannels{i};
        badChanIdx(i) = find(strcmp(thisBadChan,eeg.chanLabels)==1);
    end
    eegs(badChanIdx,:,:)=[];
    
    clear badChanIdx thisBadChan
    
    % change data shape to (trials x elects x times) for purposes of IEM script
    eegs = permute(eegs,[3,1,2]);    
    
    % convert data to total power
    fdata_total = abs(eegs).^2;

    % get nElectrodes
    nElectrodes = size(eegs,2); 
    em.nvElectrodes = nElectrodes;
    
    % determine timepoints to run analysis on
    tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time);
    
    % loop through each iteration
    for iter = 1:nIter
        
        disp(['nIteration :' num2str(iter) ' of ' num2str(nIter)])
        
        % shuffle trials
        shuffInd = randperm(length(posBin))'; % create shuffled vector of trial numbers
        shuffBin = posBin(shuffInd); % create shuffled vector of loc bins
        shuffCond = condIdx(shuffInd); % create shuffled vector of condition indices
        
        % loop through and create
        for iCond=1:4
            for bin=1:nBins
                idx = find(shuffBin == bin & shuffCond == iCond); % get index for trials belonging to the current bin && condition
                idx = idx(1:minCnt); % drop excess trials from idx for balancing purposes
                trialIdxBalancedMat(iCond,bin,:) = idx; % create balanced matrix of trials
                binIdxBalancedMat(iCond,bin,:) = repmat(bin,[1,minCnt]); % create balanced matrix of bins
                condIdxBalancedMat(iCond,bin,:) = repmat(iCond,[1,minCnt]); % create balanced matrix of conditions 
                clear idx
            end
        end
        
        % use colon operator to transform balanced matrices into vectors
        trialIdxBal = trialIdxBalancedMat(:);
        binIdxBal = binIdxBalancedMat(:);
        condIdxBal = condIdxBalancedMat(:);
        
        % create 'C' (matrix of basis functions for each trial)
        for iGenC = 1:length(binIdxBal)
            C(iGenC,:) = basisSet(binIdxBal(iGenC),:);
        end
        
        % create balanced set of EEG total power trials
        eegDataBal = fdata_total(trialIdxBal,:,:);

        clear trialIdxBalancedMat binIdxBalancedMat condIdxBalancedMat
        
        
        % loop through samples
        for t=1:nSamps
            
            disp(['Sample ' num2str(t) ' of ' num2str(nSamps)])
            
            % define timepoint t 
            toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
           
            % get data for specified timepoint
            dt = mean(eegDataBal(:,:,toi),3);
            
            % create encoding model (em) loop index
            emIdx = 1:size(dt,1);
            
            % run model
            cnt1=0;cnt2=0;cnt3=0;cnt4=0;
            for iEM=emIdx
                
                % create training and testing indexes (coulda just used
                % cvpartition here)
                trnIdx = emIdx;
                trnIdx(iEM)=[];
                tstIdx = iEM;
                
                % get C1
                C1 = C(trnIdx,:);
                
                % get B1
                B1 = dt(trnIdx,:);
                
                % get a weights matrix (w)
                tWeights = C1\B1;
                
                % get test data
                B2 = dt(tstIdx,:);
                
                % set up centering
                centerind = find(C(iEM,:)==1);
                centeredC1(iEM,:) = circshift(C(iEM,:)',5-centerind)'; % reality check - the C1 should be uniformly centered

                % test (each condition)
                if condIdxBal(iEM)==1
                    cnt1=cnt1+1;
                    C2_cond1(cnt1,:) = (tWeights'\B2')';
                    C2_cond1_centered(cnt1,:) = circshift(C2_cond1(cnt1,:)',5-centerind);
                elseif condIdxBal(iEM)==2
                    cnt2=cnt2+1;
                    C2_cond2(cnt2,:) = (tWeights'\B2')';
                    C2_cond2_centered(cnt2,:) = circshift(C2_cond1(cnt2,:)',5-centerind);
                elseif condIdxBal(iEM)==3
                    cnt3=cnt3+1;
                    C2_cond3(cnt3,:) = (tWeights'\B2')';
                    C2_cond3_centered(cnt3,:) = circshift(C2_cond1(cnt3,:)',5-centerind);
                elseif condIdxBal(iEM)==4
                    cnt4=cnt4+1;
                    C2_cond4(cnt4,:) = (tWeights'\B2')';
                    C2_cond4_centered(cnt4,:) = circshift(C2_cond1(cnt4,:)',5-centerind);
                end
                
            end
            

            %average across trials in each condition
            allTF(1,f,iter,t,:) = mean(C2_cond1_centered,1);
            allTF(2,f,iter,t,:) = mean(C2_cond2_centered,1);
            allTF(3,f,iter,t,:) = mean(C2_cond3_centered,1);
            allTF(4,f,iter,t,:) = mean(C2_cond4_centered,1);

            
            clear C2_cond1 C2_cond1_centered C2_cond2 C2_cond2_centered C2_cond3 C2_cond3_centered C2_cond4 C2_cond4_centered
           
        end
    end
end

% save tfs to structure
em.allTF = allTF;

% save data
save([iemDir '/' sprintf('sj%02d_fixed_IEM.mat',sn)],'em');
          
return