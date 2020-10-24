%{
IEM_Fixed_Model
Original Author: Joshua Foster (heavily edited by Tom Bullock, UCSB Attention Lab)
Date:11.15.19

Load preprocessed EEG and Beh data.
Merge across sessions and split into condition
Get time-freq data
Run IEM, training on all conditions and testing each condition separately

%}

function IEM_Fixed_Model_Blocked(sn)

% if nargin==0
%     subjects = [1:7,9:14,16:20,22:27,31];
% end

% set dirs
rDir = '/home/waldrop/Desktop/WTF_EYE';
eegDir = [rDir '/' 'EEG_Prepro2_Avg_Baseline'];
%bandpassedDir = [rDir '/' 'EEG_Bandpassed'];
iemDir = [rDir '/' 'IEM_Results_TT_Within_Fixed' ];

% add paths
addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts'))

% import IEM settings
em = IEM_Settings;
nChans = em.nChans;
nBins = em.nBins;
nIter = em.nIter;
%nPerms = em.nPerms; % different for within and cross plots
nBlocks = 8;%4;%em.nBlocks; EDIT: CHANGED THIS TO 4!
freqs = em.frequencies;
bands = em.bands;
times = em.time;
nFreqs = size(em.frequencies,1);
nSamps = length(em.time);
Fs = em.Fs;
basisSet = em.basisSet;

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

% FIXED MODEL EDITS HERE ONWARDS

% find min bin cnt across all conditions
minCnt = min([allConds.minCnt]) - 1;

% get position bins across all conditions
posBin = [...
    allPosBin(1).posBin,...
    allPosBin(2).posBin,...
    allPosBin(3).posBin,...
    allPosBin(4).posBin...
    ];

%Jordan Added: create index to keep track of rest and exercise position trials
trial_indx = [...
    repmat(1,1,length(allPosBin(1).posBin)),...
    repmat(2,1,length(allPosBin(2).posBin)),...
    repmat(3,1,length(allPosBin(3).posBin)),...
    repmat(4,1,length(allPosBin(4).posBin))...
    ];




% loop through frequency bands (typically alpha, theta)
for f=1%:size(freqs,1) % just alpha for now
    
    clear data bandEEG dataFilt1
    
    % grab eegs
    %eegs = allConds(iCond).eegs;
    
%     %% get no. of electrodes
%     nElectrodes = size(eeg.data,1);
%     em.nElectrodes = nElectrodes;
%     disp(['nElecrodes changed to :' num2str(nElectrodes)])
    
    % apply butterworth filter (3rd order, bandpass)
    [z1,p1] = butter(3, [freqs(f,1), freqs(f,2)]./(eeg.sampRate/2),'bandpass');
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
    
    % create bandpassed data structure and save for sub\cond\freq
    band.eeg = eegs;
    band.freqs = freqs(f,:);
    band.chanlocs = eeg.chanLabels;
    band.srate = eeg.sampRate;
    band.beh = allConds(iCond).beh;
    band.times = eeg.times;
    %save([bandpassedDir '/' sprintf('sj%02d_cond%02d_%s',sn,iCond,bands{f})],'band','-v7.3');
    
    % remove bad channels for IEM analyses
    clear badChanIdx thisBadChan
    badChanIdx = [];
    cnt=0;
    for i=1:length(eeg.badChannels)
        thisBadChan=eeg.badChannels{i};
        badChanIdx(i) = find(strcmp(thisBadChan,eeg.chanLabels)==1);
    end
    eegs(badChanIdx,:,:)=[];
    
    % generate evoked and total data mats
    % permute data to(trials x elects x times for IEM script)
    eegs = permute(eegs,[3,1,2]);
    fdata_evoked = eegs;
    fdata_total = abs(eegs).^2;
    
    % get nElectrodes
    nElectrodes = size(eegs,2); % Jordan changed 9/8/2020
    em.nElectrodes = nElectrodes;
    
    % determine timepoints to run analysis on
    tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time);
    nTimes = length(tois); % index time points for analysis.
    
    % Loop through each iteration
    for iter = 1:nIter
        
        disp(['nIteration :' num2str(iter) ' of ' num2str(nIter)])
        
        % assign trials to blocks - preallocate arrays
        blocks = nan(size(posBin));
        shuffBlocks = nan(size(posBin));
        nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
        
        % shuffle trials
        %shuffInd = randperm(nTrials)';
        shuffInd = randperm(length(posBin))'; % create shuffle index
        shuffBin = posBin(shuffInd); % shuffle trial order
        
        % added for fixed model
        shuffTrial_indx = trial_indx(shuffInd);

        
        % take the 1st nPerBin x nBlocks trials for each position bin.
        for bin = 1:nBins
            
            idx = find(shuffBin == bin); % get index for trials belonging to the current bin
            idx = idx(1:nPerBin*nBlocks); % drop excess trials
            %x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
            
            % TOM ADDED FIXED
            idx_cond1 = find(shuffBin == bin & shuffTrial_indx == 1);
            idx_cond1 = idx_cond1(1:nPerBin*nBlocks);
            
            idx_cond2 = find(shuffBin == bin & shuffTrial_indx == 2);
            idx_cond2 = idx_cond2(1:nPerBin*nBlocks);
            
            idx_cond3 = find(shuffBin == bin & shuffTrial_indx == 3);
            idx_cond3 = idx_cond3(1:nPerBin*nBlocks);
            
            idx_cond4 = find(shuffBin == bin & shuffTrial_indx == 4);
            idx_cond4 = idx_cond4(1:nPerBin*nBlocks);
            
            fixed_TrainIdx = [idx_cond1;idx_cond2;idx_cond3;idx_cond4];
            
            
            % assign one condition to each of the four blocks
            if rem(length(idx_cond1),2)
                first_blockN = round(length(idx_cond1)/2);
                second_blockN = length(idx_cond1) - first_blockN;
                x_cond1 = [ones(1,first_blockN),repmat(2,1,second_blockN)];
            else
                x_cond1 = repmat(1:2,1,length(idx_cond1)/2);
            end
            
            if rem(length(idx_cond2),2)
                first_blockN = round(length(idx_cond2)/2);
                second_blockN = length(idx_cond2) - first_blockN;
                x_cond2 = [repmat(3,1,first_blockN),repmat(4,1,second_blockN)];
            else
                x_cond2 = repmat(3:4,1,length(idx_cond2)/2);
            end
            
            if rem(length(idx_cond3),2)
                first_blockN = round(length(idx_cond3)/2);
                second_blockN = length(idx_cond3) - first_blockN;
                x_cond3 = [repmat(5,1,first_blockN),repmat(6,1,second_blockN)];
            else
                x_cond3 = repmat(5:6,1,length(idx_cond3)/2);
            end
            
            if rem(length(idx_cond4),2)
                first_blockN = round(length(idx_cond4)/2);
                second_blockN = length(idx_cond4) - first_blockN;
                x_cond4 = [repmat(7,1,first_blockN),repmat(8,1,second_blockN)];
            else
                x_cond4 = repmat(7:8,1,length(idx_cond4)/2);
            end
            
            fixed_x = [x_cond1, x_cond2, x_cond3, x_cond4];
            shuffBlocks(fixed_TrainIdx) = fixed_x;
           
            
        end
        
        % unshuffle block assignment
        blocks(shuffInd) = shuffBlocks;
        
        % FIXED MODEL EDIT!
        block_conMat = [blocks, shuffTrial_indx];
        blockCon_trialNum = [];
        for iBlock = 1:nBlocks
            for iCon = 1:4
                blockCon_trialNum(iBlock, iCon) = sum(block_conMat(:,1) == iBlock & block_conMat(:,2) == iCon);
            end
        end
        
        
        
        % save block assignment
        em.blocks(:,iter) = blocks; % block assignment
        em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
        
        % FIXED MODEL ADDED
        em.trial_Indx(:,iter) = trial_indx; %save condition assignment for each trial
        
        % Average data for each position bin across blocks
        posBins = 1:nBins;
        blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
        blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
        
        blockDat_evoked_cond1 = [];%nan(nBins*nBlocks/2,nElectrodes,nSamps);
        blockDat_total_cond1 = nan(size(blockDat_evoked_cond1));
        
        blockDat_evoked_cond2 = nan(size(blockDat_evoked_cond1));
        blockDat_total_cond2 = nan(size(blockDat_evoked_cond2));
        
        blockDat_evoked_cond3 = nan(size(blockDat_evoked_cond1));
        blockDat_total_cond3= nan(size(blockDat_evoked_cond3));
        
        blockDat_evoked_cond4 = nan(size(blockDat_evoked_cond1));
        blockDat_total_cond4 = nan(size(blockDat_evoked_cond4));
        
        
        
        labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
        blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
        c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
        bCnt = 1;
        
%         % Jordan added
%         rest_labels = [];%nan(nBins*nBlocks/2,1);
%         rest_blockNum = nan(size(rest_labels));
%         rest_c =[];%nan(nBins*nBlocks/2,nChans);
%         
%         ex_labels = nan(size(rest_labels));
%         ex_blockNum = nan(size(ex_labels));
%         ex_c = nan(size(rest_c));
        
        cond1_bCnt = 1;
        cond2_bCnt = 1;
        cond3_bCnt = 1;
        cond4_bCnt = 1;
        
        for ii = 1:nBins
            for iii = 1:nBlocks
              %  blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
               blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
                
                
                % Seperate condition data
                if (iii == 1 || iii == 2) % Rest condition
                    blockDat_evoked_cond1(cond1_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond1(cond1_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                    
                    
                    % Con labels for testing
                    cond1_labels(cond1_bCnt) = ii;
                    cond1_blockNum(cond1_bCnt) = iii;
                    cond1_c(cond1_bCnt,:) = basisSet(ii,:);
                    
                    cond1_bCnt = cond1_bCnt + 1;
                elseif (iii == 3 || iii == 4) % Exercise condition
                    blockDat_evoked_cond2(cond2_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond2(cond2_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                    
                    cond2_labels(cond2_bCnt) = ii;
                    cond2_blockNum(cond2_bCnt) = iii;
                    cond2_c(cond2_bCnt,:) = basisSet(ii,:);
                    
                    cond2_bCnt = cond2_bCnt + 1;
                    
                elseif (iii == 5 || iii == 6) % Exercise condition
                    blockDat_evoked_cond3(cond3_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond3(cond3_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                    
                    cond3_labels(cond3_bCnt) = ii;
                    cond3_blockNum(cond3_bCnt) = iii;
                    cond3_c(cond3_bCnt,:) = basisSet(ii,:);
                    
                    cond3_bCnt = cond3_bCnt + 1;
                    
                elseif (iii == 7 || iii == 8) % Exercise condition
                    blockDat_evoked_cond4(cond4_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond4(cond4_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                    
                    cond4_labels(cond4_bCnt) = ii;
                    cond4_blockNum(cond4_bCnt) = iii;
                    cond4_c(cond4_bCnt,:) = basisSet(ii,:);
                    
                    cond4_bCnt = cond4_bCnt + 1;
                    
                    
                end
                
                %rest
                %blockDat_evoked_Rest(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks_rest==iii,:,tois),1))).^2;
                %blockDat_total_Rest(bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks_rest==iii,:,tois),1)));
                
                %exercise
                %blockDat_evoked_Low(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks_low==iii,:,tois),1))).^2;
                %blockDat_total_Low(bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks_low==iii,:,tois),1)));
                
                labels(bCnt) = ii; %provide labels for training and testing
                blockNum(bCnt) = iii;
                c(bCnt,:) = basisSet(ii,:);
                %cond_blockDat(bCnt,:) = trial_indx(posBin == posBins(ii) & blocks == iii,1);
                bCnt = bCnt+1;
                
                
                
            end
        end

        
        %-------------------------------------------------------------------------
        
%         % Average data for each position bin across blocks
%         nElectrodes = size(fdata_total,2);
%         posBins = 1:nBins;
%         blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
%         blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
%         labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
%         blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
%         c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
%         bCnt = 1;
%         for ii = 1:nBins
%             for iii = 1:nBlocks
%                 blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
%                 blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
%                 labels(bCnt) = ii;
%                 blockNum(bCnt) = iii;
%                 c(bCnt,:) = basisSet(ii,:);
%                 bCnt = bCnt+1;
%             end
%         end
        
        %%% IEM STUFF %%%
        
        %==========================================================
        %{
                Run high temporal resolution IEM for "within" analysis only
        %}
        %==========================================================
        
        for t = 1:nSamps
            
            % grab data for timepoint t
            toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
            %find(toi==1)
            de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
            dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
            
            % Do forward model
            
%             cond1_blockNum=1;
%             cond2_blockNum=1;
%             cond3_blockNum=1;
%             cond4_blockNum=1;
            
            
            %Jordan added
            de_cond1 = squeeze(mean(blockDat_evoked_cond1(:,:,toi),3));
            dt_cond1 = squeeze(mean(blockDat_total_cond1(:,:,toi),3));
            
            de_cond2 = squeeze(mean(blockDat_evoked_cond2(:,:,toi),3));
            dt_cond2 = squeeze(mean(blockDat_total_cond2(:,:,toi),3));
            
            de_cond3 = squeeze(mean(blockDat_evoked_cond3(:,:,toi),3));
            dt_cond3 = squeeze(mean(blockDat_total_cond3(:,:,toi),3));
            
            de_cond4 = squeeze(mean(blockDat_evoked_cond4(:,:,toi),3));
            dt_cond4 = squeeze(mean(blockDat_total_cond4(:,:,toi),3));
            
            
            
            % do forward model
            cond1Total_i = 1; cond2Total_i = 1; cond3Total_i=1; cond4Total_i=1;
            for i=1:nBlocks % loop through blocks, holding each out as the test set
                
                trnl = labels(blockNum~=i); % training labels
                tstl = labels(blockNum==i); % test labels
                
%                 % TOM ADDED
%                 cond1_tstl = [allPosBin(1).posBin];
%                 cond2_tstl = [allPosBin(2).posBin];
%                 cond3_tstl = [allPosBin(3).posBin];
%                 cond4_tstl = [allPosBin(4).posBin];
%                 
                
                
                cond1_tstl = cond1_labels(cond1_blockNum == i);
                cond2_tstl = cond2_labels(cond2_blockNum == i);
                cond3_tstl = cond3_labels(cond3_blockNum == i);
                cond4_tstl = cond4_labels(cond4_blockNum == i);
                
                
                %-----------------------------------------------------%
                % Analysis on Total Power                             %
                %-----------------------------------------------------%
                B1 = dt(blockNum~=i,:);    % training data
                B2 = dt(blockNum==i,:);    % test data
                C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                W = C1\B1;          % estimate weight matrix
                C2 = (W'\B2')';     % estimate channel responses
                
                C2_total(f,iter,t,i,:,:) = C2;
                
                % shift eegs to common center
                n2shift = ceil(size(C2,2)/2);
                for ii=1:size(C2,1)
                    [~, shiftInd] = min(abs(posBins-tstl(ii)));
                    C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                end
                
                tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                
                % Jordan - after training on both, test on conditions
                % seperately
                
                if i == 1 || i == 2
                    B2_cond1 = dt_cond1(cond1_blockNum==i,:);
                    C2_cond1 = (W'\B2_cond1')'; % test
                    C2_total_cond1(f,iter,t,cond1Total_i,:,:) = C2_cond1;
                    
                    %shift conditions
                    n2shift = ceil(size(C2_cond1,2)/2);
                    for ii=1:size(C2_cond1,1)
                        [~, shiftInd] = min(abs(posBins-cond1_tstl(ii)));
                        C2_cond1(ii,:) = wshift('1D', C2_cond1(ii,:), shiftInd-n2shift-1);
                    end
                    
                    tf_total_cond1(f,iter,t,cond1Total_i,:) = mean(C2_cond1,1);
                    cond1Total_i = cond1Total_i + 1;
                    
                elseif i == 3 || i == 4
                    
                    B2_cond2 = dt_cond2(cond2_blockNum==i,:);
                    C2_cond2 = (W'\B2_cond2')'; % test
                    C2_total_cond2(f,iter,t,cond2Total_i,:,:) = C2_cond2;
                    
                    
                    n2shift = ceil(size(C2_cond2,2)/2);
                    for ii=1:size(C2_cond2,1)
                        [~, shiftInd] = min(abs(posBins-cond2_tstl(ii)));
                        C2_cond2(ii,:) = wshift('1D', C2_cond2(ii,:), shiftInd-n2shift-1);
                    end
                    
                    
                    tf_total_cond2 (f,iter,t,cond2Total_i,:)  = mean(C2_cond2,1);
                    cond2Total_i = cond2Total_i + 1;
                    
                elseif i == 5 || i == 6
                    
                    B2_cond3 = dt_cond3(cond3_blockNum==i,:);
                    C2_cond3 = (W'\B2_cond3')'; % test
                    C2_total_cond3(f,iter,t,cond3Total_i,:,:) = C2_cond3;
                    
                    
                    n2shift = ceil(size(C2_cond3,2)/2);
                    for ii=1:size(C2_cond3,1)
                        [~, shiftInd] = min(abs(posBins-cond3_tstl(ii)));
                        C2_cond3(ii,:) = wshift('1D', C2_cond3(ii,:), shiftInd-n2shift-1);
                    end
                    
                    
                    tf_total_cond3 (f,iter,t,cond3Total_i,:)  = mean(C2_cond3,1);
                    cond3Total_i = cond3Total_i + 1;
                    
                elseif i == 7 || i == 8
                    
                    B2_cond4 = dt_cond4(cond4_blockNum==i,:);
                    C2_cond4 = (W'\B2_cond4')'; % test
                    C2_total_cond4(f,iter,t,cond4Total_i,:,:) = C2_cond4;
                    
                    
                    n2shift = ceil(size(C2_cond4,2)/2);
                    for ii=1:size(C2_cond4,1)
                        [~, shiftInd] = min(abs(posBins-cond4_tstl(ii)));
                        C2_cond4(ii,:) = wshift('1D', C2_cond4(ii,:), shiftInd-n2shift-1);
                    end
                    
                    
                    tf_total_cond4 (f,iter,t,cond4Total_i,:)  = mean(C2_cond4,1);
                    cond4Total_i = cond4Total_i + 1;
                    
                end
                
                %-----------------------------------------------------%
                
                
                
%                 %-----------------------------------------------------%
%                 % Analysis on Evoked Power                            %
%                 %-----------------------------------------------------%
%                 B1 = de(blockNum~=i,:);    % training data
%                 B2 = de(blockNum==i,:);    % test data
%                 C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                 W = C1\B1;          % estimate weight matrix
%                 
%                 
%                 
%                 %C2 = (W'\B2')';     % estimate channel responses
%                 
%                 %C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
%                 
%                 % Save weight matrix (Mary MacLean added)
%                 %WeightsTotal(:,:,i) = W;
%                 
%                 % shift eegs to common center
%                 n2shift = ceil(size(C2,2)/2);
%                 for ii=1:size(C2,1)
%                     [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                     C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                 end
%                 
%                 tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                 
%                 
%                 %-----------------------------------------------------%
%                 % Analysis on Total Power                             %
%                 %-----------------------------------------------------%
%                 B1 = dt(blockNum~=i,:);    % training data
%                 B2 = dt(blockNum==i,:);    % test data
%                 C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                 W = C1\B1;          % estimate weight matrix
%                 C2 = (W'\B2')';     % estimate channel responses
%                 
%                 C2_total(f,iter,t,i,:,:) = C2;
%                 
%                 % shift eegs to common center
%                 n2shift = ceil(size(C2,2)/2);
%                 for ii=1:size(C2,1)
%                     [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                     C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                 end
%                 
%                 tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                 %-----------------------------------------------------%
                
            end
            % Average weights across blocks add to matrix (Mary
            % MacLean, added)
            %%allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
        end
        
        %%% IEM END %%%
        
        
        
        
    end
    
end

% collapse and save TF data
em.tf_total_cond1 = squeeze(mean(mean(tf_total_cond1,2),4));
em.tf_total_cond2 = squeeze(mean(mean(tf_total_cond2,2),4));
em.tf_total_cond3 = squeeze(mean(mean(tf_total_cond3,2),4));
em.tf_total_cond4 = squeeze(mean(mean(tf_total_cond4,2),4));

% save data
save([iemDir '/' sprintf('sj%02d_fixed_IEM.mat',sn)],'em');




% 
% 
% 
%         %% TOM ADDED
%         % average over 1000 ITS + 3 BLOCKS to reduce size of saved file!
%         tf_evoked = squeeze(mean(mean(tf_evoked,2),4)); %tfs are nBins x ms 
%         tf_total = squeeze(mean(mean(tf_total,2),4));
%         
%         %Jordan Added
%         tf_evoked_rest = squeeze(mean(mean(tf_evoked_rest,2),4));
%         tf_evoked_low =  squeeze(mean(mean(tf_evoked_low,2),4));
%         
%         tf_total_rest = squeeze(mean(mean(tf_total_rest,2),4));
%         tf_total_low =  squeeze(mean(mean(tf_total_low,2),4));
%         
%         cd([root '/Analysis_Scripts'])
%         % save data
%         if processForIEM==1
%             fName = [dRoot, 'TrainBoth/',sprintf('sj%02d_TrainBoth_changeDect_accTrials', subs(s)),name];
%             em.C2.evoked = C2_evoked;
%             em.C2.total = C2_total;
%             em.tfs.evoked = tf_evoked;
%             em.tfs.total = tf_total;
%             em.tfs.totalW = allWeightsTotal;
%             em.nBlocks = nBlocks;
%             %save(fName,'em','-v7.3');
%             em.C2_rest.evoked = C2_evoked_rest;
%             em.C2_rest.total = C2_total_rest;
%             em.C2_low.evoked = C2_evoked_low;
%             em.C2_low.total = C2_total_low;
%             em.tfs_rest.evoked = tf_evoked_rest;
%             em.tfs_rest.total = tf_total_rest;
%             em.tfs_low.evoked = tf_evoked_low;
%             em.tfs_low.total = tf_total_low;
%             parsave(fName,em,minCnt,nElectrodes);
%         else
%             % save raw hilbert files
%             fName = [hRoot,sprintf('sj%02d_TrainBoth_changeDect_accTrials', subs(s)),hilbName];
%             eegInfo.chanLabels = eeg.chanLabels;
%             eegInfo.preTime = eeg.preTime;
%             eegInfo.postTime = eeg.postTime;
%             eegInfo.sampRate = eeg.sampRate;
%             eegInfo.posBin = posBin;
%             parsave(fName, eegs, eegInfo, eegBand)
%         end
%         
%     end
%     
%     %matlabpool close
%     
% end
%        
% 
% %%%%%%%%%%%%%%%%%%%
% 
% % NOT SAVING ANYTHING YET
% 
% %%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 






% 
% % loop through conditions
% for iCond=1:4
%     
%     % set posBin and minCnt for this cond
%     posBin = allPosBin(iCond).posBin;
%     
%     for iMin=1:4
%         minCntMat(iMin)=allConds(iMin).minCnt;
%     end
%     
%     minCnt = min(minCntMat)-1; % minus 1 to give all bins variability across iterations
%     
%     %%clear eegs
%     
%     
%     % loop through frequency bands (typically alpha, theta)
%     for f=1:size(freqs,1)
%         
%         clear data bandEEG dataFilt1
%         
%         % grab eegs
%         %eegs = allConds(iCond).eegs;
%         
%         % apply butterworth filter (3rd order, bandpass)
%         [z1,p1] = butter(3, [freqs(f,1), freqs(f,2)]./(eeg.sampRate/2),'bandpass');
%         data = double(eegs);
%         bandEEG = NaN(size(data,1),size(data,2),size(data,3));
%         for x = 1:size(data,1)
%             for y = 1:size(data,3)
%                 dataFilt1 = filtfilt(z1,p1,data(x,:,y));
%                 bandEEG(x,:,y) = dataFilt1;
%             end
%         end
%         
%         % apply hilbert tranform to bandpassed data
%         eegs = [];
%         for i=1:size(bandEEG,1) % chan loop
%             i
%             for j=1:size(bandEEG,3) % trial loop
%                 eegs(i,:,j) = hilbert(squeeze(bandEEG(i,:,j)));
%             end
%         end
%         
%         % create bandpassed data structure and save for sub\cond\freq
%         band.eeg = eegs;
%         band.freqs = freqs(f,:);
%         band.chanlocs = eeg.chanLabels;
%         band.srate = eeg.sampRate;
%         band.beh = allConds(iCond).beh;
%         band.times = eeg.times;
%         %save([bandpassedDir '/' sprintf('sj%02d_cond%02d_%s',sn,iCond,bands{f})],'band','-v7.3');
%         
%         % remove bad channels for IEM analyses
%         clear badChanIdx thisBadChan
%         badChanIdx = [];
%         cnt=0;
%         for i=1:length(eeg.badChannels)
%             thisBadChan=eeg.badChannels{i};
%             badChanIdx(i) = find(strcmp(thisBadChan,eeg.chanLabels)==1);
%         end
%         eegs(badChanIdx,:,:)=[];
%         
%         % generate evoked and total data mats
%         % permute data to(trials x elects x times for IEM script)
%         eegs = permute(eegs,[3,1,2]);
%         fdata_evoked = eegs;
%         fdata_total = abs(eegs).^2;
%         
%         % determine timepoints to run analysis on
%         tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time);
%         nTimes = length(tois); % index time points for analysis.
%         
%         % Loop through each iteration
%         for iter = 1:nIter
%             
%             disp(['nIteration :' num2str(iter) ' of ' num2str(nIter)])
%             
%             % assign trials to blocks - preallocate arrays
%             blocks = nan(size(posBin));
%             shuffBlocks = nan(size(posBin));
%             nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
%             
%             % shuffle trials
%             %shuffInd = randperm(nTrials)';
%             shuffInd = randperm(length(posBin))'; % create shuffle index
%             shuffBin = posBin(shuffInd); % shuffle trial order
%             
%             % take the 1st nPerBin x nBlocks trials for each position bin.
%             for bin = 1:nBins
%                 idx = find(shuffBin == bin); % get index for trials belonging to the current bin
%                 idx = idx(1:nPerBin*nBlocks); % drop excess trials
%                 x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
%             end
%             
%             % unshuffle block assignment
%             blocks(shuffInd) = shuffBlocks;
%             
%             % save block assignment
%             em.blocks(:,iter) = blocks; % block assignment
%             em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
%             
%             %-------------------------------------------------------------------------
%             
%             % Average data for each position bin across blocks
%             nElectrodes = size(fdata_total,2);
%             posBins = 1:nBins;
%             blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
%             blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
%             labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
%             blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
%             c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
%             bCnt = 1;
%             for ii = 1:nBins
%                 for iii = 1:nBlocks
%                     blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
%                     blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
%                     labels(bCnt) = ii;
%                     blockNum(bCnt) = iii;
%                     c(bCnt,:) = basisSet(ii,:);
%                     bCnt = bCnt+1;
%                 end
%             end
%             
%             %==========================================================
%             %{
%                 Run high temporal resolution IEM for "within" analysis only
%             %}
%             %==========================================================
%             
%             for t = 1:nSamps
%                 
%                 % grab data for timepoint t
%                 toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
%                 %find(toi==1)
%                 de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
%                 dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
%                 
%                 % Do forward model
%                 for i=1:nBlocks % loop through blocks, holding each out as the test set
%                     
%                     trnl = labels(blockNum~=i); % training labels
%                     tstl = labels(blockNum==i); % test labels
%                     
%                     %-----------------------------------------------------%
%                     % Analysis on Evoked Power                            %
%                     %-----------------------------------------------------%
%                     B1 = de(blockNum~=i,:);    % training data
%                     B2 = de(blockNum==i,:);    % test data
%                     C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                     W = C1\B1;          % estimate weight matrix
%                     C2 = (W'\B2')';     % estimate channel responses
%                     
%                     C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
%                     
%                     % Save weight matrix (Mary MacLean added)
%                     %WeightsTotal(:,:,i) = W;
%                     
%                     % shift eegs to common center
%                     n2shift = ceil(size(C2,2)/2);
%                     for ii=1:size(C2,1)
%                         [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                         C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                     end
%                     
%                     tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                     
%                     
%                     %-----------------------------------------------------%
%                     % Analysis on Total Power                             %
%                     %-----------------------------------------------------%
%                     B1 = dt(blockNum~=i,:);    % training data
%                     B2 = dt(blockNum==i,:);    % test data
%                     C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                     W = C1\B1;          % estimate weight matrix
%                     C2 = (W'\B2')';     % estimate channel responses
%                     
%                     C2_total(f,iter,t,i,:,:) = C2;
%                     
%                     % shift eegs to common center
%                     n2shift = ceil(size(C2,2)/2);
%                     for ii=1:size(C2,1)
%                         [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                         C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                     end
%                     
%                     tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                     %-----------------------------------------------------%
%                     
%                 end
%                 % Average weights across blocks add to matrix (Mary
%                 % MacLean, added)
%                 %%allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
%             end
%             
%             
%             %==========================================================
%             %{
%                 Run IEM with reduced SR for cross condition
%                 training/testing and save weights.
%                 Alpha band only!
%             %}
%             %==========================================================
%             
%             if f==1
%                 
%                 % REAL DATA
%                 
%                 % training loop
%                 %%tf_total=[];
%                 for trLoop=1:40 % divide 640 samples into 40 samples of 16 (62.5ms per sample)
%                     
%                     theseSamples = (16*trLoop)-15:(16*trLoop);
%                     trainData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
%                     
%                     % testing loop
%                     for teLoop=1:40
%                         
%                         theseSamples = (16*teLoop)-15:(16*teLoop);
%                         testData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
%                         
%                         % Do forward model
%                         for i=1:nBlocks % loop through blocks, holding each out as the test set
%                             
%                             trnl = labels(blockNum~=i); % training labels
%                             tstl = labels(blockNum==i); % test labels
%                             
%                             %-----------------------------------------------------%
%                             % Analysis on Total Power                             %
%                             %-----------------------------------------------------%
%                             B1 = trainData(blockNum~=i,:);    % training data
%                             B2 = testData(blockNum==i,:);    % test data
%                             C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                             W = C1\B1;          % estimate weight matrix
%                             C2 = (W'\B2')';     % estimate channel responses
%                             
%                             %C2_total(f,iter,t,i,:,:) = C2;
%                             
%                             % shift eegs to common center
%                             n2shift = ceil(size(C2,2)/2);
%                             for ii=1:size(C2,1)
%                                 [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                                 C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                             end
%                             
%                             %tf_total(iter,trLoop,teLoop,i,:) = mean(C2,1);
%                             
%                             % save B1, B2 and C1 for cross condition TT
%                             allB1(iter,trLoop,teLoop,i,:,:)=B1;
%                             allB2(iter,trLoop,teLoop,i,:,:)=B2;
%                             %allC1(iter,trLoop,teLoop,i,:,:)=C1;
%                             
%                         end
%                         
%                     end
%                 end
%             end
%             
%             
%             %==========================================================
%             %{
%                 Run high temporal IEM permuted analyis
%             %}
%             %==========================================================
%             
%             % Loop through permutations
%             for perm = 1:10;%nPerms
%                 tic % start timing permutation loop
%                 %fprintf('Permutation %d out of %d/n',perm,nPerms);
%                 
%                 %-----------------------------------------------------------------------------
%                 % Permute trial assignment within each block
%                 %-----------------------------------------------------------------------------
%                 permedPosBin = nan(size(posBin)); % preallocate permuted position bins vector
%                 for b = 1:nBlocks % for each block..
%                     pInd = randperm(em.nTrialsPerBlock); % create a permutation index
%                     permedBins(pInd) = posBin(blocks == b); % grab block b data and permute according data according to the index
%                     permedPosBin(blocks == b) = permedBins; % put permuted data into permedPosBin
%                     permInd(f,iter,perm,b,:) = pInd; % save the permutation (permInd is saved at end of the script)
%                 end
%                 
%                 %-----------------------------------------------------------------------------
%                 
%                 % Average data for each position bin across blocks
%                 posBins = 1:nBins;
%                 blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
%                 blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
%                 labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
%                 blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
%                 c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
%                 bCnt = 1;
%                 for ii = 1:nBins
%                     for iii = 1:nBlocks
%                         blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(permedPosBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
%                         blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1));
%                         labels(bCnt) = ii;
%                         blockNum(bCnt) = iii;
%                         c(bCnt,:) = basisSet(ii,:);
%                         bCnt = bCnt+1;
%                     end
%                 end
%                 
%                 for t = 1:nSamps
%                     
%                     % grab data for timepoint t
%                     toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
%                     de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
%                     dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
%                     
%                     % Do forward model
%                     tmpeC2 = nan(nBlocks,nBins,nChans); tmptC2 = tmpeC2; % for unshifted channel responses
%                     tmpeCR = nan(nBlocks,nChans); tmptCR = nan(nBlocks,nChans); % for shifted channel respones
%                     
%                     for i=1:nBlocks % loop through blocks, holding each out as the test set
%                         
%                         trnl = labels(blockNum~=i); % training labels
%                         tstl = labels(blockNum==i); % test labels
%                         
%                         %-----------------------------------------------------%
%                         % Analysis on Evoked Power                            %
%                         %-----------------------------------------------------%
%                         B1 = de(blockNum~=i,:);    % training data
%                         B2 = de(blockNum==i,:);    % test data
%                         C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                         W = C1\B1;          % estimate weight matrix
%                         C2 = (W'\B2')';     % estimate channel responses
%                         
%                         % tmpeC2(i,:,:) = C2;
%                         
%                         % shift eegs to common center
%                         n2shift = ceil(size(C2,2)/2);
%                         for ii=1:size(C2,1)
%                             [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                             C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                         end
%                         
%                         tmpeCR(i,:) = mean(C2); % average shifted channel responses
%                         
%                         %-----------------------------------------------------%
%                         % Analysis on Total Power                             %
%                         %-----------------------------------------------------%
%                         B1 = dt(blockNum~=i,:);    % training data
%                         B2 = dt(blockNum==i,:);    % test data
%                         C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                         W = C1\B1;          % estimate weight matrix
%                         C2 = (W'\B2')';     % estimate channel responses
%                         
%                         % tmptC2(i,:,:) = C2;
%                         
%                         % shift eegs to common center
%                         n2shift = ceil(size(C2,2)/2);
%                         for ii=1:size(C2,1)
%                             [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                             C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                         end
%                         
%                         tmptCR(i,:) = mean(C2); % averaged shifted channel responses
%                         
%                         %-----------------------------------------------------%
%                         
%                     end
%                     % save data to indexed matrix
%                     % C2_evoked(f,iter,perm,t,:,:) = mean(tmpeC2);
%                     % C2_total(f,iter,perm,t,:,:) = mean(tmptC2);
%                     tf_evoked_perm(f,iter,perm,t,:) = mean(tmpeCR);
%                     tf_total_perm(f,iter,perm,t,:) = mean(tmptCR);
%                 end
%                 toc
%             end
%             
%             %clear permInd permedBins
%             
%             
%             %==========================================================
%             %{
%                 Run PERMUTED IEM with reduced SR for cross condition
%                 training/testing and save weights.
%                 Alpha band only!
%             %}
%             %==========================================================
%             
%             % Loop through permutations
%             for perm = 1:5;%nPerms
%                 tic % start timing permutation loop
%                 %fprintf('Permutation %d out of %d/n',perm,nPerms);
%                 
%                 %-----------------------------------------------------------------------------
%                 % Permute trial assignment within each block
%                 %-----------------------------------------------------------------------------
%                 permedPosBin = nan(size(posBin)); % preallocate permuted position bins vector
%                 for b = 1:nBlocks % for each block..
%                     pInd = randperm(em.nTrialsPerBlock); % create a permutation index
%                     permedBins(pInd) = posBin(blocks == b); % grab block b data and permute according data according to the index
%                     permedPosBin(blocks == b) = permedBins; % put permuted data into permedPosBin
%                     permInd(f,iter,perm,b,:) = pInd; % save the permutation (permInd is saved at end of the script)
%                 end
%                 
%                 %-----------------------------------------------------------------------------
%                 
%                 % Average data for each position bin across blocks
%                 posBins = 1:nBins;
%                 blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
%                 blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
%                 labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
%                 blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
%                 c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
%                 bCnt = 1;
%                 for ii = 1:nBins
%                     for iii = 1:nBlocks
%                         blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(permedPosBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
%                         blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1));
%                         labels(bCnt) = ii;
%                         blockNum(bCnt) = iii;
%                         c(bCnt,:) = basisSet(ii,:);
%                         bCnt = bCnt+1;
%                     end
%                 end
%                 
%                 % training loop
%                 for trLoop=1:40 % divide 640 samples into 40 samples of 16 (62.5ms per sample)
%                     
%                     theseSamples = (16*trLoop)-15:(16*trLoop);
%                     trainData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
%                     
%                     % testing loop
%                     for teLoop=1:40
%                         
%                         theseSamples = (16*teLoop)-15:(16*teLoop);
%                         testData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
%                         
%                         % Do forward model
%                         for i=1:nBlocks % loop through blocks, holding each out as the test set
%                             
%                             trnl = labels(blockNum~=i); % training labels
%                             tstl = labels(blockNum==i); % test labels
%                             
%                             %-----------------------------------------------------%
%                             % Analysis on Total Power                             %
%                             %-----------------------------------------------------%
%                             B1 = trainData(blockNum~=i,:);    % training data
%                             B2 = testData(blockNum==i,:);    % test data
%                             C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                             W = C1\B1;          % estimate weight matrix
%                             C2 = (W'\B2')';     % estimate channel responses
%                             
%                             %C2_total(f,iter,t,i,:,:) = C2;
%                             
%                             % shift eegs to common center
%                             n2shift = ceil(size(C2,2)/2);
%                             for ii=1:size(C2,1)
%                                 [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                                 C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                             end
%                             
%                             %tf_total(perm,iter,trLoop,teLoop,i,:) = mean(C2,1);
%                             
%                             % save B1, B2 and C1 for cross condition TT
%                             allB1_perm(perm,iter,trLoop,teLoop,i,:,:)=B1;
%                             allB2_perm(perm,iter,trLoop,teLoop,i,:,:)=B2;
%                             %allC1_perm(perm,iter,trLoop,teLoop,i,:,:)=C1;
%                             
%                         end
%                     end
%                 end
%             end
%             
%             %clear permInd permedBins
%             
%         end
%         
%     end
%     
%     % average high temporal resolution IEMs over block and iter to
%     % reduce file size for save
%     tf_evoked = squeeze(mean(mean(tf_evoked,2),4));
%     tf_total = squeeze(mean(mean(tf_total,2),4));
%     
%     tf_evoked_perm = squeeze(mean(mean(tf_evoked_perm,2),3));
%     tf_total_perm = squeeze(mean(mean(tf_total_perm,2),3));
%     
%     % organize high temporal res data
%     em_within.tfs.evoked = tf_evoked;
%     em_within.tfs.total = tf_total;
%     
%     % organize low temporal res data for cross training/testing
%     em.tfs_cross_alpha.allB1=allB1;
%     em.tfs_cross_alpha.allB2=allB2;
%     
%     % save high temporal res permuted data
%     em_within.tfs_perm.evoked = tf_evoked_perm;
%     em_within.tfs_perm.total = tf_total_perm;
%     
%     % save low temporal res permuted data
%     em.tfs_cross_alpha_perm.allB1=allB1_perm;
%     em.tfs_cross_alpha_perm.allB2=allB2_perm;
%     
%     save([iemDir '/' sprintf('sj%02d_cond%02d_IEM.mat',sn,iCond)],'em','em_within','minCnt','nElectrodes','-v7.3')
%     
%     clear em_within minCnt nElectrodes tf_evoked tf_total tf_evoked_perm tf_total_perm allB1 allB2 allB1_perm allB2_perm allC1 permInd perm permedBins permedPosBin B1 B2
%     
%     em.blocks = [];
%     
% end

%end

return