function SpatialIEM_TrainBoth_Balance (subjects, bands)
%==========================================================================
%{
Purpose: run spatial encoding model on single bandpass e.g. 8-12 Hz

Original Author:
Joshua J. Foster
joshua.james.foster@gmail.com
University of Chicago
August 12, 2015
Modified by Tom Bullock
UCSB Attention Lab
%}
%==========================================================================

close all

subs = subjects;


nSubs = length(subs);

%% process for IEM (1) or HILBERT (2)
processForIEM=1;

% setup directories
root = '/home/garrett/WTF_Bike';
dRoot = [root '/Data/'];

% % setup directories
%
% %root = pwd;
% out = 'AnalysisScripts';
% dRoot = [root(1:end-length(out)),'Data/'];

if processForIEM==1
    %eRoot = [root(1:end-length(out)),'EEG/'];
    eRoot = [root '/EEG/'];
else
    eRoot = [root '/HILBERT/'];
end

%hRoot = [root(1:end-length(out)),'HILBERT/'];
hRoot = [root '/HILBERT/'];

%%%CLEAN UP%%%

% choose method for finding the min number of trials per location bin
% (1=per-condition, 2=across conditions)
findMinType = 2; %********************************JORDAN EDIT, was 1 now 2**************

% load matrix of min trials per loc bin per condition
if processForIEM == 1
    minLocsData = load([root '/Analysis_Scripts/Minimum_Location_Bin_Mat_AccTrials.mat']);
elseif processForIEM== 2
    minLocsData = load([root '/Analysis_Scripts/HILBERT_Minimum_Location_Bin_Mat_AccTrials.mat']);
end

for bandpassLoop=bands % loop through different bandpass analyses e.g. alpha, theta
    
    if bandpassLoop==1
        name = '_SpatialTF_ALPHA.mat'; % name of files to be saved for IEM
        hilbName = '_Hilbert_ALPHA.mat';
        thisFreq = {'Alpha'};
        freqBandpass = [8 12];
    else
        name = '_SpatialTF_THETA.mat';
        hilbName = '_Hilbert_THETA.mat';
        thisFreq = {'Theta'};
        freqBandpass = [4 7];
    end
    
    % Loop through participants
    %matlabpool open 72
    for s = 1:nSubs
        
        % had to move all the settings insdie the parfor to get it to work...
        em = [];
        EEG = [];
        eegInfo = [];
        WeightsTotal = [];
        allWeightsTotal = [];
        
        % parameters to set
        em.nChans = 8; % # of channels;
        em.nBins = em.nChans; % # of stimulus bins, was em.nChans; using 16 since doing both conditions. > 8 is exercise data
        em.nIter = 10; % # of iterations %%% WAS SET TO 10!!!
        em.nBlocks = 4; %3; % # of blocks for cross-validation
        em.frequencies = freqBandpass; % frequency bands to analyze
        em.bands = thisFreq;
        em.window = 4;
        em.Fs = 250; % WAS 250
        em.nElectrodes = 999; % get later
        em.time = -.5*1000:1000/em.Fs:1.9961*1000; %   -500:4:2000; % time points of interest
        
        % for brevity in analysis
        nChans = em.nChans;
        nBins = em.nBins;
        nIter = em.nIter;
        nBlocks = em.nBlocks;
        freqs = em.frequencies;
        times = em.time;
        nFreqs = size(em.frequencies,1);
        nElectrodes = em.nElectrodes;
        nSamps = length(em.time);
        Fs = em.Fs;
        
        % Specify basis set
        em.sinPower = 7;
        em.x = linspace(0, 2*pi-2*pi/nBins, nBins);
        em.cCenters = linspace(0, 2*pi-2*pi/nChans, nChans);
        em.cCenters = rad2deg(em.cCenters);
        pred = sin(0.5*em.x).^em.sinPower; % hypothetical channel responses
        pred = wshift('1D',pred,5); % shift the initial basis function
        basisSet = nan(nChans,nBins);
        for c = 1:nChans
            basisSet(c,:) = wshift('1D',pred,-c); % generate circularly shifted basis functions
        end
        em.basisSet = basisSet; % save basis set to data structure
        
        sn = subs(s);
        fprintf('Subject:\t%d\n',sn)
        
        % Grab data------------------------------------------------------------
        
        % Get position bin index from behavior file
        restfName = [dRoot, sprintf('sj%02d_exerCon01_changeDect_MixModel_wBias_accTrials.mat',subs(s))];
        exfName = [dRoot, sprintf('sj%02d_exerCon02_changeDect_MixModel_wBias_accTrials.mat',subs(s))];
        
        % load files
        rest_tmp = []; rest_beh = [];
        rest_tmp = load(restfName);
        rest_beh = rest_tmp.beh;
        
        ex_tmp = []; ex_beh = [];
        ex_tmp = load(exfName);
        ex_beh = ex_tmp.beh;
        
        beh.trial.posBin = [rest_beh.trial.posBin,ex_beh.trial.posBin]; %merge data
        
        %Jordan Added: create index to keep track of rest and exercise position trials
        trial_indx = [repmat(1,1,length(rest_beh.trial.posBin)), repmat(2,1,length(ex_beh.trial.posBin))];
        
        em.posBin = beh.trial.posBin'; % add to fm structure so it's saved
        posBin = em.posBin;
        
        % Get EEG data
        restfName = [eRoot, sprintf('sj%02d_exerCon01_changeDect_EEG_accTrials.mat',subs(s))];
        exfName = [eRoot, sprintf('sj%02d_exerCon02_changeDect_EEG_accTrials.mat',subs(s))];
        
        % load file
        rest_tmp = []; rest_eeg = [];
        rest_tmp = load(restfName);
        rest_eeg = rest_tmp.eeg;
        
        ex_tmp = []; ex_eeg = [];
        ex_tmp = load(exfName);
        ex_eeg = ex_tmp.eeg;
        
        %Jordan Added
        eeg.data = cat(1,rest_eeg.data,ex_eeg.data); %merge eeg trials along first dimension (trials x channels x timepoints)
        eeg.arf.artIndCleaned = [rest_eeg.arf.artIndCleaned, ex_eeg.arf.artIndCleaned];
        
        eeg.preTime = rest_eeg.preTime;
        eeg.postTime = rest_eeg.postTime;
        eeg.chanLabels = rest_eeg.chanLabels;
        eeg.sampRate = rest_eeg.sampRate;
        
        % get n channels (to save later with TF files)
        %%% nElects = size(eeg.chanLabels,2);
        
        eegs = eeg.data(:,:,:); % get scalp EEG (drop EOG electrodes)
        artInd = eeg.arf.artIndCleaned.'; % grab artifact rejection index
        tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time); nTimes = length(tois); % index time points for analysis.
        
        % %     %% TOM EDIT TO PROCESS WHOLE EPOCH (WTF DATA EPOCHED FROM -.5 to 2)
        % %     tois = ones(1,size(eegs,3)); nTimes = length(tois);
        
        % Remove rejected trials
        eegs = eegs(~artInd,:,:);
        posBin = posBin(~artInd);
        
        %Jordan Added
        trial_indx = trial_indx(~artInd);
        trial_indx = trial_indx';
        
        em.nTrials = length(posBin); 
        nTrials = em.nTrials; % # of good trials
        
        %----------------------------------------------------------------------
        
        % Preallocate Matrices
        tf_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nChans); tf_total = tf_evoked;
        C2_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nBins,nChans); C2_total = C2_evoked;
        em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments
        
        %Jordan added
        tf_evoked_rest = []; %nan(nFreqs,nIter,nSamps,nBlocks/2,nChans); 
        tf_evoked_low = tf_evoked_rest; 
        tf_total_rest = tf_evoked_rest; tf_total_low = tf_total_rest;
        C2_evoked_rest = [];%nan(nFreqs,nIter,nSamps,nBlocks/2,nBins,nChans);
        C2_total_rest = C2_evoked_rest;
        
        C2_evoked_low = C2_evoked_rest;
        C2_total_low = C2_evoked_low;
        
        % TOM ADDED (needed to convert to EEGLAB format to use new eeglab filter)
        EEG.data = permute(eegs,[2,3,1]); % converts to chans x times x trials
        EEG.srate = Fs;
        EEG.trials = size(EEG.data,3);
        EEG.nbchan = size(EEG.data,1);
        EEG.pnts = size(EEG.data,2);
        
        % Loop through each frequency
        for f = 1:nFreqs
            tic % start timing frequency loop
            fprintf('Frequency %d out of %d\n', f, nFreqs)
            
            %% get no. of electrodes
            nElectrodes = size(eeg.data,2);
            em.nElectrodes = nElectrodes;
            disp(['nElecrodes changed to :' num2str(nElectrodes)])
            
            %% BUTTERWORTH FILTER
            filterorder = 3;
            type = 'bandpass';
            [z1,p1] = butter(filterorder, [freqs(f,1), freqs(f,2)]./(EEG.srate/2),type);
            %freqz(z1,p1,[],250)
            data = double(EEG.data);
            tempEEG = NaN(size(data,1),EEG.pnts,size(data,3));
            for x = 1:size(data,1) % loop through chans
                for y = 1:size(data,3) % loop through trials
                    dataFilt1 = filtfilt(z1,p1,data(x,:,y)); % was filtfilt
                    tempEEG(x,:,y) = dataFilt1; % tymp = chans x times x trials
                end
            end
            
            eegBand = [];
            eegBand = tempEEG;
            
            %% apply hilbert to each channel and epoch in turn (this should be correct)
            eegs = [];
            for j=1:size(tempEEG,1) % chan loop
                for i=1:size(tempEEG,3) % trial loop
                    eegs(i,j,:) = hilbert(squeeze(tempEEG(j,:,i)));
                end
            end
            
            % eegs is trials x elects x times
            fdata_evoked = eegs;
            fdata_total = abs(eegs).^2;
            
            % Loop through each iteration
            for iter = 1:nIter
                
                disp(['nIteration :' num2str(iter) ' of ' num2str(nIter)])
                %--------------------------------------------------------------------------
                % Assign trials to blocks (such that trials per position are
                % equated within blocks)
                %--------------------------------------------------------------------------
                
                % preallocate arrays
                blocks = nan(size(posBin));
                shuffBlocks = nan(size(posBin));
                
                %Jordan Added
                blocks_rest = nan(size(posBin)); blocks_low = nan(size(posBin));
                shuffBlocks_rest = nan(size(posBin));
                shuffBlocks_low = nan(size(posBin));
                
                % count number of trials within each position bin
                %clear binCnt
                binCnt = [];
                for bin = 1:nBins
                    binCnt(bin) = sum(posBin == bin);
                end
                
                % choose method of determining the location bin with the min
                % number of trials (1= per condition, 2= across conditions)
                if findMinType==1
                    minCnt = min(binCnt); % # of trials for position bin with fewest trials
                elseif findMinType==2
                    %if minCnt is based across all four conditions.  We get
                    %the values from "Minimum_Location_Bin_MAT.mat" created
                    %earlier
                    if ismember(sn,1)
                        minCnt=minLocsData.minLocBinAllConds(1);
                    elseif ismember(sn,2)
                        minCnt=minLocsData.minLocBinAllConds(2);
                    elseif ismember(sn,3)
                        minCnt=minLocsData.minLocBinAllConds(3);
                    elseif ismember(sn,4)
                        minCnt=minLocsData.minLocBinAllConds(4);
                    elseif ismember(sn,5)
                        minCnt=minLocsData.minLocBinAllConds(5);
                    elseif ismember(sn,6)
                        minCnt=minLocsData.minLocBinAllConds(6);
                    elseif ismember(sn,7)
                        minCnt=minLocsData.minLocBinAllConds(7);
                    elseif ismember(sn,8)
                        minCnt=minLocsData.minLocBinAllConds(8);
                    elseif ismember(sn,10)
                        minCnt=minLocsData.minLocBinAllConds(9);
                    elseif ismember(sn,11)
                        minCnt=minLocsData.minLocBinAllConds(10);
                    elseif ismember(sn,12)
                        minCnt=minLocsData.minLocBinAllConds(11);
                    elseif ismember(sn,13)
                        minCnt=minLocsData.minLocBinAllConds(12);
                    elseif ismember(sn,14)
                        minCnt=minLocsData.minLocBinAllConds(13);
                    elseif ismember(sn,15)
                        minCnt=minLocsData.minLocBinAllConds(14);
                    elseif ismember(sn,16)
                        minCnt=minLocsData.minLocBinAllConds(15);
                    elseif ismember(sn,17)
                        minCnt=minLocsData.minLocBinAllConds(16);
                    elseif ismember(sn,18)
                        minCnt=minLocsData.minLocBinAllConds(17);
                    elseif ismember(sn,19)
                        minCnt=minLocsData.minLocBinAllConds(18);
                    elseif ismember(sn,20)
                        minCnt=minLocsData.minLocBinAllConds(19);
                    elseif ismember(sn,21)
                        minCnt=minLocsData.minLocBinAllConds(20);
                    elseif ismember(sn,22)
                        minCnt=minLocsData.minLocBinAllConds(21);
                    elseif ismember(sn,23)
                        minCnt=minLocsData.minLocBinAllConds(22);
                    elseif ismember(sn,24)
                        minCnt=minLocsData.minLocBinAllConds(23);
                    elseif ismember(sn,25)
                        minCnt=minLocsData.minLocBinAllConds(24);
                    elseif ismember(sn,26)
                        minCnt=minLocsData.minLocBinAllConds(25);
                    elseif ismember(sn,27)
                        minCnt=minLocsData.minLocBinAllConds(26);
                    elseif ismember(sn,28)
                        minCnt=minLocsData.minLocBinAllConds(27);
                    elseif ismember(sn,29)
                        minCnt=minLocsData.minLocBinAllConds(28);
                    elseif ismember(sn,30)
                        minCnt=minLocsData.minLocBinAllConds(29);
                    elseif ismember(sn,31)
                        minCnt=minLocsData.minLocBinAllConds(30);
                    elseif ismember(sn,32)
                        minCnt=minLocsData.minLocBinAllConds(31);
                    elseif ismember(sn,33)
                        minCnt=minLocsData.minLocBinAllConds(32);
                    elseif ismember(sn,34)
                        minCnt=minLocsData.minLocBinAllConds(33);
                    elseif ismember(sn,35)
                        minCnt=minLocsData.minLocBinAllConds(34);
                    end
                end
                
                % reduce minCnt by 1 trial to ensure that ALL location bins
                % have some degree of trial randomization per IEM iteration
                % (without this, the loc bin with the minimum no. of trials
                % will not be randomized at all).
                minCnt = minCnt-1;
                
                nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
                
                % shuffle trials
                shuffInd = randperm(nTrials)'; % create shuffle index
                shuffBin = posBin(shuffInd); % shuffle trial order
                
                shuffTrial_indx = trial_indx(shuffInd);
                
                % take the 1st nPerBin x nBlocks trials for each position bin.
                for bin = 1:nBins
                    idx = find(shuffBin == bin); % get index for trials belonging to the current bin
                    idx = idx(1:nPerBin*nBlocks); % drop excess trials
                    %x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
                    

                    % Jordan Added
                    idx_rest = find(shuffBin == bin & shuffTrial_indx == 1);
                    idx_rest = idx_rest(1:nPerBin*nBlocks);
                    
                    idx_low = find(shuffBin == bin & shuffTrial_indx == 2);
                    idx_low = idx_low(1:nPerBin*nBlocks);
                    
                    con_idx = shuffTrial_indx(idx); 
                    con_idx = con_idx(1:nPerBin*nBlocks)';
                    
                    % Try to have equal number of trials for each bin from
                    % each condition during training
                    
%                     block_numTrials = max(1:nPerBin*nBlocks);
%                     
%                     % If nPerBin not evenly divisible, then randomly decide
%                     % which condition gets the extra trial
%                     if rem(block_numTrials,2)
%                         first_half = round(block_numTrials/2);
%                         second_half = block_numTrials - first_half;
%                         con_trialNum = [first_half, second_half];
%                         %shuffle to avoid bias
%                         shuff_splitNum = con_trialNum(randperm(length(con_trialNum)));
%                     else
%                         shuff_splitNum = repmat(block_numTrials/2,1,2);
%                     end
%                     
%                     idx_restTrain = idx_rest(1:shuff_splitNum(1));
%                     idx_lowTrain = idx_low(1:shuff_splitNum(2));
                    fixed_TrainIdx = [idx_rest;idx_low];
                    
                    
                    % assign rest trials to first two blocks and exercise
                    % trials to second two blocks
                    if rem(length(idx_rest),2)
                        first_blockN = round(length(idx_rest)/2);
                        second_blockN = length(idx_rest) - first_blockN;
                        x_rest = [ones(1,first_blockN),repmat(2,1,second_blockN)];
                    else
                        x_rest = repmat(1:2,1,length(idx_rest)/2);
                    end
                    
                    if rem(length(idx_low),2)
                        first_blockN = round(length(idx_low)/2);
                        second_blockN = length(idx_low) - first_blockN;
                        x_low = [repmat(3,1,first_blockN),repmat(4,1,second_blockN)];
                    else
                        x_low = repmat(3:4,1,length(idx_low)/2);
                    end
                    
                    fixed_x = [x_rest, x_low];
                    shuffBlocks(fixed_TrainIdx) = fixed_x;
                    
                    % Have test blocks that contain only trials from one
                    % condition. Need to make sure not testing on trials
                    % previously trained on
%                     x = repmat(1:nBlocks',nPerBin,1);
%                     idx_restTest = idx_rest(1:nPerBin*nBlocks);
%                     idx_lowTest = idx_low(1:nPerBin*nBlocks);
%                     shuffBlocks_rest(idx_restTest) = x; 
%                     shuffBlocks_low(idx_lowTest) = x;
                end
                
                % unshuffle block assignment
                blocks(shuffInd) = shuffBlocks;
                
                %Jordan Added
                
                block_conMat = [blocks, shuffTrial_indx];
                blockCon_trialNum = [];
                for iBlock = 1:nBlocks
                    for iCon = 1:2
                        blockCon_trialNum(iBlock, iCon) = sum(block_conMat(:,1) == iBlock & block_conMat(:,2) == iCon);
                    end
                end
                
                
%                 blocks_rest(shuffInd) = blocks(blocks == 1 | blocks == 2);
%                 blocks_low(shuffInd) = shuffBlocks_low;
               
                % save block assignment
                em.blocks(:,iter) = blocks; % block assignment
                em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
                
                %Jordan Added
                %em.blocksCon_Mat(:,:,iter) = block_conMat; %block trials condition assignment
                em.trial_Indx(:,iter) = trial_indx; %save condition assignment for each trial
%                 em.blocks_rest(:,iter) = blocks_rest;
%                 em.blocks_low(:,iter) = blocks_low;
                
                %-------------------------------------------------------------------------
                
                % Average data for each position bin across blocks
                posBins = 1:nBins;
                blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
                blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
                
                blockDat_evoked_Rest = [];%nan(nBins*nBlocks/2,nElectrodes,nSamps);
                blockDat_total_Rest = nan(size(blockDat_evoked_Rest));
                
                blockDat_evoked_Low = nan(size(blockDat_evoked_Rest));
                blockDat_total_Low = nan(size(blockDat_evoked_Low));
                
                labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
                blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
                c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
                bCnt = 1; 
                
                % Jordan added
                rest_labels = [];%nan(nBins*nBlocks/2,1);
                rest_blockNum = nan(size(rest_labels));
                rest_c =[];%nan(nBins*nBlocks/2,nChans);
                
                ex_labels = nan(size(rest_labels));
                ex_blockNum = nan(size(ex_labels));
                ex_c = nan(size(rest_c));
                
                rest_bCnt = 1;
                ex_bCnt = 1;
                for ii = 1:nBins
                    for iii = 1:nBlocks
                        blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                        blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
                        
                        
                        % Seperate condition data
                        if (iii == 1 || iii == 2) % Rest condition
                            blockDat_evoked_Rest(rest_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total_Rest(rest_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                            
                            
                            % Con labels for testing
                            rest_labels(rest_bCnt) = ii;
                            rest_blockNum(rest_bCnt) = iii;
                            rest_c(rest_bCnt,:) = basisSet(ii,:);
                            
                            rest_bCnt = rest_bCnt + 1;
                        elseif (iii == 3 || iii == 4) % Exercise condition
                            blockDat_evoked_Low(ex_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total_Low(ex_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                        
                            ex_labels(ex_bCnt) = ii;
                            ex_blockNum(ex_bCnt) = iii;
                            ex_c(ex_bCnt,:) = basisSet(ii,:);
                            
                            ex_bCnt = ex_bCnt + 1;
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
                
                for t = 1:nSamps
                    
                    % grab data for timepoint t
                    toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
                    de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
                    dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
                    
                    %Jordan added
                    de_rest = squeeze(mean(blockDat_evoked_Rest(:,:,toi),3));
                    dt_rest = squeeze(mean(blockDat_total_Rest(:,:,toi),3));
                    
                    de_low = squeeze(mean(blockDat_evoked_Low(:,:,toi),3));
                    dt_low = squeeze(mean(blockDat_total_Low(:,:,toi),3));
                    
                    
                    % Do forward model
                    restEvoked_i = 1; exEvoked_i = 1; restTotal_i = 1; exTotal_i = 1;
                    for i=1:nBlocks % loop through blocks, holding each out as the test set
                        
                        trnl = labels(blockNum~=i); % training labels
                        tstl = labels(blockNum==i); % test labels
                        
                        
                        rest_tstl = rest_labels(rest_blockNum == i);
                        ex_tstl = ex_labels(ex_blockNum == i);
                        
                        %-----------------------------------------------------%
                        % Analysis on Evoked Power                            %
                        %-----------------------------------------------------%
                        B1 = de(blockNum~=i,:);    % training data
                        B2 = de(blockNum==i,:);    % test data
                        C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                        W = C1\B1;          % estimate weight matrix
                        C2 = (W'\B2')';     % estimate channel responses
                        
                        C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
                        
                        % Save weight matrix (Mary MacLean added)
                        WeightsTotal(:,:,i) = W;
                        
                        % shift eegs to common center
                        n2shift = ceil(size(C2,2)/2);
                        for ii=1:size(C2,1)
                            [~, shiftInd] = min(abs(posBins-tstl(ii)));
                            C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                        end
                        
                        tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                        
                        % Jordan - after training on both, test on conditions
                        % seperately
                        if i == 1 || i == 2
                            B2_rest = de_rest(rest_blockNum==i,:);
                            C2_rest = (W'\B2_rest')'; % test
                            C2_evoked_rest(f,iter,t,restEvoked_i,:,:) = C2_rest;
                            
                            %shift
                            n2shift = ceil(size(C2_rest,2)/2);
                            for ii=1:size(C2_rest,1)
                                [~, shiftInd] = min(abs(posBins-rest_tstl(ii)));
                                C2_rest(ii,:) = wshift('1D', C2_rest(ii,:), shiftInd-n2shift-1);
                            end
                            
                            tf_evoked_rest(f,iter,t,restEvoked_i,:) = mean(C2_rest,1);
                            
                            restEvoked_i = restEvoked_i + 1;
                        elseif i == 3 || i == 4
                            B2_low = de_low(ex_blockNum==i,:);
                            C2_low = (W'\B2_low')'; % test
                            C2_evoked_low(f,iter,t,exEvoked_i,:,:) = C2_low;
                            
                            n2shift = ceil(size(C2_low,2)/2);
                            for ii=1:size(C2_low,1)
                                [~, shiftInd] = min(abs(posBins-ex_tstl(ii)));
                                C2_low(ii,:) = wshift('1D', C2_low(ii,:), shiftInd-n2shift-1);
                            end
                            
                            
                            tf_evoked_low(f,iter,t,i,:) = mean(C2_low,1);
                            exEvoked_i = exEvoked_i + 1;
                        end
                        
                        
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
                            B2_rest = dt_rest(rest_blockNum==i,:);
                            C2_rest = (W'\B2_rest')'; % test
                            C2_total_rest(f,iter,t,restTotal_i,:,:) = C2_rest;
                            
                            %shift conditions
                            n2shift = ceil(size(C2_rest,2)/2);
                            for ii=1:size(C2_rest,1)
                                [~, shiftInd] = min(abs(posBins-rest_tstl(ii)));
                                C2_rest(ii,:) = wshift('1D', C2_rest(ii,:), shiftInd-n2shift-1);
                            end
                            
                            tf_total_rest(f,iter,t,restTotal_i,:) = mean(C2_rest,1);
                            restTotal_i = restTotal_i + 1;
                            
                        elseif i == 3 || i == 4
                            
                            B2_low = dt_low(ex_blockNum==i,:);
                            C2_low = (W'\B2_low')'; % test
                            C2_total_low(f,iter,t,exTotal_i,:,:) = C2_low;
                            
                            
                            n2shift = ceil(size(C2_low,2)/2);
                            for ii=1:size(C2_low,1)
                                [~, shiftInd] = min(abs(posBins-ex_tstl(ii)));
                                C2_low(ii,:) = wshift('1D', C2_low(ii,:), shiftInd-n2shift-1);
                            end
                            
                            
                            tf_total_low (f,iter,t,exTotal_i,:)  = mean(C2_low,1);
                        end
                       
                        %-----------------------------------------------------%
                        
                    end
                    % Average weights across blocks add to matrix (Mary
                    % MacLean, added)
                    allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
                end
            end
            toc % stop timing the frequency loop
        end
        
        %% TOM ADDED
        % average over 1000 ITS + 3 BLOCKS to reduce size of saved file!
        tf_evoked = squeeze(mean(mean(tf_evoked,2),4)); %tfs are nBins x ms 
        tf_total = squeeze(mean(mean(tf_total,2),4));
        
        %Jordan Added
        tf_evoked_rest = squeeze(mean(mean(tf_evoked_rest,2),4));
        tf_evoked_low =  squeeze(mean(mean(tf_evoked_low,2),4));
        
        tf_total_rest = squeeze(mean(mean(tf_total_rest,2),4));
        tf_total_low =  squeeze(mean(mean(tf_total_low,2),4));
        
        cd([root '/Analysis_Scripts'])
        % save data
        if processForIEM==1
            fName = [dRoot, 'TrainBoth/',sprintf('sj%02d_TrainBoth_changeDect_accTrials', subs(s)),name];
            em.C2.evoked = C2_evoked;
            em.C2.total = C2_total;
            em.tfs.evoked = tf_evoked;
            em.tfs.total = tf_total;
            em.tfs.totalW = allWeightsTotal;
            em.nBlocks = nBlocks;
            %save(fName,'em','-v7.3');
            em.C2_rest.evoked = C2_evoked_rest;
            em.C2_rest.total = C2_total_rest;
            em.C2_low.evoked = C2_evoked_low;
            em.C2_low.total = C2_total_low;
            em.tfs_rest.evoked = tf_evoked_rest;
            em.tfs_rest.total = tf_total_rest;
            em.tfs_low.evoked = tf_evoked_low;
            em.tfs_low.total = tf_total_low;
            parsave(fName,em,minCnt,nElectrodes);
        else
            % save raw hilbert files
            fName = [hRoot,sprintf('sj%02d_TrainBoth_changeDect_accTrials', subs(s)),hilbName];
            eegInfo.chanLabels = eeg.chanLabels;
            eegInfo.preTime = eeg.preTime;
            eegInfo.postTime = eeg.postTime;
            eegInfo.sampRate = eeg.sampRate;
            eegInfo.posBin = posBin;
            parsave(fName, eegs, eegInfo, eegBand)
        end
        
    end
    
    %matlabpool close
    
end

%sendEmailToMe('SPATIAL IEM SCRIPT FINISHED PROCESSING!!')

clear all
close all
end