%{
Spatial_IEM
Original Author: Joshua Foster (heavily edited by Tom Bullock, UCSB Attention Lab)
Date:11.15.19

Load preprocessed EEG and Beh data.
Merge across sessions and split into condition
Get time-frequency data
Either 1) Run IEM or 2) just get bandpassed data ...

%}

clear
close all

% set dirs
rDir = '/Users/tombullock/Documents/Psychology/WTF_EYE';
eegDir = [rDir '/' 'EEG_Prepro2_Avg_Baseline'];
bandpassedDir = [rDir '/' 'EEG_Bandpassed'];
iemDir = [rDir '/' 'IEM_Results_TT_Within' ];
% select subjects
subs = 5;

% loop through subs
for s=1:length(subs)
    sn=subs(s);
    
    % import IEM settings
    em = IEM_Settings;
    nChans = em.nChans;
    nBins = em.nBins;
    nIter = em.nIter;
    nBlocks = em.nBlocks;
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
        for f=1:size(freqs,1)
            
            clear data bandEEG dataFilt1
            
            % grab eegs
            eegs = allConds(iCond).eegs;

            % apply butterworth filter (3rd order, bandpass)
            [z1,p1] = butter(3, [freqs(f,1), freqs(f,2)]./(eeg.sampRate/2),'bandpass');
            data = eegs;
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
           
            % create bandpassed data structure and save for sub/cond/freq
            band.eeg = eegs;
            band.freqs = freqs(f,:);
            band.chanlocs = eeg.chanLabels;
            band.srate = eeg.sampRate;
            band.beh = allConds(iCond).beh;
            band.times = eeg.times;
            save([bandpassedDir '/' sprintf('sj%02d_cond%02d_%s',sn,iCond,bands{f})],'band','-v7.3');
            
            % remove bad channels for IEM analyses
            clear badChanIdx thisBadChan
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
                
                % take the 1st nPerBin x nBlocks trials for each position bin.
                for bin = 1:nBins
                    idx = find(shuffBin == bin); % get index for trials belonging to the current bin
                    idx = idx(1:nPerBin*nBlocks); % drop excess trials
                    x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
                end
                
                % unshuffle block assignment
                blocks(shuffInd) = shuffBlocks;
                
                % save block assignment
                em.blocks(:,iter) = blocks; % block assignment
                em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
                
                %-------------------------------------------------------------------------
                
                % Average data for each position bin across blocks
                nElectrodes = size(fdata_total,2);
                posBins = 1:nBins;
                blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
                blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
                labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
                blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
                c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
                bCnt = 1;
                for ii = 1:nBins
                    for iii = 1:nBlocks
                        blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                        blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
                        labels(bCnt) = ii;
                        blockNum(bCnt) = iii;
                        c(bCnt,:) = basisSet(ii,:);
                        bCnt = bCnt+1;
                    end
                end
                
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
                    for i=1:nBlocks % loop through blocks, holding each out as the test set
                        
                        trnl = labels(blockNum~=i); % training labels
                        tstl = labels(blockNum==i); % test labels
                        
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
                        %WeightsTotal(:,:,i) = W;
                        
                        % shift eegs to common center
                        n2shift = ceil(size(C2,2)/2);
                        for ii=1:size(C2,1)
                            [~, shiftInd] = min(abs(posBins-tstl(ii)));
                            C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                        end
                        
                        tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                        
                        
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
                        %-----------------------------------------------------%
                        
                    end
                    % Average weights across blocks add to matrix (Mary
                    % MacLean, added)
                    %%allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
                end
                
                
                %==========================================================
                %{ 
                Run IEM with reduced SR for cross condition
                training/testing and save weights.
                Alpha band only!
                %}
                %==========================================================
                
                if f==1
                    
                    % REAL DATA
                    
                    % training loop
                    %%tf_total=[];
                    for trLoop=1:40 % divide 640 samples into 40 samples of 16 (62.5ms per sample)
                        
                        theseSamples = (16*trLoop)-15:(16*trLoop);
                        trainData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
                        
                        % testing loop
                        for teLoop=1:40
                            
                            theseSamples = (16*teLoop)-15:(16*teLoop);
                            testData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
                            
                            % Do forward model
                            for i=1:nBlocks % loop through blocks, holding each out as the test set
                                
                                trnl = labels(blockNum~=i); % training labels
                                tstl = labels(blockNum==i); % test labels
                                
                                %-----------------------------------------------------%
                                % Analysis on Total Power                             %
                                %-----------------------------------------------------%
                                B1 = trainData(blockNum~=i,:);    % training data
                                B2 = testData(blockNum==i,:);    % test data
                                C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                                W = C1\B1;          % estimate weight matrix
                                C2 = (W'\B2')';     % estimate channel responses
                                
                                %C2_total(f,iter,t,i,:,:) = C2;
                                
                                % shift eegs to common center
                                n2shift = ceil(size(C2,2)/2);
                                for ii=1:size(C2,1)
                                    [~, shiftInd] = min(abs(posBins-tstl(ii)));
                                    C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                                end
                                
                                %tf_total(iter,trLoop,teLoop,i,:) = mean(C2,1);
                                
                                % save B1, B2 and C1 for cross condition TT
                                allB1(iter,trLoop,teLoop,i,:,:)=B1;
                                allB2(iter,trLoop,teLoop,i,:,:)=B2;
                                %allC1(iter,trLoop,teLoop,i,:,:)=C1;
                                
                            end
                            
                        end
                    end
                end
                
                
                %==========================================================
                %{
                Run high temporal IEM permuted analyis
                %}
                %==========================================================
                
                % Loop through permutations
                nPerms=10;
                for perm = 1:nPerms
                    tic % start timing permutation loop
                    %fprintf('Permutation %d out of %d\n',perm,nPerms);
                    
                    %-----------------------------------------------------------------------------
                    % Permute trial assignment within each block
                    %-----------------------------------------------------------------------------
                    permedPosBin = nan(size(posBin)); % preallocate permuted position bins vector
                    for b = 1:nBlocks % for each block..
                        pInd = randperm(em.nTrialsPerBlock); % create a permutation index
                        permedBins(pInd) = posBin(blocks == b); % grab block b data and permute according data according to the index
                        permedPosBin(blocks == b) = permedBins; % put permuted data into permedPosBin
                        permInd(f,iter,perm,b,:) = pInd; % save the permutation (permInd is saved at end of the script)
                    end
                    
                    %-----------------------------------------------------------------------------
                    
                    % Average data for each position bin across blocks
                    posBins = 1:nBins;
                    blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
                    blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
                    labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
                    blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
                    c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
                    bCnt = 1;
                    for ii = 1:nBins
                        for iii = 1:nBlocks
                            blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(permedPosBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1));
                            labels(bCnt) = ii;
                            blockNum(bCnt) = iii;
                            c(bCnt,:) = basisSet(ii,:);
                            bCnt = bCnt+1;
                        end
                    end
                    
                    for t = 1:nSamps
                        
                        % grab data for timepoint t
                        toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
                        de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
                        dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
                        
                        % Do forward model
                        tmpeC2 = nan(nBlocks,nBins,nChans); tmptC2 = tmpeC2; % for unshifted channel responses
                        tmpeCR = nan(nBlocks,nChans); tmptCR = nan(nBlocks,nChans); % for shifted channel respones
                        
                        for i=1:nBlocks % loop through blocks, holding each out as the test set
                            
                            trnl = labels(blockNum~=i); % training labels
                            tstl = labels(blockNum==i); % test labels
                            
                            %-----------------------------------------------------%
                            % Analysis on Evoked Power                            %
                            %-----------------------------------------------------%
                            B1 = de(blockNum~=i,:);    % training data
                            B2 = de(blockNum==i,:);    % test data
                            C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                            W = C1\B1;          % estimate weight matrix
                            C2 = (W'\B2')';     % estimate channel responses
                            
                            % tmpeC2(i,:,:) = C2;
                            
                            % shift eegs to common center
                            n2shift = ceil(size(C2,2)/2);
                            for ii=1:size(C2,1)
                                [~, shiftInd] = min(abs(posBins-tstl(ii)));
                                C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                            end
                            
                            tmpeCR(i,:) = mean(C2); % average shifted channel responses
                            
                            %-----------------------------------------------------%
                            % Analysis on Total Power                             %
                            %-----------------------------------------------------%
                            B1 = dt(blockNum~=i,:);    % training data
                            B2 = dt(blockNum==i,:);    % test data
                            C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                            W = C1\B1;          % estimate weight matrix
                            C2 = (W'\B2')';     % estimate channel responses
                            
                            % tmptC2(i,:,:) = C2;
                            
                            % shift eegs to common center
                            n2shift = ceil(size(C2,2)/2);
                            for ii=1:size(C2,1)
                                [~, shiftInd] = min(abs(posBins-tstl(ii)));
                                C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                            end
                            
                            tmptCR(i,:) = mean(C2); % averaged shifted channel responses
                            
                            %-----------------------------------------------------%
                            
                        end
                        % save data to indexed matrix
                        % C2_evoked(f,iter,perm,t,:,:) = mean(tmpeC2);
                        % C2_total(f,iter,perm,t,:,:) = mean(tmptC2);
                        tf_evoked_perm(f,iter,perm,t,:) = mean(tmpeCR);
                        tf_total_perm(f,iter,perm,t,:) = mean(tmptCR);
                    end
                    toc
                end
                
                %clear permInd permedBins
                
                
                %==========================================================
                %{
                Run PERMUTED IEM with reduced SR for cross condition
                training/testing and save weights.
                Alpha band only!
                %}
                %==========================================================
                
                % Loop through permutations
                for perm = 1:nPerms
                    tic % start timing permutation loop
                    %fprintf('Permutation %d out of %d\n',perm,nPerms);
                    
                    %-----------------------------------------------------------------------------
                    % Permute trial assignment within each block
                    %-----------------------------------------------------------------------------
                    permedPosBin = nan(size(posBin)); % preallocate permuted position bins vector
                    for b = 1:nBlocks % for each block..
                        pInd = randperm(em.nTrialsPerBlock); % create a permutation index
                        permedBins(pInd) = posBin(blocks == b); % grab block b data and permute according data according to the index
                        permedPosBin(blocks == b) = permedBins; % put permuted data into permedPosBin
                        permInd(f,iter,perm,b,:) = pInd; % save the permutation (permInd is saved at end of the script)
                    end
                    
                    %-----------------------------------------------------------------------------
                    
                    % Average data for each position bin across blocks
                    posBins = 1:nBins;
                    blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
                    blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
                    labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
                    blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
                    c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
                    bCnt = 1;
                    for ii = 1:nBins
                        for iii = 1:nBlocks
                            blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(permedPosBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1));
                            labels(bCnt) = ii;
                            blockNum(bCnt) = iii;
                            c(bCnt,:) = basisSet(ii,:);
                            bCnt = bCnt+1;
                        end
                    end
                    
                    % training loop
                    for trLoop=1:40 % divide 640 samples into 40 samples of 16 (62.5ms per sample)
                        
                        theseSamples = (16*trLoop)-15:(16*trLoop);
                        trainData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
                        
                        % testing loop
                        for teLoop=1:40
                            
                            theseSamples = (16*teLoop)-15:(16*teLoop);
                            testData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
                            
                            % Do forward model
                            for i=1:nBlocks % loop through blocks, holding each out as the test set
                                
                                trnl = labels(blockNum~=i); % training labels
                                tstl = labels(blockNum==i); % test labels
                                
                                %-----------------------------------------------------%
                                % Analysis on Total Power                             %
                                %-----------------------------------------------------%
                                B1 = trainData(blockNum~=i,:);    % training data
                                B2 = testData(blockNum==i,:);    % test data
                                C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                                W = C1\B1;          % estimate weight matrix
                                C2 = (W'\B2')';     % estimate channel responses
                                
                                %C2_total(f,iter,t,i,:,:) = C2;
                                
                                % shift eegs to common center
                                n2shift = ceil(size(C2,2)/2);
                                for ii=1:size(C2,1)
                                    [~, shiftInd] = min(abs(posBins-tstl(ii)));
                                    C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                                end
                                
                                %tf_total(perm,iter,trLoop,teLoop,i,:) = mean(C2,1);
                                
                                % save B1, B2 and C1 for cross condition TT
                                allB1_perm(perm,iter,trLoop,teLoop,i,:,:)=B1;
                                allB2_perm(perm,iter,trLoop,teLoop,i,:,:)=B2;
                                %allC1_perm(perm,iter,trLoop,teLoop,i,:,:)=C1;
                                
                            end
                        end
                    end
                end
                
                %clear permInd permedBins

            end
            
        end
        
        % average high temporal resolution IEMs over block and iter to
        % reduce file size for save
        tf_evoked = squeeze(mean(mean(tf_evoked,2),4));
        tf_total = squeeze(mean(mean(tf_total,2),4));
        
        tf_evoked_perm = squeeze(mean(mean(tf_evoked_perm,2),3));
        tf_total_perm = squeeze(mean(mean(tf_total_perm,2),3));
        
        % organize high temporal res data
        em_within.tfs.evoked = tf_evoked;
        em_within.tfs.total = tf_total;
        
        % organize low temporal res data for cross training/testing
        em.tfs_cross_alpha.allB1=allB1;
        em.tfs_cross_alpha.allB2=allB2;
 
        % save high temporal res permuted data
        em_within.tfs_perm.evoked = tf_evoked_perm;
        em_within.tfs_perm.total = tf_total_perm;
    
        % save low temporal res permuted data
        em.tfs_cross_alpha_perm.allB1=allB1_perm;
        em.tfs_cross_alpha_perm.allB2=allB2_perm;
        
        save([iemDir '/' sprintf('sj%02d_cond%02d_IEM.mat',sn,iCond)],'em','em_within','minCnt','nElectrodes','-v7.3')

        clear em_within minCnt nElectrodes tf_evoked tf_total tf_evoked_perm tf_total_perm allB1 allB2 allB1_perm allB2_Perm allC1 permInd perm permedBins permedPosBin
        
        em.blocks = [];

    end
    
   % remember to clear any leftovers between subs/conds etc.
    
end















% 
% 
% 
% 
% %     for i=1:length(beh)
% %         posBin(i) = beh(i).stimLoc;   
% %     end
% 
% %==========================================================================
% %{
% Purpose: run spatial encoding model on single bandpass e.g. 8-12 Hz
% 
% Original Author:
% Joshua J. Foster
% joshua.james.foster@gmail.com
% University of Chicago
% August 12, 2015
% 
% Modified by Tom Bullock
% UCSB Attention Lab
% %}
% %==========================================================================
% 
% clear all
% close all
% 
% subs = [101:104 201:204 301:304 401:404 501:504 601:604 701:704 801:804 ...
%     901:904 1001:1004 1301:1304 1601:1604 1901:1904 2001:2004 2101:2104 ...
%     2201:2204 2301:2304 2401:2404];
% 
% nSubs = length(subs);
% 
% %% process for IEM (1) or HILBERT (2)
% processForIEM=1;
% 
% % setup directories
% root = '/home/bullock/WTF';
% dRoot = [root '/Data/'];
% 
% % % setup directories
% %
% % %root = pwd;
% % out = 'AnalysisScripts';
% % dRoot = [root(1:end-length(out)),'Data/'];
% 
% if processForIEM==1
%     %eRoot = [root(1:end-length(out)),'EEG/'];
%     eRoot = [root '/EEG/'];
% else
%     eRoot = [root(1:end-length(out)),'HILBERT/'];
% end
% 
% %hRoot = [root(1:end-length(out)),'HILBERT/'];
% hRoot = [root '/HILBERT/'];
% 
% %%%CLEAN UP%%%
% 
% % choose method for finding the min number of trials per location bin
% % (1=per-condition, 2=across conditions)
% findMinType = 2;
% 
% % load matrix of min trials per loc bin per condition
% minLocsData = load([dRoot 'Minimum_Location_Bin_Mat.mat']);
% 
% for bandpassLoop=1 % loop through different bandpass analyses e.g. alpha, theta
%     
%     if bandpassLoop==1
%         name = '_SpatialTF_ALPHA.mat'; % name of files to be saved for IEM
%         hilbName = '_Hilbert_ALPHA.mat';
%         thisFreq = {'Alpha'};
%         freqBandpass = [8 12];
%     else
%         name = '_SpatialTF_THETA.mat';
%         hilbName = '_Hilbert_THETA.mat';
%         thisFreq = {'Theta'};
%         freqBandpass = [4 7];
%     end
%     
%     % Loop through participants
%     matlabpool open 72
%     parfor s = 1:nSubs
%         
%         % had to move all the settings insdie the parfor to get it to work...
%         em = [];
%         EEG = [];
%         eegInfo = [];
%         
%         % parameters to set
%         em.nChans = 8; % # of channels
%         em.nBins = em.nChans; % # of stimulus bins
%         em.nIter = 10; % # of iterations %%% WAS SET TO 10!!!
%         em.nBlocks = 3; % # of blocks for cross-validation
%         em.frequencies = freqBandpass; % frequency bands to analyze
%         em.bands = thisFreq;
%         em.window = 4;
%         em.Fs = 256; % WAS 250
%         em.nElectrodes = 999; % get later
%         em.time = -.5*1000:1000/em.Fs:1.9961*1000; %   -500:4:2000; % time points of interest
%         
%         % for brevity in analysis
%         nChans = em.nChans;
%         nBins = em.nBins;
%         nIter = em.nIter;
%         nBlocks = em.nBlocks;
%         freqs = em.frequencies;
%         times = em.time;
%         nFreqs = size(em.frequencies,1);
%         nElectrodes = em.nElectrodes;
%         nSamps = length(em.time);
%         Fs = em.Fs;
%         
%         % Specify basis set
%         em.sinPower = 7;
%         em.x = linspace(0, 2*pi-2*pi/nBins, nBins);
%         em.cCenters = linspace(0, 2*pi-2*pi/nChans, nChans);
%         em.cCenters = rad2deg(em.cCenters);
%         pred = sin(0.5*em.x).^em.sinPower; % hypothetical channel responses
%         pred = wshift('1D',pred,5); % shift the initial basis function
%         basisSet = nan(nChans,nBins);
%         for c = 1:nChans;
%             basisSet(c,:) = wshift('1D',pred,-c); % generate circularly shifted basis functions
%         end
%         em.basisSet = basisSet; % save basis set to data structure
%         
%         sn = subs(s);
%         fprintf('Subject:\t%d\n',sn)
%         
%         % Grab data------------------------------------------------------------
%         
%         % Get position bin index from behavior file
%         fName = [dRoot, num2str(sn), '_MixModel_wBias.mat'];
%         
%         % load file
%         tmp = []; beh = [];
%         tmp = load(fName);
%         beh = tmp.beh;
%         
%         em.posBin = beh.trial.posBin'; % add to fm structure so it's saved
%         posBin = em.posBin;
%         
%         % Get EEG data
%         fName = [eRoot, num2str(sn), '_EEG.mat'];
%         
%         % load file
%         tmp = []; eeg = [];
%         tmp = load(fName);
%         eeg = tmp.eeg;
%         
%         % get n channels (to save later with TF files)
%         %%% nElects = size(eeg.chanLabels,2);
%         
%         eegs = eeg.data(:,:,:); % get scalp EEG (drop EOG electrodes)
%         artInd = eeg.arf.artIndCleaned.'; % grab artifact rejection index
%         tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time); nTimes = length(tois); % index time points for analysis.
%         
%         % %     %% TOM EDIT TO PROCESS WHOLE EPOCH (WTF DATA EPOCHED FROM -.5 to 2)
%         % %     tois = ones(1,size(eegs,3)); nTimes = length(tois);
%         
%         % Remove rejected trials
%         eegs = eegs(~artInd,:,:);
%         posBin = posBin(~artInd);
%         
%         em.nTrials = length(posBin); nTrials = em.nTrials; % # of good trials
%         
%         %----------------------------------------------------------------------
%         
%         % Preallocate Matrices
%         tf_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nChans); tf_total = tf_evoked;
%         C2_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nBins,nChans); C2_total = C2_evoked;
%         em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments
%         
%         % TOM ADDED (needed to convert to EEGLAB format to use new eeglab filter)
%         EEG.data = permute(eegs,[2,3,1]); % converts to chans x times x trials
%         EEG.srate = Fs;
%         EEG.trials = size(EEG.data,3);
%         EEG.nbchan = size(EEG.data,1);
%         EEG.pnts = size(EEG.data,2);
%         
%         % Loop through each frequency
%         for f = 1:nFreqs
%             tic % start timing frequency loop
%             fprintf('Frequency %d out of %d\n', f, nFreqs)
%             
%             %% get no. of electrodes
%             nElectrodes = size(eeg.data,2);
%             em.nElectrodes = nElectrodes;
%             disp(['nElecrodes changed to :' num2str(nElectrodes)])
%             
%             %% BUTTERWORTH FILTER
%             filterorder = 3;
%             type = 'bandpass';
%             [z1,p1] = butter(filterorder, [freqs(f,1), freqs(f,2)]./(EEG.srate/2),type);
%             %freqz(z1,p1,[],250)
%             data = double(EEG.data);
%             tempEEG = NaN(size(data,1),EEG.pnts,size(data,3));
%             for x = 1:size(data,1)
%                 for y = 1:size(data,3)
%                     dataFilt1 = filtfilt(z1,p1,data(x,:,y)); % was filtfilt
%                     tempEEG(x,:,y) = dataFilt1; % tymp = chans x times x trials
%                 end
%             end
%             
%             eegBand = [];
%             eegBand = tempEEG;
%             
%             %% apply hilbert to each channel and epoch in turn (this should be correct)
%             eegs = [];
%             for j=1:size(tempEEG,1) % chan loop
%                 for i=1:size(tempEEG,3) % trial loop
%                     eegs(i,j,:) = hilbert(squeeze(tempEEG(j,:,i)));
%                 end
%             end
%             
%             % eegs is trials x elects x times
%             fdata_evoked = eegs;
%             fdata_total = abs(eegs).^2;
%             
%             % Loop through each iteration
%             for iter = 1:nIter
%                 
%                 disp(['nIteration :' num2str(iter) ' of ' num2str(nIter)])
%                 %--------------------------------------------------------------------------
%                 % Assign trials to blocks (such that trials per position are
%                 % equated within blocks)
%                 %--------------------------------------------------------------------------
%                 
%                 % preallocate arrays
%                 blocks = nan(size(posBin));
%                 shuffBlocks = nan(size(posBin));
%                 
%                 % count number of trials within each position bin
%                 %clear binCnt
%                 binCnt = [];
%                 for bin = 1:nBins
%                     binCnt(bin) = sum(posBin == bin);
%                 end
%                 
%                 % choose method of determining the location bin with the min
%                 % number of trials (1= per condition, 2= across conditions)
%                 if findMinType==1
%                     minCnt = min(binCnt); % # of trials for position bin with fewest trials
%                 elseif findMinType==2
%                     %if minCnt is based across all four conditions.  We get
%                     %the values from "Minimum_Location_Bin_MAT.mat" created
%                     %earlier
%                     if ismember(sn,101:104)
%                         minCnt=minLocsData.minLocBinAllConds(1);
%                     elseif ismember(sn,201:204)
%                         minCnt=minLocsData.minLocBinAllConds(2);
%                     elseif ismember(sn,301:304)
%                         minCnt=minLocsData.minLocBinAllConds(3);
%                     elseif ismember(sn,401:404)
%                         minCnt=minLocsData.minLocBinAllConds(4);
%                     elseif ismember(sn,501:504)
%                         minCnt=minLocsData.minLocBinAllConds(5);
%                     elseif ismember(sn,601:604)
%                         minCnt=minLocsData.minLocBinAllConds(6);
%                     elseif ismember(sn,701:704)
%                         minCnt=minLocsData.minLocBinAllConds(7);
%                     elseif ismember(sn,801:804)
%                         minCnt=minLocsData.minLocBinAllConds(8);
%                     elseif ismember(sn,901:904)
%                         minCnt=minLocsData.minLocBinAllConds(9);
%                     elseif ismember(sn,1001:1004)
%                         minCnt=minLocsData.minLocBinAllConds(10);
%                     elseif ismember(sn,1301:1304)
%                         minCnt=minLocsData.minLocBinAllConds(11);
%                     elseif ismember(sn,1601:1604)
%                         minCnt=minLocsData.minLocBinAllConds(12);
%                     elseif ismember(sn,1901:1904)
%                         minCnt=minLocsData.minLocBinAllConds(13);
%                     elseif ismember(sn,2001:2004)
%                         minCnt=minLocsData.minLocBinAllConds(14);
%                     elseif ismember(sn,2101:2104)
%                         minCnt=minLocsData.minLocBinAllConds(15);
%                     elseif ismember(sn,2201:2204)
%                         minCnt=minLocsData.minLocBinAllConds(16);
%                     elseif ismember(sn,2301:2304)
%                         minCnt=minLocsData.minLocBinAllConds(17);
%                     elseif ismember(sn,2401:2404)
%                         minCnt=minLocsData.minLocBinAllConds(18);
%                     end
%                 end
%                 
%                 % reduce minCnt by 1 trial to ensure that ALL location bins
%                 % have some degree of trial randomization per IEM iteration
%                 % (without this, the loc bin with the minimum no. of trials
%                 % will not be randomized at all).
%                 minCnt = minCnt-1;
%                 
%                 nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
%                 
%                 % shuffle trials
%                 shuffInd = randperm(nTrials)'; % create shuffle index
%                 shuffBin = posBin(shuffInd); % shuffle trial order
%                 
%                 % take the 1st nPerBin x nBlocks trials for each position bin.
%                 for bin = 1:nBins;
%                     idx = find(shuffBin == bin); % get index for trials belonging to the current bin
%                     idx = idx(1:nPerBin*nBlocks); % drop excess trials
%                     x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
%                 end
%                 
%                 % unshuffle block assignment
%                 blocks(shuffInd) = shuffBlocks;
%                 
%                 % save block assignment
%                 em.blocks(:,iter) = blocks; % block assignment
%                 em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
%                 
%                 %-------------------------------------------------------------------------
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
%                         blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
%                         blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
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
%                         C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
%                         
%                         % Save weight matrix (Mary MacLean added)
%                         WeightsTotal(:,:,i) = W;
%                         
%                         % shift eegs to common center
%                         n2shift = ceil(size(C2,2)/2);
%                         for ii=1:size(C2,1)
%                             [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                             C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                         end
%                         
%                         tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                         
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
%                         C2_total(f,iter,t,i,:,:) = C2;
%                         
%                         % shift eegs to common center
%                         n2shift = ceil(size(C2,2)/2);
%                         for ii=1:size(C2,1)
%                             [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                             C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                         end
%                         
%                         tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                         %-----------------------------------------------------%
%                         
%                     end
%                     % Average weights across blocks add to matrix (Mary
%                     % MacLean, added)
%                     allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
%                 end
%             end
%             toc % stop timing the frequency loop
%         end
%         
%         %% TOM ADDED
%         % average over 1000 ITS + 3 BLOCKS to reduce size of saved file!
%         tf_evoked = squeeze(mean(mean(tf_evoked,2),4));
%         tf_total = squeeze(mean(mean(tf_total,2),4));
%         
%         % save data
%         if processForIEM==1
%             fName = [dRoot,num2str(sn),name];
%             em.C2.evoked = C2_evoked;
%             em.C2.total = C2_total;
%             em.tfs.evoked = tf_evoked;
%             em.tfs.total = tf_total;
%             em.tfs.totalW = allWeightsTotal;
%             em.nBlocks = nBlocks;
%             %save(fName,'em','-v7.3');
%             parsave(fName,em,minCnt,nElectrodes);
%         else
%             % save raw hilbert files
%             fName = [hRoot,num2str(sn),hilbName];
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
%     matlabpool close
%     
% end
% 
% sendEmailToMe('SPATIAL IEM SCRIPT FINISHED PROCESSING!!')
% 
% clear all
% close all
