%{
Spatial_IEM
Original Author: Joshua Foster (heavily edited by Tom Bullock, UCSB Attention Lab)
Date:11.15.19

Load preprocessed EEG and Beh data.
Merge across sessions and split into condition
Get time-frequency data
Either 1) Run IEM or 2) just get bandpassed data ...

%}

function IEM_Cross_Session_Gen(sn)

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
%bandpassedDir = [rDir '/' 'EEG_Bandpassed'];
iemDir = [rDir '/' 'IEM_Results_Cross_Session_Gen' ];
% select subjects
%subjects = 4;

addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts'))

% loop through subs
%for s=1:length(subjects)
%    sn=subjects(s);

% import IEM settings
em = IEM_Settings;
nChans = em.nChans;
nBins = em.nBins;
nIter = em.nIter;
nPerms = em.nPerms;
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
    
    clear eegs beh allConds allPosBin
    
    load([eegDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sn,iSession)])
    
    % reject artifact trials from eeg and beh
    eegs = eeg.data;
    eegs = eegs(:,:,~eeg.arf.artIndCleaned);
    beh = allBehStruct(~eeg.arf.artIndCleaned);
    
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
                
% %                 %==========================================================
% %                 %{
% %                 Run high temporal resolution IEM for "within" analysis only
% %                 %}
% %                 %==========================================================
% %                 
% %                 for t = 1:nSamps
% %                     
% %                     % grab data for timepoint t
% %                     toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
% %                     %find(toi==1)
% %                     de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
% %                     dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
% %                     
% %                     % Do forward model
% %                     for i=1:nBlocks % loop through blocks, holding each out as the test set
% %                         
% %                         trnl = labels(blockNum~=i); % training labels
% %                         tstl = labels(blockNum==i); % test labels
% %                         
% %                         %-----------------------------------------------------%
% %                         % Analysis on Evoked Power                            %
% %                         %-----------------------------------------------------%
% %                         B1 = de(blockNum~=i,:);    % training data
% %                         B2 = de(blockNum==i,:);    % test data
% %                         C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
% %                         W = C1\B1;          % estimate weight matrix
% %                         C2 = (W'\B2')';     % estimate channel responses
% %                         
% %                         C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
% %                         
% %                         % Save weight matrix (Mary MacLean added)
% %                         %WeightsTotal(:,:,i) = W;
% %                         
% %                         % shift eegs to common center
% %                         n2shift = ceil(size(C2,2)/2);
% %                         for ii=1:size(C2,1)
% %                             [~, shiftInd] = min(abs(posBins-tstl(ii)));
% %                             C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
% %                         end
% %                         
% %                         tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
% %                         
% %                         
% %                         %-----------------------------------------------------%
% %                         % Analysis on Total Power                             %
% %                         %-----------------------------------------------------%
% %                         B1 = dt(blockNum~=i,:);    % training data
% %                         B2 = dt(blockNum==i,:);    % test data
% %                         C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
% %                         W = C1\B1;          % estimate weight matrix
% %                         C2 = (W'\B2')';     % estimate channel responses
% %                         
% %                         C2_total(f,iter,t,i,:,:) = C2;
% %                         
% %                         % shift eegs to common center
% %                         n2shift = ceil(size(C2,2)/2);
% %                         for ii=1:size(C2,1)
% %                             [~, shiftInd] = min(abs(posBins-tstl(ii)));
% %                             C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
% %                         end
% %                         
% %                         tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
% %                         %-----------------------------------------------------%
% %                         
% %                     end
% %                     % Average weights across blocks add to matrix (Mary
% %                     % MacLean, added)
% %                     %%allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
% %                 end
                
                
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
                
                
%                 %==========================================================
%                 %{
%                 Run high temporal IEM permuted analyis
%                 %}
%                 %==========================================================
%                 
%                 % Loop through permutations
%                 for perm = 1:nPerms
%                     tic % start timing permutation loop
%                     %fprintf('Permutation %d out of %d/n',perm,nPerms);
%                     
%                     %-----------------------------------------------------------------------------
%                     % Permute trial assignment within each block
%                     %-----------------------------------------------------------------------------
%                     permedPosBin = nan(size(posBin)); % preallocate permuted position bins vector
%                     for b = 1:nBlocks % for each block..
%                         pInd = randperm(em.nTrialsPerBlock); % create a permutation index
%                         permedBins(pInd) = posBin(blocks == b); % grab block b data and permute according data according to the index
%                         permedPosBin(blocks == b) = permedBins; % put permuted data into permedPosBin
%                         permInd(f,iter,perm,b,:) = pInd; % save the permutation (permInd is saved at end of the script)
%                     end
%                     
%                     %-----------------------------------------------------------------------------
%                     
%                     % Average data for each position bin across blocks
%                     posBins = 1:nBins;
%                     blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
%                     blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
%                     labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
%                     blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
%                     c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
%                     bCnt = 1;
%                     for ii = 1:nBins
%                         for iii = 1:nBlocks
%                             blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(permedPosBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
%                             blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1));
%                             labels(bCnt) = ii;
%                             blockNum(bCnt) = iii;
%                             c(bCnt,:) = basisSet(ii,:);
%                             bCnt = bCnt+1;
%                         end
%                     end
%                     
%                     for t = 1:nSamps
%                         
%                         % grab data for timepoint t
%                         toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
%                         de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
%                         dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
%                         
%                         % Do forward model
%                         tmpeC2 = nan(nBlocks,nBins,nChans); tmptC2 = tmpeC2; % for unshifted channel responses
%                         tmpeCR = nan(nBlocks,nChans); tmptCR = nan(nBlocks,nChans); % for shifted channel respones
%                         
%                         for i=1:nBlocks % loop through blocks, holding each out as the test set
%                             
%                             trnl = labels(blockNum~=i); % training labels
%                             tstl = labels(blockNum==i); % test labels
%                             
%                             %-----------------------------------------------------%
%                             % Analysis on Evoked Power                            %
%                             %-----------------------------------------------------%
%                             B1 = de(blockNum~=i,:);    % training data
%                             B2 = de(blockNum==i,:);    % test data
%                             C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                             W = C1\B1;          % estimate weight matrix
%                             C2 = (W'\B2')';     % estimate channel responses
%                             
%                             % tmpeC2(i,:,:) = C2;
%                             
%                             % shift eegs to common center
%                             n2shift = ceil(size(C2,2)/2);
%                             for ii=1:size(C2,1)
%                                 [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                                 C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                             end
%                             
%                             tmpeCR(i,:) = mean(C2); % average shifted channel responses
%                             
%                             %-----------------------------------------------------%
%                             % Analysis on Total Power                             %
%                             %-----------------------------------------------------%
%                             B1 = dt(blockNum~=i,:);    % training data
%                             B2 = dt(blockNum==i,:);    % test data
%                             C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                             W = C1\B1;          % estimate weight matrix
%                             C2 = (W'\B2')';     % estimate channel responses
%                             
%                             % tmptC2(i,:,:) = C2;
%                             
%                             % shift eegs to common center
%                             n2shift = ceil(size(C2,2)/2);
%                             for ii=1:size(C2,1)
%                                 [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                                 C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                             end
%                             
%                             tmptCR(i,:) = mean(C2); % averaged shifted channel responses
%                             
%                             %-----------------------------------------------------%
%                             
%                         end
%                         % save data to indexed matrix
%                         % C2_evoked(f,iter,perm,t,:,:) = mean(tmpeC2);
%                         % C2_total(f,iter,perm,t,:,:) = mean(tmptC2);
%                         tf_evoked_perm(f,iter,perm,t,:) = mean(tmpeCR);
%                         tf_total_perm(f,iter,perm,t,:) = mean(tmptCR);
%                     end
%                     toc
%                 end
%                 
%                 %clear permInd permedBins
                
                
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
                    %fprintf('Permutation %d out of %d/n',perm,nPerms);
                    
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
        
%         % average high temporal resolution IEMs over block and iter to
%         % reduce file size for save
%         tf_evoked = squeeze(mean(mean(tf_evoked,2),4));
%         tf_total = squeeze(mean(mean(tf_total,2),4));
%         
%         tf_evoked_perm = squeeze(mean(mean(tf_evoked_perm,2),3));
%         tf_total_perm = squeeze(mean(mean(tf_total_perm,2),3));
%         
%         % organize high temporal res data
%         em_within.tfs.evoked = tf_evoked;
%         em_within.tfs.total = tf_total;
        
        % organize low temporal res data for cross training/testing
        em.tfs_cross_alpha.allB1=allB1;
        em.tfs_cross_alpha.allB2=allB2;
        
%         % save high temporal res permuted data
%         em_within.tfs_perm.evoked = tf_evoked_perm;
%         em_within.tfs_perm.total = tf_total_perm;
        
        % save low temporal res permuted data
        em.tfs_cross_alpha_perm.allB1=allB1_perm;
        em.tfs_cross_alpha_perm.allB2=allB2_perm;
        
        save([iemDir '/' sprintf('sj%02d_cond%02d_se%02d_IEM.mat',sn,iCond,iSession)],'em','minCnt','nElectrodes','-v7.3')
        
        clear em_within minCnt nElectrodes tf_evoked tf_total tf_evoked_perm tf_total_perm allB1 allB2 allB1_perm allB2_perm allC1 permInd perm permedBins permedPosBin B1 B2
        
        em.blocks = [];
        
    end
    
end

return