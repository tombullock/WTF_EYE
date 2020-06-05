%function SpatialEM_AllFs(sn)
%
% Run spatial encoding model across all frequency bands (see em.freqs).
%
% The model is not run at every samplev point. Instead, the sample rate for
% the analysis is specified by 'em.sampRate'.
%
% Cheked by Dave Sutterer 8.13.2015
%
% Joshua J. Foster
% joshua.james.foster@gmail.com
% University of Chicago
% August 12, 2015

clear all
close all

% %% classification settings
% thisSig = 'Feat'; %'Snr' or 'Ori'
% myFeat = 'Pwr'; % Pwr or Phase or PwrPhase or Coeff
% %nElecsMax=[56	22	49	36	57	58	39	50	37	57	60	60	23	39	21	53	52	16]; %number of electrodes it takes to get to max class acc collapsed across condition (after ranking)
% permuteLabels=0;

% sourceFolder = '/home/bullock/Foster_WM/Spectra';
% destFolder = '/home/bullock/Foster_WM/Classification';

%classType='loc'; %loc or ori
%classFreqs = 17:33;%33:57;

% setup directories
root = '/home/garrett/WTF_Bike/';
dRoot = [root 'Data/'];
outDir = [root 'SNR/'];
hRoot = [root 'HILBERT/'];

% load matrix of min trials per loc bin per condition
minLocsData = load([root '/Analysis_Scripts/HILBERT_Minimum_Location_Bin_Mat_AccTrials.mat']);
subjects = [1:8,10:35];


for iCon=1:2
    
    % Loop through participants
    for iSub = 1:length(subjects)
        sjNum = subjects(iSub);
        fprintf('Subject:\t%d\n',sjNum)
        
        % designate em structure at beginning, and clear it
        em = [];
        fdata_all = [];
        accData = [];
        accData_dt = [];
        accData_de = [];
        
        name = 'SNR_Hilbert_noAlpha.mat'; % name of files to be saved
        
        % parameters to set
        em.nChans = 8; % # of channels
        em.nBins = 8; % # of position bins
        em.nIter = 1; % # of iterations
        em.nBlocks = 3; % # of blocks for cross-validation
        fqs = 4:30;     % range of frequency bands
        for fbnd = 1:length(fqs)
            em.freqs(fbnd,:) = [fqs(fbnd) fqs(fbnd)+1]; % frequency bands to analyze
        end
        em.window = 4;
        em.Fs = 250;
        em.nFreqs = size(em.freqs,1);
        %em.nElectrodes = 20;
        em.sampRate =  4; %19.5312; % downsampled sample rate (in ms)
        %em.time = -500:4:2000; % specificy time window of interest
        
        %TOM ADDED
        em.time = -.5*1000:1000/em.Fs:1.9961*1000; %   -500:4:2000; % time points of interest
        
        em.dtime = -500:em.sampRate:2000; % downsampled time points
        em.stepSize = em.sampRate/(1000/em.Fs); % number of samples the classifer jumps with each shift
        em.nSamps = length(em.time); %  # of samples
        em.dSamps = length(em.dtime); % # of samples at downsampled rate
        
        em.dSamps = em.dSamps-1; %%%NOT SURE WHY, BUT THIS WORKS...
        
        
        % for brevity in analysis
        nChans = em.nChans;
        nBins = em.nBins;
        nIter = em.nIter;
        nBlocks = em.nBlocks;
        freqs = em.freqs;
        dtimes = em.dtime;
        nFreqs = size(em.freqs,1);
        %nElectrodes = em.nElectrodes;
        nSamps = em.nSamps;
        dSamps = em.dSamps;
        stepSize = em.stepSize;
        Fs = em.Fs;
        window = em.window;
        
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
        
        
        % Grab data------------------------------------------------------------
        
        % Get position bin index from behavior file
        fName = [dRoot sprintf( 'sj%02d_exerCon%02d_changeDect_MixModel_wBias_accTrials.mat',sjNum,iCon)];
        tmpA = load(fName);
        beh = tmpA.beh;
        tmpA.beh = [];
        em.posBin = beh.trial.posBin'; % add to class structure so it's saved
        posBin = em.posBin;
        
        % Get EEG data
        fName = [hRoot sprintf('sj%02d_exerCon%02d_changeDect_EEG_accTrials.mat',sjNum,iCon)];
        tmpA = load(fName);
        eeg = tmpA.eeg;
        tmpA.eeg = [];
        eegs = eeg.data(:,:,:); % get scalp EEG (drop EOG electrodes)
        artInd = eeg.arf.artIndCleaned.'; % grab artifact rejection index
        %tois = ismember(eeg.preTime:4:eeg.postTime,em.time); nTimes = length(tois); % index time points for analysis.
        tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time); nTimes = length(tois); % index time points for analysis.
        
        
        
        % GET nElectrodes (tom added here)
        nElectrodes = size(eegs,2);
        eegs = double(eegs); % also added this here
        
        % Remove rejected trials
        eegs = eegs(~artInd,:,:);
        posBin = posBin(~artInd);
        
        em.nTrials = length(posBin); nTrials = em.nTrials; % # of good trials
        
        %----------------------------------------------------------------------
        
        % Preallocate Matrices
        tf_evoked = nan(nFreqs,nIter,dSamps,nBlocks,nChans); tf_total = tf_evoked;
        em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments
        
        %--------------------------------------------------------------------------
        % Create block assignment for each iteration
        %--------------------------------------------------------------------------
        % trials are assigned to blocks so that # of trials per position are equated within blocks
        % this is done before the frequency loop so that the same blocks assignments are used for all freqs
        
        for iter = 1:nIter
            
            % preallocate arrays
            blocks = nan(size(posBin));
            shuffBlocks = nan(size(posBin));
            
            % count number of trials within each position bin
            binCnt=[];
            for bin = 1:nBins
                binCnt(bin) = sum(posBin == bin);
            end
            
            minCnt = min(binCnt); % # of trials for position bin with fewest trials
            nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
            
            % shuffle trials
            shuffInd = randperm(nTrials)'; % create shuffle index
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
            
        end
        %--------------------------------------------------------------------------
        
        % do Hilbert on Alpha only (8-12 Hz)
        % Filter Data
        fdata_evoked = nan(nTrials,nElectrodes,nTimes);
        %fdata_total = nan(nTrials,nElectrodes,nTimes);
        for c = 1:nElectrodes
            [b,a] = butter(3,[8 12]./Fs, 'stop');
            fdata_evoked(:,c,:) = hilbert(filtfilt(b,a,squeeze(eegs(:,c,:)))')';
        end
        
        % trim filtered data to remove times that are not of interest (after filtering to avoid edge artifacts)
        fdata_evoked = fdata_evoked(:,:,tois);
        
        %average over all trials to reduce the size of the matrix
        rawHilbert = squeeze(std(fdata_evoked,1));
        
        % save data
        fName = [outDir sprintf('sj%02d_exerCon%02d_%s',sjNum,iCon,name)];
        %em.hilb = fdata_all;
        %em.hilb.total = fdata_total;
        %save(fName,'em','-v7.3');
        save(fName,'rawHilbert');
        
        
    end
    
    
    
end