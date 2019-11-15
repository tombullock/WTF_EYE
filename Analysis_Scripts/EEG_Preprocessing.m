%{ 
Pre-process WTF EEG data
Purpose: From raw EEG > Epoched data
Author: Tom Bullock, UCSB Attention Lab
Date Created: 4.4.16
Date Last Modified: 3.3.18

Processed using EEGLAB v13.0.1
Must have already merged behavioral data (EXPLAIN)
%}

clear all
close all

% seed rng
for rngShuffle =1:100
    rng('shuffle')
end

% which subject number(s) to process?
subjectNumbers = [101:104 201:204 301:304 401:404 501:504 601:604 701:704 ...
    801:804 901:904 1001:1004 1301:1304 1601:1604 1901:1904 2001:2004 ...
    2101:2104 2201:2204 2301:2304 2401:2404];

%{
Subject Notes
302 - imported from bdf not edf
303 - exception added so don't try and import from edf/bdf again
604 - ditto
901 - eyelink edf file missing
904 - eyeling edf file corrupted
1902 - high number of dropped triggers in bdf
%}
    
% decine which processing steps to run
filterData=0;
synchEpochData=1;
artifactCorrection = 1; % artifact correction is a stage within this loop (keep switched on)
artifactRejectionThreshold=0;    disableThresholdRejection=0; % keep this off

% set directories for EEG (raw, filtered), Eye and Behavioral data
cdTmp = cd;
EEG_ftFolder = '/home/bullock/Foster_WM/EEG_ft'; %set EEG_ft directory
eyeSynchEEGfolder = '/home/bullock/Foster_WM/Eye_Synch_Folder'; %set eye-synching folder
eyeFolder = '/home/bullock/Foster_WM/EYE_raw'; % set raw eye data folder
unbrokenTrialsFolder = '/home/bullock/Foster_WM/Beh_Unbroken_Trials'; % set dir to store unbroken trials (for beh analysis)

%{
decide which prestimulus preprocessing to do: 
0=average baseline - use this for main spectral analyses, 
1= pre-stim baseline - use this for ERPs, EOG etc.
%}
whichPreprocessing=0;

if whichPreprocessing==0 % average baseline 
    EEG_ft_epFolder = '/home/bullock/Foster_WM/EEG_ft_ep_all';              %set EEG_ft_ep directory
    thisBaseline = [-500 2000];                                             %set baseline correction (whole epoch)
    doHilbert=1; %%%%?????   
elseif whichPreprocessing==1 % pre-stim baseline 
    EEG_ft_epFolder = '/home/bullock/Foster_WM/PRESTIM_BL_EEG_ft_ep_all';   %set EEG_ft_ep directory
    thisBaseline = [-500 0];                                                % set baseline correction (pre-stim)
    doHilbert=0;                                                            % turn hilbert off  
end

% create some blank variables (these are for synch consistency checks)
misMatchLog = [];

% open cluster pool (max 144)
%matlabpool open 72

%loop through subjects
for subjLoop=1:size(subjectNumbers,2)
    
    brkTrialVector=[]; EEG = [];                                            % create empty vars
    
    %% set subject numbers and define exceptions (e.g. split EEG files due to bathroom breaks)
    sjNum = subjectNumbers(subjLoop);
    disp(['PROCESSING SUBJECT: ' num2str(sjNum)])
    if sjNum==101||sjNum==303||sjNum==603||sjNum==604||sjNum==801||sjNum==1501||sjNum==1904||sjNum==1901||sjNum==2002||sjNum==2003||sjNum==2103
        mergedFile=1;
    else
        mergedFile=0;
    end
    
    %% import/filter data routine
    if filterData==1
        
        %cd into the raw data (edf) directorty
        cd '/home/bullock/Foster_WM/EEG_edfs'
        
        % merged or complete files (see "A_Merge_Broken_EEG_Files.m" script)
        if mergedFile==0
            if sjNum~=302 % exception - bdf won't export to edf
                d = dir(sprintf('sj%d.edf',sjNum));
            else
                d = dir(sprintf('sj%d.bdf',sjNum));
            end
            EEG = pop_biosig(d(1).name); 
        elseif mergedFile==1
            % gets filename for this subject & loads raw .set file
            d = dir(sprintf('sj%d.set',sjNum));
            EEG = pop_loadset(d(1).name);
        end
        
        % reference to mastoids
        EEG = pop_reref( EEG, [65 66],'keepref','on');
        
        % filter
        EEG = pop_eegfiltnew(EEG,0,80); % DO THIS FOR EOG
        
        % add channel locations
        EEG=pop_chanedit(EEG, 'lookup','/home/bullock/matlab_2011b/TOOLBOXES/eeglab13_0_1b/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp'); % cluster
        %EEG=pop_chanedit(EEG,'lookup','/Users/tombullock1/Documents/MATLAB/ML_TOOLBOXES/eeglab13_0_1b/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp'); % local
        
        % *subject exceptions*
        if sjNum==302 % eyetracker froze at start of bl08, had to restart display stuff
            EEG.event(3002:3010) = [];
            EEG.urevent(3002:3010) = [];
        end
        if sjNum==2101 % adds a 1 as the first event code (this was missing in the EEG data)
            EEG = pop_editeventvals(EEG,'insert',{1 [] [] []},'changefield',{1 'type' 1},'changefield',{1 'latency' 1});
        end
        
        % save filtered data
        EEG = pop_saveset(EEG,'filename',sprintf('%s_ft.set',d(1).name(1:end-4)),'filepath',EEG_ftFolder);
        
        %cd back to main directory
        cd '/home/bullock/Foster_WM'
        
    end
    
    %% epoch EEG data and synch with MAT file (trial info)
    if synchEpochData==1
        
        % gets filtered data filename for this subject
        d = dir([EEG_ftFolder '/' sprintf('sj%d_ft.set',sjNum)]);
        EEG = pop_loadset('filename',d(1).name,'filepath',EEG_ftFolder);
        
        % epoch files (first pass, epoch around fixation point because we
        % know that will be present in all trials, event broken ones)
        EEG = pop_epoch(EEG,{102},[-.5 4]);
        
        % diagnore missing triggers (only activate if you run into issues)
%         EEG= pop_epoch(EEG,{201 201 202 203 204 205 206 211 212 213 214 215 216},[0 .01])  
%         EEG = pop_rmbase(EEG,[]); % remove baseline from whole epoch
        
        %this segment prevents multiple target triggers being coded into an epoch
        %in the event of a broken trial (which messes up the later epoching around the targets)
        for i=1:length(EEG.epoch)
            
            %if merged file then events are strings, if not merged they are
            %numbers...
            if mergedFile==1
                %output vector indicating target stims in epoch
                tmpVec =ismember( EEG.epoch(i).eventtype,{'201','202','203','204','205','206','207','208'});
            elseif mergedFile==0
                %output vector indicating target stims in epoch
                tmpVec = ismember(cell2mat(EEG.epoch(i).eventtype),[201,202,203,204,205,206,207,208]);
            end
            
            %if tmpVec indicates more than one target code in epoch
            if sum(tmpVec)>1
                thisCode = EEG.epoch(i).eventtype(tmpVec);
                EEG.epoch(i).eventtype = thisCode(1); %broken trial to scrap all triggers except the first single target one
            end
        end
        i = []; tmpVec = []; thisCode = [];
        
        % subject exception (last block mat data lost for sj 603)
        if sjNum==603
            EEG = pop_select(EEG,'trial',1:993);
        end
        
        % load behavioral data from matlab (.mat file, must be merged)
        trialMat = [];
        trialMat = load(['/home/bullock/Foster_WM/TRIAL_DATA/MERGED_DATA' '/' sprintf('sj%d_all_beh.mat',sjNum)]);
        
        % subject exception
        if sjNum==1601
            trialMat.allTrialData(1291:1292) = [];
        end
        
        % check for EEG/MAT consistency, break script if mismatched (note: 
        % can only activate "break" if not running parfor)
        if length(trialMat.allTrialData) == length(EEG.epoch)
            disp('EEG AND MAT FILES MATCHING LENGTH!!!')
        else
            disp('MISMATCH BETWEEN EEG AND MAT FILES!!!')
            %break
        end
        
        % add trialinfo structure to the EEG
        EEG.trialInfo = trialMat.allTrialData;
        trialMat.allTrialData = [];
        
        % idenify broken trials (create vector)
        brkCounter=0;
        for i=1:length(EEG.trialInfo)
            if EEG.trialInfo(i).brokenTrial == 0
                brkCounter=brkCounter+1;
                brkTrialVector(brkCounter) = i;
                EEG.newTrialInfo(brkCounter) = EEG.trialInfo(i);
            end
        end
        
        % clears old EEG.trialInfo to avoid confusion at later stages!
        EEG.trialInfo = [];
        
        % selects unbroken trials only
        EEG = pop_select(EEG,'trial',brkTrialVector);
        
        % subject exception 
        if sjNum==1601
            EEG = pop_select(EEG,'trial',1:957);
            EEG.newTrialInfo(958) = [];
        end
        
        % check EEG.epoch and EEG.newTrialInfo contain same number of
        % trials (another consistency check)
        if length(EEG.epoch)==length(EEG.newTrialInfo)
            disp(['EEG NOW CONTAINS ' num2str(length(EEG.epoch)) ' UNBROKEN TRIALS'])
        else
            disp('EEG.epoch or EEG.newTrialInfo DOES NOT CONTAIN CORRECT NO. TRIALS!  ABORT!')
            %break
        end
        
        % re-eoch data to [-.5 to 2.5] around target onset
        EEG = pop_epoch(EEG,{201 202 203 204 205 206 207 208},[-.5 2.5]); 
        
        % subject exception
        if sjNum==303
            EEG.newTrialInfo(265) =[];
        end
        
        % use this to diagnose mismatching data (occasional dropped epoch
        % code will result in EEG/BEH files not synching properly)
        findMismatchMat = [];
        if ismember(sjNum,[101,303,603,604,801,1501,1904,1901,2002,2003,2103])   % if codes are strings          
            for i=1:size(EEG.epoch,2)
                findMismatchMat(i,:) = [str2num(EEG.epoch(i).eventtype{1}), EEG.newTrialInfo(i).locTrigger];
            end
            % if 255 is first event, use second event code for comparison
            if findMismatchMat(i,1)==255|| findMismatchMat(i,1)==223
                findMismatchMat(i,:) = [EEG.epoch(i).eventtype{2}, EEG.newTrialInfo(i).locTrigger];
            end       
        else % if event codes are numerical     
            for i=1:size(EEG.epoch,2)
                findMismatchMat(i,:) = [EEG.epoch(i).eventtype{1}, EEG.newTrialInfo(i).locTrigger];  
                % if 255 is first event, use second event code for comparison
                if findMismatchMat(i,1)==255|| findMismatchMat(i,1)==223
                    findMismatchMat(i,:) = [EEG.epoch(i).eventtype{2}, EEG.newTrialInfo(i).locTrigger];
                end    
            end   
        end
        
        % gen. a third column to check epoch and trial values
        for i=1:size(EEG.epoch,2)
            findMismatchMat(i,3) = findMismatchMat(i,1) - findMismatchMat(i,2);
        end
        if sum(findMismatchMat(:,3))==0
            thisMatch=1;
            disp('EPOCH CODES CONSISTENT WITH MAT CODES!!!');
        else
            thisMatch=0;
            disp(['EPOCH CODES INCONSISTENT WITH MAT CODES FOR SJ ' num2str(sjNum)  ' -ABORT!!!'])
            %break
        end
        
        % create a mismatch log
        misMatchLog(subjLoop,:) = [sjNum, thisMatch];
        
        % remove baseline (from whole epoch/baseline - see earlier setting)
        EEG = pop_rmbase(EEG,thisBaseline); % remove baseline from whole epoch
        
        % do CRLS artifact correction (regresses out eye-closure shifts that 
        % would otherwise cause issues with later threshold based art. rej )
        if artifactCorrection==1
            EEG = pop_crls_regression( EEG, [67:72], 1, 0.9999, 0.01,[]); 
        end
        
        %save the synched, epoched data
        EEG = pop_saveset(EEG,'filename',sprintf('%s_ep.set',d(1).name(1:end-4)),'filepath',EEG_ft_epFolder);
        
        %save trial data for behavioral analysis only (this excludes
        %broken trials)
        trialInfoUnbroken = EEG.newTrialInfo;
        parsave([unbrokenTrialsFolder '/' sprintf('sj%d_newBeh.mat',sjNum)],trialInfoUnbroken)
        trialInfoUnbroken = [];
        
    end
    
end

% save mismatch checker
save('misMatchLog.mat','misMatchLog')

% close matlab pool
matlabpool close

% run a function to find the fewest trials across all four conditions for
% each subject (only run if I also run doHilbert)
if doHilbert==1
    findFewestTrials % saves a mat called "Minimum_Location_Bin_Mat.mat" in the main dir
end

clear all
close all