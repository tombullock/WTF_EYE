%{
EEG_Preprocessing2
Author: Tom Bullock, UCSB Attention Lab
Date: 11.14.19

Remove baseline (depending on analysis)
Remove bad channels (try cleanrawdata)
Do artifact correction (see how this copes with EOG in eyes-moved conds)
Do threshold based artifact rejection
Reshape data for further processing

%}

function EEG_Preprocessing2(subjects)

if nargin==0
    subjects = input('Input Subject(s): ')
end

% add paths
cd('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts');
cd('/home/waldrop/Desktop/WTF_EYE/Dependancies/eeglab14_1_2b');
eeglab
%clear
close all
cd('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts');

% choose subjects
%subjects = [4];%[1:6];

% choose sessions
theseSessions=1:2;

% show rejected trials? (0=no,1=yes)
showRejTrials=0;

% set dirs
rDir = '/home/waldrop/Desktop/WTF_EYE';
eegPrepro1Dir = [rDir '/' 'EEG_Prepro1'];

% choose analysis type
whichPreprocessing=0;
if whichPreprocessing==0 % average baseline (use for main spectral\\IEM analyses)
    destDir = [rDir '/' 'EEG_Prepro2_Avg_Baseline'];
    thisBaseline = [-500 2000]; 
    interpolateBadChannels=0;
elseif whichPreprocessing==1 % pre-stim baseline (use for ERPs, EOG etc.)
    destDir = [rDir '/' 'EEG_Prepro2_Prestim_Baseline'];
    thisBaseline = [-500 0]; 
    interpolateBadChannels=1;
end

% loop through subs
for iSub=1:length(subjects)
    sjNum=subjects(iSub);
    
    % loop through sessions
    for iSession=theseSessions
        
        % load data
        load([eegPrepro1Dir '/' sprintf('sj%02d_se%02d_wtfEye_prepro1.mat',sjNum,iSession) ]); 
        
        % remove baseline
        EEG = pop_rmbase(EEG,thisBaseline); 
        
        % remove eye activity
        EEG = pop_crls_regression( EEG, [67:72], 1, 0.9999, 0.01,[]); 
   
        % remove reference channels (1&2) and EOG (3:8)
        EEG = pop_select(EEG,'nochannel',{'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8'});

        % list bad channels via visual inspection (both sessions)!
        EEG.original_chanlocs = EEG.chanlocs;
        eeg.chanlocsOriginal = EEG.chanlocs;
        badChannels = badChannelInfo(sjNum);
        eeg.badChannels=badChannels;
        
        % remove bad channels from EEG
        EEG = pop_select(EEG,'nochannel',badChannels);

        % threshold based artifact rejection (just store an index of removed trials)
        if showRejTrials==0
            EEG  = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -150 ,150, 0,2, 1, 0); % only [0-2secs]
            eeg.arf.artIndCleaned = EEG.reject.rejthresh;
        else
            EEG  = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -150 ,150, 0,2, 1, 0); % only [0-2secs]
            eeg.arf.artIndCleaned = EEG.reject.rejthresh;
            plotRejThr=trial2eegplot(EEG.reject.rejthresh,EEG.reject.rejthreshE,EEG.pnts,EEG.reject.rejthreshcol);
            rejE=plotRejThr;
            %Draw the data.
            m=get(0,'MonitorPositions');
            eegplot(EEG.data,...
                'eloc_file',EEG.chanlocs, ...
                'srate',EEG.srate,...
                'events',EEG.event,...
                'winrej',rejE,...
                'spacing',100,...
                'winlength',5,...
                'title',sprintf('sj%02d_se%02d Plot',sjNum,iSession),...
                'position',m(1,:));
            
            % option to move forwards and complete processing OR exit from
            % script, enter bad channels and re-run
            disp('Are there any bad channels?  "Y" = Yes, "N" = no')
            fig=gcf;
            pressloop=1;
            while pressloop
                was_a_key = waitforbuttonpress;
                if was_a_key && strcmp(get(fig, 'CurrentKey'), 'y')
                    pressloop=0;
                    disp('EXITING SCRIPT - INPUT BAD CHANNELS AND RERUN!')
                    return
                elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'n')
                    pressloop=0;
                end
            end
            close;
        end
        
        % interpolate bad channels (can always reject later for IEM)
        EEG = pop_interp(EEG,EEG.original_chanlocs, 'spherical');
        
        % log percent trials rejected
        eeg.pcTrialsRej = round([sum(EEG.reject.rejthresh)/length(EEG.reject.rejthresh)*100]);
        %pcTrialsRej = eeg.pcTrialRej;

        % NO! permute EEG data to trials x chans x samples
        %eeg.data = permute(EEG.data,[3,1,2]);
        eeg.data = EEG.data;
        
        % extract channel labels
        for iChans = 1:length(EEG.chanlocs)
            eeg.chanLabels{iChans} = EEG.chanlocs(iChans).labels;
        end
        
        % add some extra info for ease of data processing
        eeg.preTime = EEG.xmin*1000;
        eeg.postTime = EEG.xmax*1000;
        eeg.sampRate = EEG.srate;
        eeg.pnts = EEG.pnts;
        eeg.times = EEG.times;
        
        % save the behavior data
        disp('Saving data...');
        
        % save data
        save([destDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum,iSession)],'eeg','allBehStruct','demographicInfo','-v7.3')
        
        clear eeg EEG badChannels allBehStruct
        
    end

end

return