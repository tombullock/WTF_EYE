%{
EEG_Preprocessing1
Author: Tom Bullock
Date: 11.12.19

Import raw bdfs
Re-ref to mastoids
Add channel locations
Downsample to 256 Hz
Filter
Save as mat file

% REMEMBER TO BASELINE CORRECT IN SUBSEQUENT ANALYSIS SCRIPTS! %

%}

% add paths
cd('/Users/tombullock/Documents/Psychology/WTF_EYE/Analysis_Scripts');
cd('/Users/tombullock/Documents/MATLAB/eeglab14_1_2b');
eeglab
clear
close all
cd('/Users/tombullock/Documents/Psychology/WTF_EYE/Analysis_Scripts');

% choose subjects
subjects = [5];%[1:6];

% were data pre-merged?
preMergedDataSubjects = [1];

% set dirs
rDir = '/Users/tombullock/Documents/Psychology/WTF_EYE';
eegRawDir = [rDir '/' 'EEG_Raw'];
eegPrepro1Dir = [rDir '/' 'EEG_Prepro1'];
behDir = [rDir '/' 'Beh_Data_Processed'];

% % choose analysis type
% whichPreprocessing=0;
% if whichPreprocessing==0 % average baseline (use for main spectral analyses)
%     EEG_ft_epFolder = '/home/bullock/Foster_WM/EEG_ft_ep_all';              %set EEG_ft_ep directory
%     thisBaseline = [-500 2000];                                             %set baseline correction (whole epoch)
% elseif whichPreprocessing==1 % pre-stim baseline (use for ERPs, EOG etc.)
%     EEG_ft_epFolder = '/home/bullock/Foster_WM/PRESTIM_BL_EEG_ft_ep_all';   %set EEG_ft_ep directory
%     thisBaseline = [-500 0];                                                % set baseline correction (pre-stim)                                                           
% end

% loop through subs
for iSub=1:length(subjects)
    sjNum=subjects(iSub);
    
    % load behavioral data
    load([behDir '/' sprintf('sj%02d_allBeh.mat',sjNum)])
    
    % loop through sessions
    for iSession=1:2
        
        % was data file pre-merged?
        if ismember(sjNum,preMergedDataSubjects)
            mergedFile=1;
        else
            mergedFile=0;
        end
        
        % import/load data
        if mergedFile
            load([eegRawDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum,iSession) ]); 
        else
            EEG=pop_biosig([eegRawDir '/' sprintf('sj%02d_se%02d_wtfEye.bdf',sjNum,iSession)]);
        end
        
        % downsample data
        if sjNum==1 && iSession==1
            disp('Data already resampled')
        else
            EEG = pop_resample(EEG,256);
        end
        
        % reference to mastoids
        EEG = pop_reref( EEG, [65 66],'keepref','on');
        
        % filter
        EEG = pop_eegfiltnew(EEG,0,80); % DO THIS FOR EOG
        
        % add channel locations
        EEG=pop_chanedit(EEG, 'lookup','/Users/tombullock/Documents/MATLAB/eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'); 
        
        % if merged file, convert EEG event codes from strings to nums
        if mergedFile
            for i=1:length(EEG.event)
                EEG.event(i).type = str2double(EEG.event(i).type);
            end
        end
        
        % epoch files around fixation point (include ALL trials)
        EEG = pop_epoch(EEG,{102},[-.5 4]);
        
        % prevent multiple targ triggers being coded into an epoch if trial
        % broken (which screws up later epoching around targets)
        for i=1:length(EEG.epoch)
            tmpVec = ismember(cell2mat(EEG.epoch(i).eventtype),[201,202,203,204,205,206,207,208]);
            %if tmpVec indicates more than one target code in epoch
            if sum(tmpVec)>1
                thisCode = EEG.epoch(i).eventtype(tmpVec);
                EEG.epoch(i).eventtype = thisCode(1); %broken trial to scrap all triggers except the first single target one
            end
        end
      
        % organize trial data in cronological order
        for iCond=1:4
            if iCond==1
                allBehStruct(1:size(masterStruct(iSession,cbOrder(iCond)).allTrialData,2)) = masterStruct(iSession,cbOrder(iCond)).allTrialData;
            else
                allBehStruct(length(allBehStruct)+1:length(allBehStruct)+length(masterStruct(iSession,cbOrder(iCond)).allTrialData))=masterStruct(iSession,cbOrder(iCond)).allTrialData;     
            end
        end
        
        % create broken trial vector
        cnt=0;
        for i=1:length(allBehStruct)          
           if allBehStruct(i).brokenTrial
               cnt=cnt+1;
              brokenTrialVec(cnt)=i;
           end    
        end
        
        % remove broken trials from EEG and trialData
        EEG = pop_select(EEG,'notrial',brokenTrialVec);
        allBehStruct(brokenTrialVec) = [];
        
        % re-epoch data to [-.5 to 2.5] around target onset
        EEG = pop_epoch(EEG,{201 202 203 204 205 206 207 208},[-.5 2.5]);
        
        % check that EEG and trial data match by checking event codes
        clear consistencyCheck
        for i=1:length(EEG.epoch)
           consistencyCheck(i,:) = [EEG.epoch(i).eventtype{1},allBehStruct(i).locTrigger,EEG.epoch(i).eventtype{1}-allBehStruct(i).locTrigger];          
        end
        
        if sum(consistencyCheck(:,3))~=0
            disp(sprintf('EEG/TRIAL DATA MISMATCH! INVESTIGATE sj%02d_se%02d!!',sjNum,iSession))
            return
        else
            disp('EEG/TRIAL MATCH!! YAY!!')
        end
        
        % save the synced, epoched data
        save([eegPrepro1Dir '/' sprintf('sj%02d_se%02d_wtfEye_prepro1.mat',sjNum,iSession)],'EEG','allBehStruct','demographicInfo')
        
        clear brokenTrialVec  cnt consistencyCheck EEG thisCode tmpVec allBehStruct 
        
    end
end
        
        
        
        
            
            
        
        
        
        

        
        
        

        
      
%         if iBlock == 1
%                     allTrialData(1:size(trialInfo,2)) = trialInfo;
%                 else
%                     allTrialData(length(allTrialData)+1:length(allTrialData)+length(trialInfo))=trialInfo;
%                 end
%         
%         
%         
%         
%         % sync eeg and beh data
%         allBehStruct = [
%             masterStruct(iSession,cbOrder(1));...
%             masterStruct(iSession,cbOrder(2));...
%             masterStruct(iSession,cbOrder(3));...
%             masterStruct(iSession,cbOrder(4))];
%             
        
        
%         % save data
%         %save([
%     end
%     
%     
% end





% %% set subject numbers and define exceptions (e.g. split EEG files due to bathroom breaks)
%     sjNum = subjectNumbers(subjLoop);
%     disp(['PROCESSING SUBJECT: ' num2str(sjNum)])
%     if sjNum==101||sjNum==303||sjNum==603||sjNum==604||sjNum==801||sjNum==1501||sjNum==1904||sjNum==1901||sjNum==2002||sjNum==2003||sjNum==2103
%         mergedFile=1;
%     else
%         mergedFile=0;
%     end
%     
%     %% import/filter data routine
%     if filterData==1
%         
%         %cd into the raw data (edf) directorty
%         cd '/home/bullock/Foster_WM/EEG_edfs'
%         
%         % merged or complete files (see "A_Merge_Broken_EEG_Files.m" script)
%         if mergedFile==0
%             if sjNum~=302 % exception - bdf won't export to edf
%                 d = dir(sprintf('sj%d.edf',sjNum));
%             else
%                 d = dir(sprintf('sj%d.bdf',sjNum));
%             end
%             EEG = pop_biosig(d(1).name); 
%         elseif mergedFile==1
%             % gets filename for this subject & loads raw .set file
%             d = dir(sprintf('sj%d.set',sjNum));
%             EEG = pop_loadset(d(1).name);
%         end
%         
%         % reference to mastoids
%         EEG = pop_reref( EEG, [65 66],'keepref','on');
%         
%         % filter
%         EEG = pop_eegfiltnew(EEG,0,80); % DO THIS FOR EOG
%         
%         % add channel locations
%         EEG=pop_chanedit(EEG, 'lookup','/home/bullock/matlab_2011b/TOOLBOXES/eeglab13_0_1b/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp'); % cluster
%         %EEG=pop_chanedit(EEG,'lookup','/Users/tombullock1/Documents/MATLAB/ML_TOOLBOXES/eeglab13_0_1b/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp'); % local
%         
%         % *subject exceptions*
%         if sjNum==302 % eyetracker froze at start of bl08, had to restart display stuff
%             EEG.event(3002:3010) = [];
%             EEG.urevent(3002:3010) = [];
%         end
%         if sjNum==2101 % adds a 1 as the first event code (this was missing in the EEG data)
%             EEG = pop_editeventvals(EEG,'insert',{1 [] [] []},'changefield',{1 'type' 1},'changefield',{1 'latency' 1});
%         end
%         
%         % save filtered data
%         EEG = pop_saveset(EEG,'filename',sprintf('%s_ft.set',d(1).name(1:end-4)),'filepath',EEG_ftFolder);
%         
%         %cd back to main directory
%         cd '/home/bullock/Foster_WM'
%         
%     end
%     


