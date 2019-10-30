%==========================================================================
%{
Beh_Merge_Blocks
Purpose: this merges the blocks for each individual subject number and
version.  Currently only does one sj/cond at a time. 
Author: Tom Bullock, UCSB Attention Lab
Date last modified: 03.04.18
%}
%==========================================================================

function A_behavior_master(subjectNumber,condition)

cdTmp = cd;

% this part merges the blocks
for iSub = subjectNumber 
    for iBlock = 1:15
        load(sprintf('sj%02d_bl%02d_cond%02d_fr.mat',iSub,iBlock,condition));  % load file
        % this part merges all blocks into one large mat
        if iBlock == 1;
            %allTrialData(1:size(dataMat,1),:) = dataMat;
            allTrialData(1:size(trialInfo,2)) = trialInfo;
        else
            %allTrialData((size(allTrialData,1)+1:size(allTrialData,1)+size(dataMat,1)),:) = dataMat;
            allTrialData(length(allTrialData)+1:length(allTrialData)+length(trialInfo))=trialInfo; 
        end  
    end
    % save stuff
    save(['/home/bullock/Foster_WM/TRIAL_DATA/MERGED_DATA' '/' sprintf('sj%d%02d_all_beh.mat',iSub,condition)],'allTrialData');
end

return