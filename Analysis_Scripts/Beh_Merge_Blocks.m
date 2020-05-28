%{
Beh_Merge_Blocks
Purpose: merge blocks
Author: Tom Bullock, UCSB Attention Lab
Date last modified: 10.30.19

%}

function Beh_Merge_Blocks(subjects)

%which subject(s)
subjects=[7];

% set dirs
rDir='D:\\WTF_EYE';
behDir=[rDir '\\' 'Beh_Data'];
behDirProcessed=[rDir '\\' 'Beh_Data_Processed'];

% rename raw files
Beh_Rename_Files(behDir)

for iSub=1:length(subjects)
    sjNum=subjects(iSub);
    
    clear allTrialData
    
    for iSession=1:2
        
        % find CB order
        load([behDir '\\' sprintf('sj%02d_se%02d_bl01_mc01_ec02.mat',sjNum,iSession)]);
        cbOrder=trialInfo(1).counterbalancingOrder;
        clear trialInfo
        
        for iCond=cbOrder
            
            clear allTrialData
            
%             if      iCond==1; mc=1;ec=2;
%             elseif  iCond==2; mc=1;ec=6;
%             elseif  iCond==3; mc=2;ec=2;
%             elseif  iCond==4; mc=2;ec=6;
%             end
            
            if      iCond==1; mc=1;ec=2;
            elseif  iCond==2; mc=2;ec=2;
            elseif  iCond==3; mc=1;ec=6;
            elseif  iCond==4; mc=2;ec=6;
            end
            
            clear demographicInfo
            
            for iBlock=1:4
                
                clear trialInfo 
                
                load([behDir '\\' sprintf('sj%02d_se%02d_bl%02d_mc%02d_ec%02d.mat',sjNum,iSession,iBlock,mc,ec) ])
                
                if iBlock == 1
                    allTrialData(1:size(trialInfo,2)) = trialInfo;
                else
                    allTrialData(length(allTrialData)+1:length(allTrialData)+length(trialInfo))=trialInfo;
                end
                
            end
            
            % compile into master structure
            masterStruct(iSession,iCond).allTrialData=allTrialData;
            
        end
        
    end
    
    save([behDirProcessed '\\' sprintf('sj%02d_allBeh.mat',sjNum) ],'masterStruct','demographicInfo','cbOrder')
    
end

return