clear
close all

% set dirs
rDir='D:\\WTF_EYE';
behDir=[rDir '\\' 'Beh_Data'];

%which subjects/session
sjNum=14;
iSession=1;

%which block
iBlock=4;

%which condition
mc=2;
ec=2;

% load data
load(sprintf([behDir '\\sj%02d_se%02d_bl%02d_mc%02d_ec%02d'],sjNum,iSession,iBlock,mc,ec))

%change block number in trialInfo struct
for iTrial=1:length(trialInfo)
    
    badBlock = trialInfo(iTrial).thisBlock;
    if badBlock==1    %change block if this is not what it originally stated
        trialInfo(iTrial).thisBlock = iBlock;
    end
end

save([behDir '\\' sprintf('sj%02d_se%02d_bl%02d_mc%02d_ec%02d',sjNum,iSession,iBlock,mc,ec)],'demographicInfo','trialInfo','thisTimeDate')
