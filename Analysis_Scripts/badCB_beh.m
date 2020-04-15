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
mc=1;
ec=2;

% load data
load(sprintf([behDir '\\sj%02d_se%02d_bl%02d_mc%02d_ec%02d'],sjNum,iSession,iBlock,mc,ec))

for iTrial=1:length(trialInfo)
    badCB = trialInfo(iTrial).counterbalancingOrder;
    
    if isequal(badCB,[1,4,3])    %change block numbers if not what it originally stated
        trialInfo(iTrial).counterbalancingOrder = [2,1,4,3];
    end
    
end

save([behDir '\\' sprintf('sj%02d_se%02d_bl%02d_mc%02d_ec%02d',sjNum,iSession,iBlock,mc,ec) ],'demographicInfo','trialInfo','thisTimeDate')




