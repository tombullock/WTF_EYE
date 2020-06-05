%{
IEM_Cross_TT
Author: Tom Bullock (using code borrowed from Joshua Foster)
Date: 11.18.19

Train and test within and across all conditions

To Do: build in the permuted data analysis
%}

function IEM_Cross_Session_TT(sjNum)

addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts'))

%clear
%close all

% set source dir
sourceDir = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_Cross_Session_Gen';
destDir  = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Cross_Session_Gen';

% select subjects
%%subjects = [10:14,16:20,22:27,31];

% loop through all training\testing combinations
%for iSub=1:length(subjects)

%sjNum=subjects(iSub);

em=[]; tf_total=[]; tf_total_cross=[];

% loop through 88 train/test generalizations (train se01 test se02 / train
% se02 test se01)
for k=1:8
    
    disp(['Subject Number: ' num2str(sjNum) ' Combination Number: ' num2str(k)])
        
    if      k==1; trainCond=1; testCond=1; trainSession=1; testSession=2;
    elseif  k==2; trainCond=2; testCond=2; trainSession=1; testSession=2;
    elseif  k==3; trainCond=3; testCond=3; trainSession=1; testSession=2;
    elseif  k==4; trainCond=4; testCond=4; trainSession=1; testSession=2;
    elseif  k==5; trainCond=1; testCond=1; trainSession=2; testSession=1;
    elseif  k==6; trainCond=2; testCond=2; trainSession=2; testSession=1;
    elseif  k==7; trainCond=3; testCond=3; trainSession=2; testSession=1;
    elseif  k==8; trainCond=4; testCond=4; trainSession=2; testSession=1;
    end
    
    trainSet = load([sourceDir '/' sprintf('sj%02d_cond%02d_se%02d_IEM.mat',sjNum,trainCond,trainSession)]);
    testSet = load([sourceDir '/' sprintf('sj%02d_cond%02d_se%02d_IEM.mat',sjNum,testCond,testSession)]);
    
    
    % set up some vars
    allB1=trainSet.em.tfs_cross_alpha.allB1;
    allB2=testSet.em.tfs_cross_alpha.allB2;
    
    nChans = trainSet.em.nChans;
    nBins = trainSet.em.nBins;
    nIter = trainSet.em.nIter;
    nBlocks = trainSet.em.nBlocks;
    basisSet = trainSet.em.basisSet;
    
    % iterate through x times
    for iter=1:nIter
        
        % Average data for each position bin across blocks
        posBins = 1:nBins;
        %blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
        %blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
        labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
        blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
        c = nan(nBins*nBlocks,nChans);           %%TOM                % predicted channel responses for averaged data
        
        bCnt = 1;
        for ii = 1:nBins
            for iii = 1:nBlocks
                %blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                %blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
                labels(bCnt) = ii;
                blockNum(bCnt) = iii;
                c(bCnt,:) = basisSet(ii,:);
                bCnt = bCnt+1;
            end
        end
        
        for trLoop=1:40
            
            for teLoop=1:40
                
                for i=1:nBlocks
                    
                    trnl = [1;1;2;2;3;3;4;4;5;5;6;6;7;7;8;8;];
                    tstl = [1;2;3;4;5;6;7;8;];
                    
                    % IEM
                    B1 = squeeze(allB1(iter,trLoop,teLoop,i,:,:));
                    B2 = squeeze(allB2(iter,trLoop,teLoop,i,:,:));
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
                    
                    tf_total(iter,trLoop,teLoop,i,:) = mean(C2,1);
                    
                end
            end
        end
    end
    
    % average over iterations and blocks to reduce file
    tf_total_cross(k,:,:,:) = squeeze(mean(mean(tf_total,1),4));
    
    clear tf_total
    
    
    
    %% RUN PERMUTED CROSS TT %%
    % set up some vars
    allB1=trainSet.em.tfs_cross_alpha_perm.allB1;
    allB2=testSet.em.tfs_cross_alpha_perm.allB2;
    
    nChans = trainSet.em.nChans;
    nBins = trainSet.em.nBins;
    nIter = trainSet.em.nIter;
    nPerm = trainSet.em.nPerms;
    nBlocks = trainSet.em.nBlocks;
    basisSet = trainSet.em.basisSet;
    
    % go through perms
    for perm=1:nPerm
        % iterate through x times
        for iter=1:nIter
            
            % Average data for each position bin across blocks
            posBins = 1:nBins;
            %blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
            %blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
            labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
            blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
            c = nan(nBins*nBlocks,nChans);           %%TOM                % predicted channel responses for averaged data
            
            bCnt = 1;
            for ii = 1:nBins
                for iii = 1:nBlocks
                    %blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    %blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
                    labels(bCnt) = ii;
                    blockNum(bCnt) = iii;
                    c(bCnt,:) = basisSet(ii,:);
                    bCnt = bCnt+1;
                end
            end
            
            for trLoop=1:40
                
                for teLoop=1:40
                    
                    for i=1:nBlocks
                        
                        trnl = [1;1;2;2;3;3;4;4;5;5;6;6;7;7;8;8;];
                        tstl = [1;2;3;4;5;6;7;8;];
                        
                        % IEM
                        B1 = squeeze(allB1(perm,iter,trLoop,teLoop,i,:,:));
                        B2 = squeeze(allB2(perm,iter,trLoop,teLoop,i,:,:));
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
                        
                        tf_total(perm,iter,trLoop,teLoop,i,:) = mean(C2,1);
                        
                    end
                end
            end
        end
    end
    
    % average over iterations and blocks to reduce file
    tf_total_cross_perm(k,:,:,:) = squeeze(mean(mean(mean(tf_total,1),2),5));
    
    clear tf_total
    
    
end


% calculate and save slopes (real then perm)
thisX = 0:45:180; % WTF use real angular values
for permReal=1:2
    
    if permReal==1; theseData = tf_total_cross;
    elseif permReal==2; theseData = tf_total_cross_perm;
    end
    
    for iCross=1:8
        for trSamps=1:40
            for teSamps=1:40
                %dat=rDat.total(iCross,trSamps,teSamps,:);
                dat = theseData(iCross,trSamps,teSamps,:);
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                
                % for AUC?
                %fit = trapz(d);
                
                % typical slope calc.
                fit = polyfit(x,d,1);
                
                if      permReal==1; em.rSl.total(iCross,trSamps,teSamps)=fit(1);
                elseif  permReal==2; em.pSl.total(iCross,trSamps,teSamps)=fit(1);
                end
            end
        end
    end
end

% save em
em.tfs.total_real = tf_total_cross;
em.tfs.total_perm = tf_total_cross_perm;

save([destDir '/' sprintf('sj%02d_IEM_Alpha_Cross_Session_TT.mat',sjNum)],'em');

return