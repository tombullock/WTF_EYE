%{
Plot_IEM_Cross_TT
Author: Tom Bullock, UCSB Attention Lab
Date: 11.18.19

To Do: build in perm data analysis/comparison

%}

clear 
close all

% which directories?
rDir = '/home/waldrop/Desktop/WTF_EYE';
sourceDir = [rDir '/' 'IEM_Results_TT_Cross'];
destDir = [rDir '/' 'Data_Compiled'];

% which subjects?
subs = [1:7,9:14,16:20,22:27,31];
%

% compile slope data
for iSub=1:length(subs)
    
    sjNum=subs(iSub);
    
    % load data
    load([sourceDir '/' sprintf('sj%02d_IEM_Alpha_CrossTT.mat',sjNum)])
    allTF_real_total(iSub,:,:,:)=em.rSl.total; %subs x cross x tr x te
    allTF_perm_total(iSub,:,:,:)=em.pSl.total; 
    clear em 
    
end

% create a t-test mask
tTestMaskOnOff=1;

if tTestMaskOnOff
    ttestAll = zeros(size(allTF_real_total,2),size(allTF_real_total,3),size(allTF_real_total,4));
    
    for iCross=1:size(allTF_real_total,2)
        for iTr=1:size(allTF_real_total,3)
            for iTe=1:size(allTF_real_total,4)
                ttestAll(iCross,iTr,iTe)=ttest(allTF_real_total(:,iCross,iTr,iTe),allTF_perm_total(:,iCross,iTr,iTe));
            end
        end
    end
    
end




% zero out any ns values

%
%allTF_real_total


   

% plot all train/test combos
h=figure('units','normalized','outerposition',[0 0 1 1]);
xNew=-500:2501/40:2000;
yNew=xNew;
cLims = [0 .004];
plotIdx = 0;
for k=[ 1,6,11,16,...
        2,5,9,13,...
        3,7,10,14,...
        4,8,12,15]
    
    plotIdx = plotIdx+1;
    
    % set up plot titles
    if      k==1; trainCond=1; testCond=1;
    elseif  k==2; trainCond=1; testCond=2;
    elseif  k==3; trainCond=1; testCond=3;
    elseif  k==4; trainCond=1; testCond=4;
    elseif  k==5; trainCond=2; testCond=1;
    elseif  k==6; trainCond=2; testCond=2;
    elseif  k==7; trainCond=2; testCond=3;
    elseif  k==8; trainCond=2; testCond=4;
    elseif  k==9; trainCond=3; testCond=1;
    elseif  k==10; trainCond=3; testCond=2;
    elseif  k==11; trainCond=3; testCond=3;
    elseif  k==12; trainCond=3; testCond=4;
    elseif  k==13; trainCond=4; testCond=1;
    elseif  k==14; trainCond=4; testCond=2;
    elseif  k==15; trainCond=4; testCond=3;
    elseif  k==16; trainCond=4; testCond=4;
    end
    
    if      trainCond==1; trainLabel = 'Spatial-Fix';
    elseif  trainCond==2; trainLabel = 'Color-Fix';
    elseif  trainCond==3; trainLabel = 'Spatial-Move';
    elseif  trainCond==4; trainLabel = 'Color-Move';
    end
    
    if      testCond==1; testLabel = 'Spatial-Fix';
    elseif  testCond==2; testLabel = 'Color-Fix';
    elseif  testCond==3; testLabel = 'Spatial-Move';
    elseif  testCond==4; testLabel = 'Color-Move';
    end
    
    tl = ['Train: ' trainLabel ' \\ ' 'Test: ' testLabel];
    
    % create plots
    %subplot(4,4,k)
    subplot(4,4,plotIdx)
    
    
    dataForPlot = squeeze(mean(allTF_real_total(:,k,:,:),1));
    
    if tTestMaskOnOff
        for iTr=1:size(dataForPlot,1)
            for iTe=1:size(dataForPlot,2)
                
                if squeeze(ttestAll(k,iTr,iTe))==0
                    dataForPlot(iTr,iTe)=0;
                end
                
            end
        end
    end
    
    imagesc(xNew,yNew,dataForPlot,cLims);% ,cLims
    
%    imagesc(xNew,yNew,squeeze(mean(allTF_real_total(:,k,:,:),1)),cLims) %cLims
    
    ylabel('Training')
    xlabel('Testing')
    pbaspect([1,1,1])
    title(tl)
    colorbar
    
    % target onset
    vline(0,'k--')
    hline(0,'k--')
    
    % target offset
    vline(250,'k--')
    hline(250,'k--')
    
    % first eye movement (R or L)
    vline(500,'w--')
    hline(500,'w--')
    
    % second eye movement (R or L)
    vline(1000,'w--')
    hline(1000,'w--')
    
    % third eye movement (back to fix)
    vline(1500,'w--')
    hline(1500,'w--')
    
    clear dataForPlot
    
    
end







% 
% xNew=-500:2501/40:2000;
% yNew=xNew;
% %cLims = [-.0003 .0007];
% %cLims = [0 .0007]
% 
% h=figure;
% 
% for i=1:8 %size(allTF_real_total,2)
% 
%     if i==2; tl='Eyes Closed'; m=13;
%     elseif i==1; tl='Eyes Open'; m=14;
%     elseif i==4; tl='Eyes Closed Masked'; m=15;
%     elseif i==3; tl='Eyes Open Masked'; m=16;
%     elseif i==6; tl='Train EC, Test EO'; m=1;
%     elseif i==5; tl='Train EO, Test EC'; m=4;
%     elseif i==8; tl='Train ECM, Test EOM'; m=9;
%     elseif i==7; tl='Train EOM, Test ECM'; m=12;
%     end
%     
%     subplot(2,4,i)
%     imagesc(xNew,yNew,squeeze(mean(allTF_real_total(:,m,:,:),1))) %cLims
%     ylabel('Training')
%     xlabel('Testing')
%     pbaspect([1,1,1])
%     title(tl)
%     cbar
%     
%     vline(0,'k--')
%     hline(0,'k--')
%     vline(250,'k--')
%     hline(250,'k--')
%     
%     if i==3||i==4||i==7||i==8
%     vline(300,'k--')
%     hline(300,'k--')
%     end
%     %pause(1)
%     
% end
% 
% % isolate 8 important comparisons (see above plot)
% %allTF_mat = (allTF_real_total(:,[13,14,15,16,1,4,9 12],:,:));
% 
% parsave([destDir '/' 'ALL_TT_WITHIN_BETWEEN.mat'],allTF_real_total,allTF_perm_total)
% 
% %% DISPLAY ALL 12 additional training/testing cross plots
% h=figure;
% for i=1:12
%     
%     if i==1; tl='TR-EC, TE-EO';
%     elseif i==2; tl='TR-EC, TE-ECM';
%     elseif i==3; tl='TR-EC, TE-EOM';
%     elseif i==4; tl='TR-EO, TE-EC';
%     elseif i==5; tl='TR-EO, TE-ECM';
%     elseif i==6; tl='TR-EO, TE-EOM';
%     elseif i==7; tl='TR-ECM, TE-EC';
%     elseif i==8; tl='TR-ECM, TE-EO';
%     elseif i==9; tl='TR-ECM, TE-EOM';
%     elseif i==10; tl='TR-EOM, TE-EC';
%     elseif i==11; tl='TR-EOM, TE-EO';
%     elseif i==12; tl='TR-EOM, TE-ECM';
%     end 
% 
%     subplot(4,3,i)
%     imagesc(xNew,yNew,squeeze(mean(allTF_real_total(:,i,:,:),1))) %cLims
%     ylabel('Training')
%     xlabel('Testing')
%     pbaspect([1,1,1])
%     title(tl)
%     %cbar
%     vline(0,'k--')
%     hline(0,'k--')
%     vline(250,'k--')
%     hline(250,'k--')
%    
% end