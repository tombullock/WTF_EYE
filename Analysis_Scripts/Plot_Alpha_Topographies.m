%{
Plot_Alpha_Topographies
Author: Tom Bullock
Date: 04.22.20

%}

clear
close all

%set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% load data
load([sourceDir '/' 'Bandpassed_Data_ALpha.mat'])

% adding path
addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Dependancies/eeglab14_1_2b/functions/resources'))

% create times vector
times = linspace(-500,2500,768);

% set occipital channels for norm
%chans = [23 22 21 25 26 27 20 29 30 31    32+25, 32+31, 32+32, 32+26, 32+30, 32+27, 32+28] ;
poElects = [25 26 30 62 63    22 21 20 31 57 58 59]; % AND P-ELECTS
oElects = [27 29 64];
theseChans = [poElects,oElects];

% baseline correct
allBand = allBand - mean(allBand(:,:,:,:,128));

% convert to power
%allBand = abs(allBand).^2;

% condition loop
for iCond=1
    
    thisTimeSegment = 129:192; % 0-250 ms
    
    % may need to do this for each loc separately???
    thisMax = squeeze(max(max(max(max(allBand(:,iCond,:,theseChans,thisTimeSegment))))));
    thisMin = squeeze(min(min(min(min(allBand(:,iCond,:,theseChans,thisTimeSegment))))));

    
    % normalize the bandpassed data
    allBandNorm = (allBand-thisMin)/(thisMax-thisMin);
    
    % loop through locations and plot
    for iLoc=1:8
        
        % Tom up to here - need to create a new spline file then loop through
        % and plot each loc as a headplot, for each condition (04.23.20)
        
        
        theseMapLimits = [0 1];%'maxmin'; 
        showElects = 'off';
        
        load('/home/waldrop/Desktop/WTF_EYE/Dependancies/eeglab14_1_2b/functions/resources/mheadnew.mat')
        
        theseData = squeeze(mean(mean(allBandNorm(:,iCond,iLoc,:,thisTimeSegment),1),5));
        
        headplot(theseData,'wtfSpline.spl','view','back','electrodes',showElects,'maplimits',theseMapLimits) %'cbar',0,,
        
    
    
    end
    
    
end






% 
% 
% ***************************************************************************
% 
% 
% 
% %% load data
% load('/home/bullock/WTF_DATA_FOSTER_SCRIPTS/Data_Compiled/ALL_HILBERT.mat')
% 
% %% plot dir
% plotDir = '/home/bullock/WTF_DATA_FOSTER_SCRIPTS/TOPO_PLOTS';
% 
% %% plot target topos (1) or retention topos (2)?
% plotTimePeriod = 1;
% 
% %% baseline correct hilbert mat (18*4*64*768)
% hilbMatBL = mean(hilbMat(:,:,:,1:128),4);
% 
% for i=1:size(hilbMat,1)
%     for j=1:size(hilbMat,2)
%         for k=1:size(hilbMat,3)
%             for m=1:size(hilbMat,4)
%             
%                 hilbMatBC(i,j,k,m) = hilbMat(i,j,k,m) - hilbMatBL(i,j,k);
%             
%             end
%         end
%     end
% end
% 
% %% use baseline corrected data
% hilbMat =[];
% hilbMat = hilbMatBC;
% 
% 
% 
% 
% 
% 
% %% define channels for plotting average hilbert waves
% chans = [23 22 21 25 26 27 20 29 30 31    32+25, 32+31, 32+32, 32+26, 32+30, 32+27, 32+28] ;
% 
% times = [1:640];
% 
% 
% % errorbar(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,1,chans,times),1),3)), std(squeeze(mean(hilbMat(:,1,chans,times),3)))./sqrt(size(hilbMat,1)),'r'); hold on
% % errorbar(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,2,chans,times),1),3)), std(squeeze(mean(hilbMat(:,2,chans,times),3)))./sqrt(size(hilbMat,1)),'c')
% % errorbar(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,3,chans,times),1),3)), std(squeeze(mean(hilbMat(:,3,chans,times),3)))./sqrt(size(hilbMat,1)),'g')
% % errorbar(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,4,chans,times),1),3)), std(squeeze(mean(hilbMat(:,4,chans,times),3)))./sqrt(size(hilbMat,1)),'m')
% h=figure;
% 
% title('Alpha Power (induced, Hilbert, O, PO, P electrodes')
% plot(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,1,chans,times),1),3)),'k','LineWidth',2); hold on % ec
% plot(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,2,chans,times),1),3)),'c','LineWidth',2) % eo
% plot(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,3,chans,times),1),3)),'g','LineWidth',2) % ec-m
% plot(-500:2501/640:2000,squeeze(mean(mean(hilbMat(:,4,chans,times),1),3)),'m','LineWidth',2) % eo-m
% 
% xlabel('time (ms)')
% ylabel('power uV"^2')
% set(gca,'LineWidth',2, 'box','off','FontSize',20)
% 
% %% plot topos for each location (generated spline file already)
% 
% if plotTimePeriod ==1
%     theseSamples = 179:192; % 198-250ms post stim onset TARGET
% else
%     theseSamples = 320:640; % 750 - 2000 ms post stim onset RETENTION (1000=384 idx)
% end
% 
% for thisCond=1:4 % EC, EO, EC-M, EO-M
%     % which condition?
%     %thisCond = 4;
%     
%     % try putting all on the same scale (-13 to 13 perhaps)?
%     %,'maplimits',[-6 6]
%     
%     thisEEG1=[]; thisEEG2=[]; thisEEG3=[]; thisEEG4=[]; thisEEG5=[]; thisEEG6=[]; thisEEG7=[]; thisEEG8=[];
%     
%     % define data of interest
%     thisEEG1 = squeeze(mean(mean(abs(locsMat.loc1(:,thisCond,:,theseSamples)).^2,4),1));
%     thisEEG2 = squeeze(mean(mean(abs(locsMat.loc2(:,thisCond,:,theseSamples)).^2,4),1));
%     thisEEG3 = squeeze(mean(mean(abs(locsMat.loc3(:,thisCond,:,theseSamples)).^2,4),1));
%     thisEEG4 = squeeze(mean(mean(abs(locsMat.loc4(:,thisCond,:,theseSamples)).^2,4),1));
%     thisEEG5 = squeeze(mean(mean(abs(locsMat.loc5(:,thisCond,:,theseSamples)).^2,4),1));
%     thisEEG6 = squeeze(mean(mean(abs(locsMat.loc6(:,thisCond,:,theseSamples)).^2,4),1));
%     thisEEG7 = squeeze(mean(mean(abs(locsMat.loc7(:,thisCond,:,theseSamples)).^2,4),1));
%     thisEEG8 = squeeze(mean(mean(abs(locsMat.loc8(:,thisCond,:,theseSamples)).^2,4),1));
%     
%     
%     % normalize amps between 0 and 1, based only on PO and O elecs
%     %pElects = [22 21 20 31 57 58 59 ];
%     poElects = [25 26 30 62 63    22 21 20 31 57 58 59]; % AND P-ELECTS
%     oElects = [27 29 64];
%     maxVal = max([thisEEG1([oElects, poElects]), thisEEG2([oElects, poElects]),thisEEG3([oElects, poElects]),thisEEG4([oElects, poElects]),...
%         thisEEG5([oElects, poElects]), thisEEG6([oElects, poElects]),thisEEG7([oElects, poElects]),thisEEG8([oElects, poElects])]);
%     minVal = min([thisEEG1([oElects, poElects]), thisEEG2([oElects, poElects]),thisEEG3([oElects, poElects]),thisEEG4([oElects, poElects]),...
%         thisEEG5([oElects, poElects]), thisEEG6([oElects, poElects]),thisEEG7([oElects, poElects]),thisEEG8([oElects, poElects])]);
%     
%     thisEEG1 = (thisEEG1 - minVal(1))/(maxVal(1)-minVal(1));
%     thisEEG2 = (thisEEG2 - minVal(2))/(maxVal(2)-minVal(2));
%     thisEEG3 = (thisEEG3 - minVal(3))/(maxVal(3)-minVal(3));
%     thisEEG4 = (thisEEG4 - minVal(4))/(maxVal(4)-minVal(4));
%     thisEEG5 = (thisEEG5 - minVal(5))/(maxVal(5)-minVal(5));
%     thisEEG6 = (thisEEG6 - minVal(6))/(maxVal(6)-minVal(6));
%     thisEEG7 = (thisEEG7 - minVal(7))/(maxVal(7)-minVal(7));
%     thisEEG8 = (thisEEG8 - minVal(8))/(maxVal(8)-minVal(8));
%     
%     
%     theseMapLimits = [0 1];%'maxmin'; % [0 8] works fairly well for most, but not great...
%     showElects = 'off';
%     
% %     % plot data
% %     h=figure;
% %     subplot(2,4,4)
% %     headplot(thisEEG1,'wtfSpline.spl','view','back','cbar',0,'maplimits',theseMapLimits,'electrodes',showElects)
% %     title('0-45 degs')
% %     subplot(2,4,3)
% %     headplot(thisEEG2,'wtfSpline.spl','view','back','cbar',0,'maplimits',theseMapLimits,'electrodes',showElects)
% %     title('46-90 degs')
% %     subplot(2,4,2)
% %     headplot(thisEEG3,'wtfSpline.spl','view','back','cbar',0,'maplimits',theseMapLimits,'electrodes',showElects)
% %     title('91-135 degs')
% %     subplot(2,4,1)
% %     headplot(thisEEG4,'wtfSpline.spl','view','back','cbar',0','maplimits',theseMapLimits,'electrodes',showElects)
% %     title('136-180 degs')
% %     
% %     subplot(2,4,5)
% %     headplot(thisEEG5,'wtfSpline.spl','view','back','cbar',0,'maplimits',theseMapLimits,'electrodes',showElects)
% %     title('181-225 degs')
% %     subplot(2,4,6)
% %     headplot(thisEEG6,'wtfSpline.spl','view','back','cbar',0,'maplimits',theseMapLimits,'electrodes',showElects)
% %     title('226-270 degs')
% %     subplot(2,4,7)
% %     headplot(thisEEG7,'wtfSpline.spl','view','back','cbar',0,'maplimits',theseMapLimits,'electrodes',showElects)
% %     title('271-315 degs')
% %     subplot(2,4,8)
% %     headplot(thisEEG8,'wtfSpline.spl','view','back','cbar',0,'maplimits',theseMapLimits,'electrodes',showElects)
% %     title('316-360 degs')
%     
%     
%     % plot individual heads to figs for paper/AI
%     for i=1:8
%         
%         if i==1; theseData = thisEEG1;
%         elseif i==2; theseData = thisEEG2;
%         elseif i==3; theseData = thisEEG3;
%         elseif i==4; theseData = thisEEG4;
%         elseif i==5; theseData = thisEEG5;
%         elseif i==6; theseData = thisEEG6;
%         elseif i==7; theseData = thisEEG7;
%         elseif i==8; theseData = thisEEG8;
%         end
%         
%         if plotTimePeriod==1; thisName = 'targ';
%         else thisName = 'ret';
%         end
%           
%         h=figure;
%         headplot(theseData,'wtfSpline.spl','view','back','maplimits',theseMapLimits,'electrodes',showElects) %'cbar',0,
%         set(gca,'LooseInset',get(gca,'TightInset'))
%         saveas(h,[plotDir '/' num2str(thisName) '_' 'cond' num2str(thisCond) '_' 'topoLoc' num2str(i) '.eps'],'epsc')
%         close all
%         %saveas(h,[plotDir '/' num2str(thisName) '_' 'cond' num2str(thisCond) '_' 'topoLoc' num2str(i) '.pdf'],'pdf')
%     end
%     
%     
%     
% end

