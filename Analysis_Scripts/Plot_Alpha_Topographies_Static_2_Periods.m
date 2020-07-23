%{
Plot_Alpha_Topographies
Author: Tom Bullock
Date: 04.22.20

To-Do: 

Figure out the normalizing (look at previous version of script in
more detail).

Make gifs to show alpha shifts over the course of the trial.

%}

% load EEGLAB (if needed)
% cd('/Users/tombullock/Documents/MATLAB/eeglab2019_1')
% eeglab
clear
close all
% cd /Users/tombullock/Documents/Psychology/WTF_EYE/Analysis_Scripts_Local

%set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
plotDir = '/home/waldrop/Desktop/WTF_EYE/Plots';

% load data
load([sourceDir '/' 'Bandpassed_Data_ALpha.mat'])

% plot settings
theseMapLimits = [0,1];%'maxmin';
showElects = 'off';

% create times vector
times = dataInfo.times;
chanlocs = dataInfo.chanlocs;

% set posterior channels for normalization
poElects = [25 26 30 62 63];
pElects = [22 21 20 31 57 58 59];
oElects = [27 29 64];
theseChans = [poElects,oElects,pElects];

% set time segments for analysis
thisTime=-500; % start baseline period
[~, baseStart] = min(abs(thisTime-times));
thisTime=0; % probe onset
[~, probeOn] = min(abs(thisTime-times));
thisTime=200; % probe
[~, probeFinal50ms] = min(abs(thisTime-times));
thisTime=250; % probe offset
[~, probeOff] = min(abs(thisTime-times));
thisTime=500; % EM1 cue
[~, EM1cue] = min(abs(thisTime-times));
thisTime=1000; % EM2 cue
[~, EM2cue] = min(abs(thisTime-times));
thisTime=1250;
[~, halfRet] = min(abs(thisTime-times));
thisTime=1500; % EM3 cue
[~, EM3cue] = min(abs(thisTime-times));
thisTime=2000; % end retention period
[~, retEnd] = min(abs(thisTime-times));

% baseline correct to mean of whole baseline period
%allBand = allBand - mean(allBand(:,:,:,:,baseStart:probeOn-1),5);

% convert bandpassed data from complex numbers to to power
%allBand = abs(allBand).^2;
    
% loop through and plot differnet time segments
for timeSegmentToPlot=1:2
    
    
    % time segment settings
    if timeSegmentToPlot==1
        thisTimeSegment = probeFinal50ms:probeOff;
        thisTimeSegmentLabel = 'Stimulus On (0.2-0.25 s)';
        thisTimeSegmentTitle = 'Probe';
    elseif timeSegmentToPlot==2
        thisTimeSegment = EM1cue:retEnd;
        thisTimeSegmentLabel = 'Retention (0.5-2 s)';
        thisTimeSegmentTitle = 'Retention';
%     elseif timeSegmentToPlot==3
%         thisTimeSegment = halfRet:retEnd;
%         thisTimeSegmentLabel = 'Late Retention (1.25-2s)';
%         thisTimeSegmentTitle = 'Late Retention';        
    end
    
    
    
    % condition loop
    h=figure('units','normalized','outerposition',[0.1 0.1 .9 1]);
    
    for iCond=1:4
        
        % loop through locations and plot
        for iLoc=[1:8];%[4,3,2,1,5,6,7,8]
            
            subplot(100,100,100)
 
            %% NOTE THAT EACH LOCATION IS INDEPENDENTLY NORMALIZED
            
            % find max and min across all locations for this condition/time segment
            thisMax = max(squeeze(mean(mean(mean(allBand(:,iCond,iLoc,theseChans,thisTimeSegment)),5),1)));
            thisMin = min(squeeze(mean(mean(mean(allBand(:,iCond,iLoc,theseChans,thisTimeSegment)),5),1)));
            
            % normalize the bandpassed data
            allBandNorm = (allBand-thisMin)/(thisMax-thisMin);
            
            % set data
            theseData = squeeze(mean(mean(allBandNorm(:,iCond,iLoc,:,thisTimeSegment),1),5));
            
            %headplot(theseData,'wtfSpline.spl','view','back','maplimits',theseMapLimits,'electrodes',showElects) %'cbar',0,
            ax(iLoc)=headplot(theseData,'wtfSpline.spl','view','back','electrodes',showElects,'maplimits',theseMapLimits) %'cbar',0,
            
            % set placement for each condition's topos
%             if iCond==1
%                 centerX=.3;
%                 centerY=.6;
%             elseif iCond==2
%                 centerX=.6;
%                 centerY=.6;
%             elseif iCond==3
%                 centerX=.3;
%                 centerY=.15;
%             elseif iCond==4
%                 centerX=.6;
%                 centerY=.15;
%             end
            
            if iCond==1
                centerX=.1;
                centerY=.6;
            elseif iCond==2
                centerX=.35;
                centerY=.6;
            elseif iCond==3
                centerX=.6;
                centerY=.6;
            elseif iCond==4
                centerX=.85;
                centerY=.6;
            end
            centerX = centerX-.05;
            
            % set figure sizes
            sizeX=.15;
            sizeY=.15;
%             
%             if      iLoc==1; xyShift = [.08,.07];
%             elseif  iLoc==2; xyShift = [.03,.14];
%             elseif  iLoc==3; xyShift = [-.03,.14];
%             elseif  iLoc==4; xyShift = [-.08,.07];
%             elseif  iLoc==5; xyShift = [-.08,-.07];
%             elseif  iLoc==6; xyShift = [-.03,-.14];
%             elseif  iLoc==7; xyShift = [.03,-.14];
%             elseif  iLoc==8; xyShift = [.08,-.07];
%             end

            thisX1 = .08;
            thisX2 = .03;
            thisY1 = .07;
            thisY2 = .14;
            
            if      iLoc==1; xyShift = [thisX1,thisY1];
            elseif  iLoc==2; xyShift = [thisX2,thisY2];
            elseif  iLoc==3; xyShift = [-thisX2,thisY2];
            elseif  iLoc==4; xyShift = [-thisX1,thisY1];
            elseif  iLoc==5; xyShift = [-thisX1,-thisY1];
            elseif  iLoc==6; xyShift = [-thisX2,-thisY2];
            elseif  iLoc==7; xyShift = [thisX2,-thisY2];
            elseif  iLoc==8; xyShift = [thisX1,-thisY1];
            end
            
            ax(iLoc).Position = [centerX+xyShift(1),centerY+xyShift(2),sizeX,sizeY];
            
            % add title to each plot
            %title(['Loc ' num2str(iLoc)],'FontSize',24)
            
        end
        
        % add condition labels to plot
        if      iCond==1; condTitle = 'S/F';
        elseif  iCond==2; condTitle = 'C/F';
        elseif  iCond==3; condTitle = 'S/M';
        elseif  iCond==4; condTitle = 'C/M';
        end
        %annotation('textbox', [centerX+.04,centerY+.11, 0,0], 'string', condTitle,'FontSize',36);
        
        clear allBandNorm
        
    end
    
    % add title
    if timeSegmentToPlot==1
        figureTitlePosition=[.25,.97,0,0];
    elseif timeSegmentToPlot==2
        figureTitlePosition=[.25,.97,0,0];
    elseif timeSegmentToPlot==3
        figureTitlePosition=[.25,.97,0,0];        
    end
    %annotation('textbox',figureTitlePosition, 'string', thisTimeSegmentLabel,'FontSize',36,'FitBoxToText','on','LineStyle','none');
    
    if timeSegmentToPlot==1
        saveas(h,[plotDir '/' 'Alpha_Topos_' thisTimeSegmentTitle '.eps'],'epsc')
    elseif timeSegmentToPlot==2
        saveas(h,[plotDir '/' 'Alpha_Topos_' thisTimeSegmentTitle '.eps'],'epsc')
%     elseif timeSegmentToPlot==3
%         saveas(h,[plotDir '/' 'Alpha_Topos_' thisTimeSegmentTitle '.eps'],'epsc')        
    end
    
    clear ax theseData allBandNorm
    
    
end