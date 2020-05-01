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
cd('/home/waldrop/Desktop/WTF_EYE/Dependancies/eeglab2019_1')
eeglab
clear
close all
cd /home/waldrop/Desktop/WTF_EYE/Analysis_Scripts

%set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
gifDir = '/home/waldrop/Desktop/WTF_EYE/Gifs';

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
thisTime=1500; % EM3 cue
[~, EM3cue] = min(abs(thisTime-times));
thisTime=2000; % end retention period
[~, retEnd] = min(abs(thisTime-times));

% baseline correct to mean of whole baseline period
allBand = allBand - mean(allBand(:,:,:,:,baseStart:probeOn-1),5);

% convert bandpassed data from complex numbers to to power
allBand = abs(allBand).^2;
    




% % time segment settings
% if timeSegmentToPlot==1
%     thisTimeSegment = probeFinal50ms:probeOff;
%     thisTimeSegmentLabel = 'Probe (200-250ms)';
%     thisTimeSegmentTitle = 'Probe';
% elseif timeSegmentToPlot==2
%     thisTimeSegment = EM1cue:retEnd;
%     thisTimeSegmentLabel = 'Retention (500 - 2000 ms)';
%     thisTimeSegmentTitle = 'Retention';
% end

plotCnt=0;
for thisTimeIdx = baseStart:4:retEnd
    
    % open figure
    h=figure('units','normalized','outerposition',[0.1 0.1 .9 1]);
    
    % each frame averages over 4 samples
    thisTimeSegment = thisTimeIdx:thisTimeIdx+3;
    
    plotCnt=plotCnt+1;
    
    for iCond=1:4
        
        % loop through locations and plot
        for iLoc=1:8
            
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
            if iCond==1
                centerX=.3;
                centerY=.6;
            elseif iCond==2
                centerX=.6;
                centerY=.6;
            elseif iCond==3
                centerX=.3;
                centerY=.15;
            elseif iCond==4
                centerX=.6;
                centerY=.15;
            end
            
            % set figure sizes
            sizeX=.15;
            sizeY=.15;
            
            if      iLoc==1; xyShift = [.08,.07];
            elseif  iLoc==2; xyShift = [.03,.14];
            elseif  iLoc==3; xyShift = [-.03,.14];
            elseif  iLoc==4; xyShift = [-.08,.07];
            elseif  iLoc==5; xyShift = [-.08,-.07];
            elseif  iLoc==6; xyShift = [-.03,-.14];
            elseif  iLoc==7; xyShift = [.03,-.14];
            elseif  iLoc==8; xyShift = [.08,-.07];
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
        annotation('textbox', [centerX+.04,centerY+.11, 0,0], 'string', condTitle,'FontSize',48);
        
        clear allBandNorm
        clear ax theseData allBandNorm
    end
    

    
    % add time counter 
    
    thisTimeLabel = round(times(thisTimeSegment(1)));
    
    
    an = annotation('textbox',[.45,.47,0,0], 'string', [num2str(thisTimeLabel) ' ms'],'FontSize',48,'FitBoxToText','on','LineStyle','none');
    
    if plotCnt==1
        gif([gifDir '/' 'Alpha_Topos_SCAM.gif'],'frame',gcf,'DelayTime',1/5);
    else
        gif('frame',gcf);
    end
    
    % clear annotations for next plot iter
    delete(findall(gcf,'type','annotation'))
    
    close
    
end




% if timeSegmentToPlot==1
%     saveas(h,[plotDir '/' 'Alpha_Topos_' thisTimeSegmentTitle '.png'],'png')
% elseif timeSegmentToPlot==2
%     saveas(h,[plotDir '/' 'Alpha_Topos_' thisTimeSegmentTitle '.png'],'png')
% end




