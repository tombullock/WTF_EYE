%{
Plot_Trials_Experiment_Schematic_With_Alpha_Topos
Author: Tom Bullock
Date: 04.29.20

To-Do: 

Figure out the normalizing (look at previous version of script in
more detail).

Make gifs to show alpha shifts over the course of the trial.

%}

% load EEGLAB (if needed)
%cd('/Users/tombullock/Documents/MATLAB/eeglab2019_1')
%eeglab
clear
close all
%cd /Users/tombullock/Documents/Psychology/WTF_EYE/Analysis_Scripts_Local

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

% loop through locations and plot
for iLoc=[1,4,8,5]
    for thisTimeIdx = baseStart:4:retEnd
        
        % open figure
        h=figure('units','normalized','outerposition',[0.2 0.1 .8 .5]);
        
        % each frame averages over 4 samples
        thisTimeSegment = thisTimeIdx:thisTimeIdx+3;
        
        plotCnt=plotCnt+1;
        
        for iCond=1%:4
            
            subplot(1,2,1)
            
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
            
            %ax(iLoc).Position = [centerX+xyShift(1),centerY+xyShift(2),sizeX,sizeY];
            
            % add title to each plot
            %title(['Loc ' num2str(iLoc)],'FontSize',24)
            
            % plot the experiment screen
            subplot(1,2,2)
            
            % draw rectangle (monitor frame)
            rectangle('Position',[0 0 1 1],'Curvature',0.2,'LineWidth',6,'FaceColor',[128,128,128]./255);
            
            axis off
            
            % draw fixation dot
            f = rectangle('Position',[.48,.48,.04,.05],'Curvature',1,'FaceColor',[60,60,60]/255);
            
            % draw probe
            probeRightHighPos = [.65,.65,.16,.16];
            probeLeftHighPos = [.2,.65,.16,.16];
            probeLeftLowPos = [.2,.2,.16,.16];
            probeRightLowPos = [.65,.2,.16,.16];
            if      iLoc==1; probePos = probeRightHighPos;
            elseif  iLoc==4; probePos = probeLeftHighPos;
            elseif  iLoc==5; probePos = probeLeftLowPos;
            elseif  iLoc==8; probePos = probeRightLowPos;
            end
            
            % generate probe shape at appropriate time in trial
            if times(thisTimeIdx)>0 && times(thisTimeIdx)<250
                p = rectangle('Position',probePos,'Curvature',1,'FaceColor',[1,1,255]/255);
            end
            
            
            %end
            
            %         % add condition labels to plot
            %         if      iCond==1; condTitle = 'S/F';
            %         elseif  iCond==2; condTitle = 'C/F';
            %         elseif  iCond==3; condTitle = 'S/M';
            %         elseif  iCond==4; condTitle = 'C/M';
            %         end
            %         annotation('textbox', [centerX+.04,centerY+.11, 0,0], 'string', condTitle,'FontSize',48);
            %
            %         clear allBandNorm
            %         clear ax theseData allBandNorm
        end
        
        
        
        % add time counter
        thisTimeLabel = round(times(thisTimeSegment(1)));
        an = annotation('textbox',[.38,.55,0,0], 'string', [num2str(thisTimeLabel) ' ms'],'FontSize',48,'FitBoxToText','on','LineStyle','none');
        
        % add frame to gif
        if plotCnt==1
            gif([gifDir '/' 'Trial_Example_With_Topo.gif'],'frame',gcf,'DelayTime',1/15);
        else
            gif('frame',gcf);
        end
        
        % clear annotations for next plot iter
        delete(findall(gcf,'type','annotation'))
        
        close
        
    end
    
    % plot response screen
    for ITI=1:15
        
        % plot figure
        h=figure('units','normalized','outerposition',[0.2 0.1 .8 .5]);
        
        % plot the experiment screen
        subplot(1,2,2)
        
        % draw rectangle (monitor frame)
        rectangle('Position',[0 0 1 1],'Curvature',0.2,'LineWidth',6,'FaceColor',[128,128,128]./255);
        
        axis off
        
        % draw fixation dot
        f = rectangle('Position',[.48,.48,.04,.05],'Curvature',1,'FaceColor',[60,60,60]/255);

        % draw response wheel
        rw = rectangle('Position',[.2,.2,.6,.6],'Curvature',1,'EdgeColor',[64,64,64]/255,'LineWidth',12);
        
        % add frame to gif
        gif('frame',gcf);
        
        close
        
    end
    
    % plot ITI
    for ITI=1:10
        
        % plot figure
        h=figure('units','normalized','outerposition',[0.2 0.1 .8 .5]);
        
        % plot the experiment screen
        subplot(1,2,2)
        
        % draw rectangle (monitor frame)
        rectangle('Position',[0 0 1 1],'Curvature',0.2,'LineWidth',6,'FaceColor',[128,128,128]./255);
        
        axis off
        
        % draw fixation dot
        %f = rectangle('Position',[.48,.48,.04,.05],'Curvature',1,'FaceColor',[60,60,60]/255);

        % add frame to gif
        gif('frame',gcf);
        
        close
        
    end
    
end