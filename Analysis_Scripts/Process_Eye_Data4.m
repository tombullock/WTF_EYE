%{
EYE_IEM_Analysis_Start
Author:Tom Bullock
Date: 04.20.20

Need to figure out the timing of the target onset (get from behavioral
file) then mark on plots the onset of the fixation movements then use these
as a guide for marking the actual eye-movements. This may need to be done
manually (or, semi manually) unfortunately.  Once we can identify the
moment at which the eye saccades to the new location, we can then parse
the brain activity from that moment onwards and feed into the IEM.  This
should cut out temporal variability that may be preventing the CTFs from
re-emerging in the regular analyses.
%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Sync_Final';
destDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Fixation_Mats';
plotDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Fixation_Plots';

% select subject and condition
sjNum=31;
thisCondition=3;


% load data
load([sourceDir '/' sprintf('sj%02d_eye_beh_sync_pro.mat',sjNum)])

% merge sessions
for iCond=1:4
    eyeStruct(iCond).ex = [eyeProcessed(1,iCond).ex;eyeProcessed(2,iCond).ex];
    eyeStruct(iCond).ey = [eyeProcessed(1,iCond).ey;eyeProcessed(2,iCond).ey];
    eyeStruct(iCond).pa = [eyeProcessed(1,iCond).pa;eyeProcessed(2,iCond).pa];
    eyeStruct(iCond).trialData = [eyeProcessed(1,iCond).trialData,eyeProcessed(2,iCond).trialData];
    eyeStruct(iCond).times = eyeProcessed(1,iCond).times;
end


% relabel
iCond=thisCondition;


% compute euclidian distance from baseline eye position and plot


% trial loop
nTrials = length(eyeStruct(iCond).trialData);

% get times
times = eyeStruct.times;

for iTrial = 239:nTrials
    
    figure('units','normalized','outerposition',[0 0 2 2])
    
    
    % get 200 ms pre-target baseline
    bx = nanmean(eyeStruct(iCond).ex(iTrial,find(times==-250):find(times==0)));
    by = nanmean(eyeStruct(iCond).ey(iTrial,find(times==-250):find(times==0)));
    
    % get eye data for trial
    x = eyeStruct(iCond).ex(iTrial,:);
    y = eyeStruct(iCond).ey(iTrial,:);
    
    % compute euclidian dist
    for j=1:length(x)
        e(j) = sqrt(((bx-x(j))^2) + ((by-y(j))^2));
    end
    
    % calculate differential
    for t=1:length(e)-1
        d(t) = e(t+1)-e(t);
    end
    
    % create plots
    subplot(3,1,1)
    px=plot(times,x,'LineWidth',4); %hold on
    
    ax = gca;
    ax.XLim = [-250,2000];
    line([500,500],ax.YLim,'linestyle','--','color','k')
    line([1000,1000],ax.YLim,'linestyle','--','color','k')
    line([1500,1500],ax.YLim,'linestyle','--','color','k')
    ylabel('Y Position')
    
    subplot(3,1,2)
    py=plot(times,y,'LineWidth',4); %hold on
    
    ay = gca;
    ay.XLim = [-250,2000];
    line([500,500],ay.YLim,'linestyle','--','color','k')
    line([1000,1000],ay.YLim,'linestyle','--','color','k')
    line([1500,1500],ay.YLim,'linestyle','--','color','k')
    ylabel('Y Position')
    
    % plot differential
    subplot(3,1,3)
    pd = plot(times(1:length(times)-1),d,'LineWidth',2);
    
    ad = gca;
    ad.XLim = [-250,2000];
    line([500,500],ad.YLim,'linestyle','--','color','k')
    line([1000,1000],ad.YLim,'linestyle','--','color','k')
    line([1500,1500],ad.YLim,'linestyle','--','color','k')
    ylabel('Euclidian Distance Differential')
    
    % add interative drawing
    fixObj.fix1_start = drawpoint('color','g','label','F1','LabelVisible','Hover');
    fixObj.fix1_end = drawpoint('color','r','label','F1','LabelVisible','Hover');
    
    fixObj.fix2_start = drawpoint('color','g','label','F2','LabelVisible','Hover');
    fixObj.fix2_end = drawpoint('color','r','label','F2','LabelVisible','Hover');
    
    fixObj.fix3_start = drawpoint('color','g','label','F3','LabelVisible','Hover');
    fixObj.fix3_end = drawpoint('color','r','label','F3','LabelVisible','Hover');
    
    
    
    % replot the x and y subplots with patches to reflet eye fixes
    goodPlacement=0;
    cnt=0;
    while goodPlacement==0
        
        % create plots
        subplot(3,1,1)
        px=plot(times,x,'LineWidth',4); %hold on
        
        ax = gca;
        ax.XLim = [-250,2000];
        line([500,500],ax.YLim,'linestyle','--','color','k')
        line([1000,1000],ax.YLim,'linestyle','--','color','k')
        line([1500,1500],ax.YLim,'linestyle','--','color','k')
        ylabel('X Position')
        
        % title
        title(['Trial ' num2str(iTrial) '  U = update, C = confirm, B = bad trial'],'FontSize',36)
        
        % plot patch1
        fix1_start = fixObj.fix1_start.Position(1);
        fix1_end = fixObj.fix1_end.Position(1);
        patchX=[fix1_start,fix1_end,fix1_end,fix1_start];
        patchY = [ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        patch(patchX,patchY,'g','FaceAlpha',.2);
        
        % plot patch2
        fix2_start = fixObj.fix2_start.Position(1);
        fix2_end = fixObj.fix2_end.Position(1);
        patchX=[fix2_start,fix2_end,fix2_end,fix2_start];
        patchY = [ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        patch(patchX,patchY,'g','FaceAlpha',.2);
        
        % plot patch3
        fix3_start = fixObj.fix3_start.Position(1);
        fix3_end = fixObj.fix3_end.Position(1);
        patchX=[fix3_start,fix3_end,fix3_end,fix3_start];
        patchY = [ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        patch(patchX,patchY,'g','FaceAlpha',.2);
        
        % create plots
        subplot(3,1,2)
        py=plot(times,y,'LineWidth',4); %hold on
        
        ay = gca;
        ay.XLim = [-250,2000];
        line([500,500],ay.YLim,'linestyle','--','color','k')
        line([1000,1000],ay.YLim,'linestyle','--','color','k')
        line([1500,1500],ay.YLim,'linestyle','--','color','k')
        ylabel('Y Position')
        
        % plot patch1
        fix1_start = fixObj.fix1_start.Position(1);
        fix1_end = fixObj.fix1_end.Position(1);
        patchX=[fix1_start,fix1_end,fix1_end,fix1_start];
        patchY = [ay.YLim(1),ay.YLim(1),ay.YLim(2),ay.YLim(2)];
        patch_fix1 = patch(patchX,patchY,'g','FaceAlpha',.2);
        
        % plot patch2
        fix2_start = fixObj.fix2_start.Position(1);
        fix2_end = fixObj.fix2_end.Position(1);
        patchX=[fix2_start,fix2_end,fix2_end,fix2_start];
        patchY = [ay.YLim(1),ay.YLim(1),ay.YLim(2),ay.YLim(2)];
        patch(patchX,patchY,'g','FaceAlpha',.2);
        
        % plot patch3
        fix3_start = fixObj.fix3_start.Position(1);
        fix3_end = fixObj.fix3_end.Position(1);
        patchX=[fix3_start,fix3_end,fix3_end,fix3_start];
        patchY = [ay.YLim(1),ay.YLim(1),ay.YLim(2),ay.YLim(2)];
        patch(patchX,patchY,'g','FaceAlpha',.2);
        
        % get fig
        fig=gcf;
        
        % check the RA is happy with the placement of the items
        disp('PRESS "C" TO CONFIRM TRIAL.  PRESS "N" TO RE-DO.')
        pressloop=1;
        badTrial=0;
        while pressloop
            was_a_key = waitforbuttonpress;
            if was_a_key && strcmp(get(fig, 'CurrentKey'), 'c')
                pressloop=0;
                goodPlacement=1;
                badTrial=0;
            elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'u')
                pressloop=0;
                goodPlacement=0;
            elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'b')
                pressloop=0;
                goodPlacement=1;
                badTrial=1;
            end
        end
        
    end
    
    fixStruct(iTrial).fix1_start = round(fix1_start);
    fixStruct(iTrial).fix1_end = round(fix1_end);
    fixStruct(iTrial).fix2_start = round(fix2_start);
    fixStruct(iTrial).fix2_end = round(fix2_end);
    fixStruct(iTrial).fix3_start = round(fix3_start);
    fixStruct(iTrial).fix3_end = round(fix3_end);
    fixStruct(iTrial).badTrial = badTrial;
    
    % clear a bunch of stuff
    clear allClean allClean_idx allBlocked_idx postEM_idx preEM_idx
    clear em1_idx fix1_start fix1_end fix2_start fix2_end fix3_start fix3_end patchX patchY
    
    % save figure
    %saveas(fig,[plotDir '/' sprintf('sj%02d_cond%02d_trial%02d.jpeg',sjNum,iCond,iTrial)],'jpeg')
    
    % close figure
    close
    
end

% save trial data [corrected, was wrong for sj31cond03
trialData = eyeStruct(iCond).trialData;


% save data
save([destDir '/' sprintf('sj%02d_cond%02d.mat',sjNum,iCond)],'fixStruct','trialData','times')
