%{
EYE_Movement_Analysis
Author: Tom Bullock, UCSB Attention Lab
Date: 06.04.20

Compute metric to determine whether subjects moved their eyes to (or close
to) the fixation dot in each new location.  

Find Euclidian Distance of closest eye gaze point fixation dot on for each
eye movement cue (center>newLoc>newLoc>center).  Maybe just plot for the
two newLoc positions?  Calculate an average euclidian distance error over
both positions and over all trials.  Average over subs and compute error
bars.  Done.

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Sync_Final';
destDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Euclidian_Error_mats';
%plotDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Fixation_Plots';

% select subject and condition
sjNum=2;
%thisCondition=3;

doPlots=0;


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

% loop through conditions
for iCond=3:4
    
    % % relabel
    % iCond=thisCondition;
    
    
    % trial loop
    nTrials = length(eyeStruct(iCond).trialData);
    
    % get times
    times = eyeStruct.times;
    
    % loop through trials
    for iTrial=1:nTrials
        
        % get coords for fixation dot
        xFix = eyeStruct(iCond).trialData(iTrial).eyeLocX;
        yFix = eyeStruct(iCond).trialData(iTrial).eyeLocY;
        
        % screen center coords
        centerX=1024/2;
        centerY=768/2;
        
        % get fixation dot new "move" pixel locations
        fixLoc1 = [(centerX + xFix), (centerY + yFix)];
        fixLoc2 = [(centerX - xFix), (centerY - yFix)];
        
        % get eye data for trial
        ex = eyeStruct(iCond).ex(iTrial,:);
        ey = eyeStruct(iCond).ey(iTrial,:);
        
        % get eye movements for each movement window
        theseTimes = [500,1000];
        theseTimesIdx = (find(times==theseTimes(1)):find(times==theseTimes(2)));
        
        ex1 = ex(theseTimesIdx);
        ey1 = ey(theseTimesIdx);
        
        theseTimes = [1000,1500];
        theseTimesIdx = (find(times==theseTimes(1)):find(times==theseTimes(2)));
        
        ex2 = ex(theseTimesIdx);
        ey2 = ey(theseTimesIdx);
        
        % get means over em1 and em2 windows
        mean_em1 = [mean(ex1),mean(ey1)];
        mean_em2 = [mean(ex2),mean(ey2)];
        
        
        % get euclidian distances of em1/em2 from fixLoc1/fixLoc2
        for j=1:length(ex1)
            
            % compare em1 to both pairs of fixLocs
            e1(j) = sqrt((fixLoc1(1)-ex1(j))^2 + (fixLoc1(2)-ey1(j))^2);
            e2(j) = sqrt((fixLoc2(1)-ex1(j))^2 + (fixLoc2(2)-ey1(j))^2);
            
            % compare em2 to both paris of fixLocs
            e3(j) = sqrt((fixLoc1(1)-ex2(j))^2 + (fixLoc1(2)-ey2(j))^2);
            e4(j) = sqrt((fixLoc2(1)-ex2(j))^2 + (fixLoc2(2)-ey2(j))^2);
            
            
        end
        
        % find min euclidian distance of EMs to fixLoc1 and fixLoc2
        [min_em1_fixLoc1, min_em1_fixLoc1_idx] = min(e1);
        [min_em1_fixLoc2, min_em1_fixLoc2_idx] = min(e2);
        [min_em2_fixLoc1, min_em2_fixLoc1_idx] = min(e3);
        [min_em2_fixLoc2, min_em2_fixLoc2_idx] = min(e4);
        
        %
        if min_em1_fixLoc1<min_em1_fixLoc2
            firstFixLoc = 1;
            fixLoc1_color = 'g';
            this_em1_ed = min_em1_fixLoc1;
            this_em1_pix = [ex1(min_em1_fixLoc1_idx),ey1(min_em1_fixLoc1_idx)];
        else
            firstFixLoc = 2;
            fixLoc1_color = 'b';
            this_em1_ed = min_em1_fixLoc2;
            this_em1_pix = [ex1(min_em1_fixLoc2_idx),ey1(min_em1_fixLoc2_idx)];
        end
        
        if firstFixLoc == 1
            fixLoc2_color = 'b';
            this_em2_ed = min_em2_fixLoc2;
            this_em2_pix = [ex2(min_em2_fixLoc2_idx),ey2(min_em2_fixLoc2_idx)];
        else
            fixLoc2_color = 'g';
            this_em2_ed = min_em2_fixLoc1;
            this_em2_pix = [ex2(min_em2_fixLoc1_idx),ey2(min_em2_fixLoc1_idx)];
        end
        
        % save euclidian distance errors to a struct
        euclidErrorStruct(iCond).em1(iTrial) = this_em1_ed;
        euclidErrorStruct(iCond).em2(iTrial) = this_em2_ed;
        
        
        
        
        %     if strcmp(fixLoc1_color,'g')
        %         fixLoc2_color = 'b';
        %
        %     else
        %         fixLoc2_color = 'g';
        %         this_em2 = min_em2_fixLoc2;
        %     end
        
        
        if doPlots==1
            
            % plot fixation dot locs
            plot(fixLoc1(1),fixLoc1(2),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor',fixLoc1_color,'MarkerSize',12,'LineWidth',3); hold on
            plot(fixLoc2(1),fixLoc2(2),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor',fixLoc2_color,'MarkerSize',12,'LineWidth',3);
            
            %     % plot mean ems
            %     plot(mean_em1(1),mean_em1(2),'LineStyle','none','Marker','o','MarkerFaceColor','g','MarkerSize',16);
            %     plot(mean_em2(1),mean_em2(2),'LineStyle','none','Marker','o','MarkerFaceColor','b','MarkerSize',16);
            
            % plot min ems
            plot(this_em1_pix(1),this_em1_pix(2),'LineStyle','none','Marker','o','MarkerFaceColor','g','MarkerSize',16);
            plot(this_em2_pix(1),this_em2_pix(2),'LineStyle','none','Marker','o','MarkerFaceColor','b','MarkerSize',16);
            
            % plot central fixation
            plot(centerX,centerY,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',12,'LineWidth',3);
            
            set(gca,'xlim',[0,1024],'ylim',[0,768])
            
            pause(5)
            close;
            
        end
        
        clear e1 e2 e3 e4 ex ex1 ex2 ey ey1 ey2 
        clear mean_em1 mean_em2 
        clear min_em1_fixLoc1_idx min_em1_fixLoc2_idx min_em2_fixLoc1_idx min_em2_fixLoc2_idx
        clear min_em1_fixLoc1 min_em1_fixLoc2 min_em2_fixLoc1 min_em2_fixLoc2
        clear this_em1_ed this_em1_pix this_em2_ed this_em2_pix
        clear xFix yFix
        
    end
    
    
    
end

save([destDir '/' sprintf('sj%02d_euclid_error.mat',sjNum)],'euclidErrorStruct')





% 
% 
% 
% for iTrial = 239:nTrials
%     
%     figure('units','normalized','outerposition',[0 0 2 2])
%     
%     
%     % get 200 ms pre-target baseline
%     bx = nanmean(eyeStruct(iCond).ex(iTrial,find(times==-250):find(times==0)));
%     by = nanmean(eyeStruct(iCond).ey(iTrial,find(times==-250):find(times==0)));
%     
%     % get eye data for trial
%     x = eyeStruct(iCond).ex(iTrial,:);
%     y = eyeStruct(iCond).ey(iTrial,:);
%     
%     % compute euclidian dist
%     for j=1:length(x)
%         e(j) = sqrt(((bx-x(j))^2) + ((by-y(j))^2));
%     end
%     
%     % calculate differential
%     for t=1:length(e)-1
%         d(t) = e(t+1)-e(t);
%     end
%     
%     % create plots
%     subplot(3,1,1)
%     px=plot(times,x,'LineWidth',4); %hold on
%     
%     ax = gca;
%     ax.XLim = [-250,2000];
%     line([500,500],ax.YLim,'linestyle','--','color','k')
%     line([1000,1000],ax.YLim,'linestyle','--','color','k')
%     line([1500,1500],ax.YLim,'linestyle','--','color','k')
%     ylabel('Y Position')
%     
%     subplot(3,1,2)
%     py=plot(times,y,'LineWidth',4); %hold on
%     
%     ay = gca;
%     ay.XLim = [-250,2000];
%     line([500,500],ay.YLim,'linestyle','--','color','k')
%     line([1000,1000],ay.YLim,'linestyle','--','color','k')
%     line([1500,1500],ay.YLim,'linestyle','--','color','k')
%     ylabel('Y Position')
%     
%     % plot differential
%     subplot(3,1,3)
%     pd = plot(times(1:length(times)-1),d,'LineWidth',2);
%     
%     ad = gca;
%     ad.XLim = [-250,2000];
%     line([500,500],ad.YLim,'linestyle','--','color','k')
%     line([1000,1000],ad.YLim,'linestyle','--','color','k')
%     line([1500,1500],ad.YLim,'linestyle','--','color','k')
%     ylabel('Euclidian Distance Differential')
%     
%     % add interative drawing
%     fixObj.fix1_start = drawpoint('color','g','label','F1','LabelVisible','Hover');
%     fixObj.fix1_end = drawpoint('color','r','label','F1','LabelVisible','Hover');
%     
%     fixObj.fix2_start = drawpoint('color','g','label','F2','LabelVisible','Hover');
%     fixObj.fix2_end = drawpoint('color','r','label','F2','LabelVisible','Hover');
%     
%     fixObj.fix3_start = drawpoint('color','g','label','F3','LabelVisible','Hover');
%     fixObj.fix3_end = drawpoint('color','r','label','F3','LabelVisible','Hover');
%     
%     
%     
%     % replot the x and y subplots with patches to reflet eye fixes
%     goodPlacement=0;
%     cnt=0;
%     while goodPlacement==0
%         
%         % create plots
%         subplot(3,1,1)
%         px=plot(times,x,'LineWidth',4); %hold on
%         
%         ax = gca;
%         ax.XLim = [-250,2000];
%         line([500,500],ax.YLim,'linestyle','--','color','k')
%         line([1000,1000],ax.YLim,'linestyle','--','color','k')
%         line([1500,1500],ax.YLim,'linestyle','--','color','k')
%         ylabel('X Position')
%         
%         % title
%         title(['Trial ' num2str(iTrial) '  U = update, C = confirm, B = bad trial'],'FontSize',36)
%         
%         % plot patch1
%         fix1_start = fixObj.fix1_start.Position(1);
%         fix1_end = fixObj.fix1_end.Position(1);
%         patchX=[fix1_start,fix1_end,fix1_end,fix1_start];
%         patchY = [ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
%         patch(patchX,patchY,'g','FaceAlpha',.2);
%         
%         % plot patch2
%         fix2_start = fixObj.fix2_start.Position(1);
%         fix2_end = fixObj.fix2_end.Position(1);
%         patchX=[fix2_start,fix2_end,fix2_end,fix2_start];
%         patchY = [ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
%         patch(patchX,patchY,'g','FaceAlpha',.2);
%         
%         % plot patch3
%         fix3_start = fixObj.fix3_start.Position(1);
%         fix3_end = fixObj.fix3_end.Position(1);
%         patchX=[fix3_start,fix3_end,fix3_end,fix3_start];
%         patchY = [ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
%         patch(patchX,patchY,'g','FaceAlpha',.2);
%         
%         % create plots
%         subplot(3,1,2)
%         py=plot(times,y,'LineWidth',4); %hold on
%         
%         ay = gca;
%         ay.XLim = [-250,2000];
%         line([500,500],ay.YLim,'linestyle','--','color','k')
%         line([1000,1000],ay.YLim,'linestyle','--','color','k')
%         line([1500,1500],ay.YLim,'linestyle','--','color','k')
%         ylabel('Y Position')
%         
%         % plot patch1
%         fix1_start = fixObj.fix1_start.Position(1);
%         fix1_end = fixObj.fix1_end.Position(1);
%         patchX=[fix1_start,fix1_end,fix1_end,fix1_start];
%         patchY = [ay.YLim(1),ay.YLim(1),ay.YLim(2),ay.YLim(2)];
%         patch_fix1 = patch(patchX,patchY,'g','FaceAlpha',.2);
%         
%         % plot patch2
%         fix2_start = fixObj.fix2_start.Position(1);
%         fix2_end = fixObj.fix2_end.Position(1);
%         patchX=[fix2_start,fix2_end,fix2_end,fix2_start];
%         patchY = [ay.YLim(1),ay.YLim(1),ay.YLim(2),ay.YLim(2)];
%         patch(patchX,patchY,'g','FaceAlpha',.2);
%         
%         % plot patch3
%         fix3_start = fixObj.fix3_start.Position(1);
%         fix3_end = fixObj.fix3_end.Position(1);
%         patchX=[fix3_start,fix3_end,fix3_end,fix3_start];
%         patchY = [ay.YLim(1),ay.YLim(1),ay.YLim(2),ay.YLim(2)];
%         patch(patchX,patchY,'g','FaceAlpha',.2);
%         
%         % get fig
%         fig=gcf;
%         
%         % check the RA is happy with the placement of the items
%         disp('PRESS "C" TO CONFIRM TRIAL.  PRESS "N" TO RE-DO.')
%         pressloop=1;
%         badTrial=0;
%         while pressloop
%             was_a_key = waitforbuttonpress;
%             if was_a_key && strcmp(get(fig, 'CurrentKey'), 'c')
%                 pressloop=0;
%                 goodPlacement=1;
%                 badTrial=0;
%             elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'u')
%                 pressloop=0;
%                 goodPlacement=0;
%             elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'b')
%                 pressloop=0;
%                 goodPlacement=1;
%                 badTrial=1;
%             end
%         end
%         
%     end
%     
%     fixStruct(iTrial).fix1_start = round(fix1_start);
%     fixStruct(iTrial).fix1_end = round(fix1_end);
%     fixStruct(iTrial).fix2_start = round(fix2_start);
%     fixStruct(iTrial).fix2_end = round(fix2_end);
%     fixStruct(iTrial).fix3_start = round(fix3_start);
%     fixStruct(iTrial).fix3_end = round(fix3_end);
%     fixStruct(iTrial).badTrial = badTrial;
%     
%     % clear a bunch of stuff
%     clear allClean allClean_idx allBlocked_idx postEM_idx preEM_idx
%     clear em1_idx fix1_start fix1_end fix2_start fix2_end fix3_start fix3_end patchX patchY
%     
%     % save figure
%     %saveas(fig,[plotDir '/' sprintf('sj%02d_cond%02d_trial%02d.jpeg',sjNum,iCond,iTrial)],'jpeg')
%     
%     % close figure
%     close
%     
% end
% 
% % save trial data [corrected, was wrong for sj31cond03
% trialData = eyeStruct(iCond).trialData;
% 
% 
% % % save data
% % save([destDir '/' sprintf('sj%02d_cond%02d.mat',sjNum,iCond)],'fixStruct','trialData','times')
% % 
% 
