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

% select subjects (multiple)
subjects = [1:7,9,10,11,12,16,17,19,20,22:27,31];

for iSub=1:length(subjects)
    
    sjNum=subjects(iSub)
    
    % select subject number
    %sjNum=23;
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
                this_em1_ed_time = min_em1_fixLoc1_idx;
            else
                firstFixLoc = 2;
                fixLoc1_color = 'b';
                this_em1_ed = min_em1_fixLoc2;
                this_em1_pix = [ex1(min_em1_fixLoc2_idx),ey1(min_em1_fixLoc2_idx)];
                this_em1_ed_time = min_em1_fixLoc2_idx;
            end
            
            if firstFixLoc == 1
                fixLoc2_color = 'b';
                this_em2_ed = min_em2_fixLoc2;
                this_em2_pix = [ex2(min_em2_fixLoc2_idx),ey2(min_em2_fixLoc2_idx)];
                this_em2_ed_time = min_em2_fixLoc2_idx;
            else
                fixLoc2_color = 'g';
                this_em2_ed = min_em2_fixLoc1;
                this_em2_pix = [ex2(min_em2_fixLoc1_idx),ey2(min_em2_fixLoc1_idx)];
                this_em2_ed_time = min_em2_fixLoc1_idx;
            end
            
            % save euclidian distance errors to a struct
            euclidErrorStruct(iCond).em1(iTrial) = this_em1_ed;
            euclidErrorStruct(iCond).em2(iTrial) = this_em2_ed;
            
            euclidErrorStruct(iCond).em1_idx(iTrial) = this_em1_ed_time;
            euclidErrorStruct(iCond).em2_idx(iTrial) = this_em2_ed_time;
            
            
            
            
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
    
    % save data
    save([destDir '/' sprintf('sj%02d_euclid_error.mat',sjNum)],'euclidErrorStruct')
    
    clear centerX centerY euclidErrorStruct eyeProcessed eyeStruct iTrial nTrials theseTimes theseTimesIdx times this_em1_ed_time this_em2_ed_time
    
end
