%{
EEG_Compute_ERPs
Author: Tom Bullock
Date:05.28.20

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/EEG_Processed';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
plotDir = '/home/waldrop/Desktop/WTF_EYE/Plots';

addpath(genpath('/home/waldrop/Desktop/WTF_EYE/Analysis_Scripts'))

% select subs
subjects = [1:7,9:14,16:20,22:27,31];

% select elects
P1_elects = [25,26,27,29,30,62:64];

% select times (-200 to 500 ms locked around stim onset)
timeIdx = 78:257;

% sub loop
for iSub=1:length(subjects)

    sjNum=subjects(iSub)
    
    load([sourceDir '/' sprintf('sj%02d_eeg.mat',sjNum)])
    
    % condition loop
    for iCond=1:4
        
        avgERP(iSub,iCond,:) =  squeeze(mean(mean(allConds(iCond).eegs(P1_elects,timeIdx,:),1),3));
    
    end
    
end

% plot 

% baseline correct
avgERP_bl = avgERP - mean(avgERP(:,:,1:51),3);

for iCond=1:4
    
    if       iCond==1; thisColor = 'r';
    elseif   iCond==2; thisColor = 'b';
    elseif   iCond==3; thisColor = 'g';
    elseif   iCond==4; thisColor = 'm';
    end
    
    thisMean = squeeze(mean(avgERP_bl(:,iCond,:),1));
    thisSEM= squeeze(std(avgERP_bl(:,iCond,:),0,1)./sqrt(size(avgERP_bl,1)));
    
%     shadedErrorBar(linspace(-200,500,180),thisMean,thisSEM,...
%         {'color',thisColor,...
%         'LineWidth',3}); hold on
    
    plot(linspace(-200,500,180),thisMean,...
        'color',thisColor,...
        'LineWidth',3); hold on
        
    set(gca,...
        'ylim',[-2,4],...
        'linewidth',1.5,...
        'box','off',...
        'fontsize',24)
    
end

line([0,0],[-2,4],'linewidth',2,'linestyle','--','color','k')

title('Grand Average ERPs (PO/O Electrodes)')
xlabel('Time (ms)')
ylabel('Amplitude (uV)')

legend('Spatial/Fix','Color/Fix','Spatial/Move','Color/Move','fontsize',24,'location','northwest')





