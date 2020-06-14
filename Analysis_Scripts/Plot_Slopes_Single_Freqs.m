%{
SCAM Plot Slopes for Within IEM
Author: Tom Bullock
Date: 04.22.20

Updated to plot BF thresholded bars to indicate sig
%}

clear 
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
destDir = '/home/waldrop/Desktop/WTF_EYE/Plots';
sourceDirBayesResults = '/home/waldrop/Desktop/WTF_EYE/R_Results';

% load data
load([sourceDir '/' 'IEM_Slopes_Within.mat']);

% load BF results
load([sourceDirBayesResults '/' 'BF_Results_IEM_WIthin_Alpha_RvP.mat'])

% set freq band (1=alpha,2=theta)
f=1;

%% set up plot
% line widths
observedLineWidth = 4;
nullLineWidth = 2;
axisLineWidth = 1.5;

% font sizes
axisFontSize = 18;
titleFontSize = 22;
xLabelFontSize = 24;
yLabelFontSize = 24;

% colors
critTcolor = [.5 .5 .5];
obsTcolor = [0.04 0.52 .78];

% axis limits and labels
thisSR = 256; % sample rate
xNewTick = [1, thisSR*.25, thisSR*.5, thisSR*.75, thisSR*1, thisSR*1.25, thisSR*1.5, thisSR*1.75, thisSR*2, thisSR*2.25, thisSR*2.5];
xNewTickLabel = ['-0.5';'-.25';'  0 ';' .25';' 0.5';' .75';'  1 ';'1.25';' 1.5';'1.75';'  2 '];
thisYlim = [-.0017 .004];%'auto';

% plot total slope
thisYpos = -.0006; % start position for stats plots
h=figure('Units','Normalized','OuterPosition',[0,0,1,.75]); % change these numbers to change plot size
for iCond=1:4
    
    if iCond==1; thisLineStyle='-'; thisColor = 'r'; % S/F
    elseif iCond==2; thisLineStyle='-'; thisColor = 'b'; % C/F
    elseif iCond==3; thisLineStyle='-'; thisColor = 'g'; % S/M
    elseif iCond==4; thisLineStyle='-'; thisColor = 'm'; % C/M
    end
    
    % plot real data
    realLineIdx(iCond)=plot(squeeze(mean(rSl_total(:,iCond,f,:),1)),'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor); hold on

    if iCond==1; thisLineStyle='--'; thisColor = 'r'; % S/F
    elseif iCond==2; thisLineStyle='--'; thisColor = 'b'; % C/F
    elseif iCond==3; thisLineStyle='--'; thisColor = 'g'; % S/M
    elseif iCond==4; thisLineStyle='--'; thisColor = 'm'; % C/M
    end
    
    % plot permuted data
    plot(squeeze(mean(pSl_total(:,iCond,f,:),1)),'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor); hold on

    % plot stats
    bfTreshold = 3;
    thisYpos = thisYpos-.0002;
    for b=1:size(bfGrandOutput,1)
        if bfGrandOutput(b,iCond)>bfTreshold
            line([b,b+1],[thisYpos,thisYpos],...
                'LineWidth',4,...
                'LineStyle',thisLineStyle,...
                'Color',thisColor)
        end
    end
  
    
    
    % plot settings
    box('off')
    set(gca,'xTick',xNewTick,'xticklabel',xNewTickLabel,'xLim',[1 640],'yLim',thisYlim,'fontsize',axisFontSize,'LineWidth',axisLineWidth);
    vline(thisSR*.5,'k--')
    vline(thisSR*.75,'k--')
    vline(thisSR*1,'k--')
    vline(thisSR*1.5,'k--')
    vline(thisSR*2,'k--')
    xlabel('Time (s)','FontSize',xLabelFontSize)
    ylabel('Slope','FontSize',yLabelFontSize)
    box('off')
    pbaspect([3 1 1])
    %title('CRF Slopes','FontSize',titleFontSize)
end

% set legend
h_legend = legend(realLineIdx,'Spatial/Fix','Color/Fix','Spatial/Move','Color/Move');
set(h_legend,'FontSize',16);

% save figure
saveas(h,[destDir '/' 'CTF_Slope_Alpha_Within.jpeg'],'jpeg')
saveas(h,[destDir '/' 'CTF_Slope_Alpha_Within.eps'],'epsc')