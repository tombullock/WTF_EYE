%{
SCAM Plot Slopes for Within IEM
Author: Tom Bullock
Date: 04.15.20

%}

clear 
close all

% set dirs
sourceDir = '/Users/tombullock/Documents/Psychology/WTF_EYE/Data_Compiled';
destDir = '/Users/tombullock/Documents/Psychology/WTF_EYE/Plots';

% load data
load([sourceDir '/' 'IEM_Slopes_Within.mat']);

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
xNewTickLabel = ['-500';'-250';'  0 ';' 250';' 500';' 750';'1000';'1250';'1500';'1750';'2000'];
thisYlim = [-.001 .005];%'auto';

% plot total slope
h=figure('Units','Normalized','OuterPosition',[0,0,.9,.65]); % change these numbers to change plot size
for iCond=1:4
    
    if iCond==1; thisLineStyle='r'; % S/F
    elseif iCond==2; thisLineStyle='b'; % C/F
    elseif iCond==3; thisLineStyle='g'; % S/M
    elseif iCond==4; thisLineStyle='m'; % C/M
    end
    
    % plot real data
    realLineIdx(iCond)=plot(squeeze(mean(rSl_total(:,iCond,f,:),1)),thisLineStyle,'LineWidth',observedLineWidth); hold on

    if iCond==1; thisLineStyle='r--'; % S/F
    elseif iCond==2; thisLineStyle='b--'; % C/F
    elseif iCond==3; thisLineStyle='g--'; % S/M
    elseif iCond==4; thisLineStyle='m--'; % C/M
    end
    
    % plot permuted data
    plot(squeeze(mean(pSl_total(:,iCond,f,:),1)),thisLineStyle,'LineWidth',observedLineWidth); hold on

    % plot settings
    box('off')
    set(gca,'xTick',xNewTick,'xticklabel',xNewTickLabel,'xLim',[1 640],'yLim',thisYlim,'fontsize',axisFontSize,'LineWidth',axisLineWidth);
    vline(thisSR*.5,'g--')
    vline(thisSR*.75,'r--')
    vline(thisSR*1.25,'b--')
    xlabel('Time (ms)','FontSize',xLabelFontSize)
    ylabel('Slope','FontSize',yLabelFontSize)
    box('off')
    pbaspect([3 1 1])
    title('SCAM - CTF Slopes','FontSize',titleFontSize)
end

% set legend
h_legend = legend(realLineIdx,'Spatial/Fix','Color/Fix','Spatial/Move','Color/Move');
set(h_legend,'FontSize',20);

% save figure
saveas(h,[destDir '/' 'CTF_Slope_Alpha_Within.jpeg'],'jpeg')