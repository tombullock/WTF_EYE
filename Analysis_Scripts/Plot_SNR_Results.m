%{
Plot_SNR_Results
Author: Tom Bullock
Date: 05.30.20

%}

clear
close all


sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

load([sourceDir '/' 'SNR_Data.mat'])

% baseline correct
snrMat = snrMat - mean(snrMat(:,:,:,77:128),4);

% plot
for iCond=1:4
    
    if iCond==1; thisLineStyle='-'; thisColor = 'r'; % S/F
    elseif iCond==2; thisLineStyle='-'; thisColor = 'b'; % C/F
    elseif iCond==3; thisLineStyle='-'; thisColor = 'g'; % S/M
    elseif iCond==4; thisLineStyle='-'; thisColor = 'm'; % C/M
    end
    
%     plot(linspace(-200,2000,564),squeeze(mean(mean(snrMat(:,iCond,:,77:640),1),3)),...
%         'color',thisColor,...
%         'linewidth',3); hold on
    
    theseDataMean = squeeze(mean(mean(snrMat(:,iCond,:,77:640),1),3));
    theseDataSEM = squeeze(std(mean(snrMat(:,iCond,:,77:640),3),0,1))./sqrt(size(snrMat,1));
    
    shadedErrorBar(linspace(-200,2000,564),theseDataMean,theseDataSEM,...
        {'color',thisColor,...
        'linewidth',3}); hold on
        
    
end

set(gca,...
    'xlim',[-200,2000],...
    'linewidth',1.5,...
    'fontsize',24);

ylabel('SNR (BL Corrected)')
xlabel('Time (ms')
title('SNR - Total Alpha vs. All Non-Alpha Freqs')

%legend('S/F','C/F','S/M','C/M','location','northwest')