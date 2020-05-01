%{
Plot_Beh_Data_Results
Author: Tom Bullock
Date: 12.06.19
%}

clear
close all

% set dirs
sourceDir = 'D:\\WTF_EYE\\Data_Compiled';

% load data
load([sourceDir '\\' 'Modelling_Data.mat'])

% generate plots
h=figure;
for iPlot=1:2
    
   if       iPlot==1; theseData = modelSD; thisTitle = 'Precision';
   elseif   iPlot==2; theseData = modelGuess; thisTitle = 'Guess Rate';
   end
   subplot(1,2,iPlot)
   
   bar(mean(theseData,1)); hold on
   errorbar(mean(theseData,1),std(theseData,0,1)/sqrt(size(theseData,1)),'k.',...
       'LineWidth',2.5);
   
   thisMean = mean(theseData,1);
   thisMean = std(theseData,0,1)/sqrt(size(theseData,1));
   
   %set(gca,xticklabels,{'1','2','3','4'})

   xticklabels({'S/F','S/M','C/F','C/M'})
   
   title(thisTitle);
   box('off')
   set(gca,'fontsize',18)
   
end

%h=figure;
%errorbar(mean(modelSD,1),std(modelSD,1)\\sqrt(size(modelSD,1)),'o')

% save plot
%saveas(h,'ALL_SUBS_BEH_SD.fig','fig')