%{
Plot_IEM_Within_Surf_Fixed
Author: Tom Bullock
Date: 10.23.20

Plot fixed model CRFs as surf plots.

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Within_Fixed'; 

% select subjects
subjects = [1:7,9,14,16:20,22:27,31];

for iSub=1:length(subjects)
    sjNum = subjects(iSub);
    
   load([sourceDir '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)])
   
   allData(iSub,:,:,:) = squeeze(mean(em.allTF,3)); % avg across iterations, only 1 freq (alpha) currently (allTF = cond x freq x iter x sample x chan)
   
end

h=figure('units','normalized','outerposition',[0 0 1 1]);

for iPlot=1:4
    
    subplot(2,2,iPlot)
    
    surf(squeeze(mean(allData(:,iPlot,:,:),1)))
    
    set(gca,...
        'ytick',linspace(1,640,6),...
        'yticklabel',-500:500:2000);
    
    title(['Cond' num2str(iPlot)])
    
end
