%{
Eye_Movement_Analysis_Plot
Author: Tom Bullock
Date: 06.05.20

plot euclidian errors for conds 3 and 4 for ems 1 and 2

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Euclidian_Error_mats';

% select subjects
subjects = [2,31];


for iSub=1:length(subjects)
    
    sjNum = subjects(iSub);
    
    % load data
    load([sourceDir '/' sprintf('sj%02d_euclid_error.mat',sjNum)])
    
    % get means for each em and each cond
    cond3(iSub,1) = nanmean([euclidErrorStruct(3).em1]);
    cond3(iSub,2) = nanmean([euclidErrorStruct(3).em2]);
    
    cond4(iSub,1) = nanmean([euclidErrorStruct(4).em1]);
    cond4(iSub,2) = nanmean([euclidErrorStruct(4).em2]);
    
end

% generate plots
for iPlot=1:2
    
    if iPlot==1; theseData = cond3; thisTitle = 'Spatial Move';
    elseif iPlot==2; theseData = cond4; thisTitle = 'Color Move';
    end
    
    thisMean = mean(theseData);
    thisSEM = std(theseData,0,1)./sqrt(size(theseData,1));
    
    subplot(1,2,iPlot);
    bar(thisMean,'b'); hold on
    
    errorbar(thisMean,thisSEM,...
        'LineStyle','none',...
        'LineWidth',2,...
        'color','k')
    
    set(gca,...
        'box','off',...
        'linewidth',1.5,...
        'xlim',[.5,2.5],...
        'xTickLabel',{'EM1','EM2'})
    
    title(thisTitle);
        
    
end





