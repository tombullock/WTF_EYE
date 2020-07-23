%{
Eye_Movement_Analysis_Plot
Author: Tom Bullock
Date: 06.05.20

plot euclidian errors for conds 3 and 4 for ems 1 and 2

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Euclidian_Error_mats';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% select subjects (multiple)
subjects = [1:7,9,10,11,12,16,17,19,20,22:27,31];

h=figure('units','normalized','OuterPosition',[0.5703 0.1100 0.3347 0.8150]);
for iSub=1:length(subjects)
    
    sjNum = subjects(iSub);
    
    % load data
    load([sourceDir '/' sprintf('sj%02d_euclid_error.mat',sjNum)])
    
    % get means for each em and each cond
    cond3(iSub,1) = nanmean([euclidErrorStruct(3).em1]);
    cond3(iSub,2) = nanmean([euclidErrorStruct(3).em2]);
    
    cond4(iSub,1) = nanmean([euclidErrorStruct(4).em1]);
    cond4(iSub,2) = nanmean([euclidErrorStruct(4).em2]);
    
    % get times
    cond3_idx(iSub,1) = nanmean([euclidErrorStruct(3).em1_idx]);
    cond3_idx(iSub,2) = nanmean([euclidErrorStruct(3).em2_idx]);
    
    cond4_idx(iSub,1) = nanmean([euclidErrorStruct(4).em1_idx]);
    cond4_idx(iSub,2) = nanmean([euclidErrorStruct(4).em2_idx]);
    
end

cond3_visAngle = visAngleCalculate_pix(cond3,120);
cond4_visAngle = visAngleCalculate_pix(cond4,120);






% convert lat index to wide format for R
for iData=1:2
    
    % get data 
    if iData==1
        theseData = [cond3_visAngle,cond4_visAngle]; % euclid error
    else
        theseData = [cond3_idx,cond4_idx]; % euclid error time
    end
    
    nSubs = size(theseData,1);

    % convert into R "long" format
    dumMem = [zeros(nSubs,1), ones(nSubs,1), zeros(nSubs,1), ones(nSubs,1)];
    dumEyes = [zeros(nSubs,1), zeros(nSubs,1), ones(nSubs,1), ones(nSubs,1)];
    
    % create column vector of subjects
    sjNums = 1:nSubs;
    sjNums = [sjNums, sjNums, sjNums, sjNums]';
    
    % use colon operator function to get into long format
    if iData==1
        euclid_error_LF = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
    elseif iData==2
        euclid_time_LF = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
    end
    
end


save([destDir '/' 'EYE_Euclid_Error.mat'],'cond3_visAngle','cond4_visAngle','cond3_idx','cond4_idx','euclid_error_LF','euclid_time_LF','-v7.3')





% generate plots for dist
subplot(1,2,1);
for iPlot=1:2
    
    if iPlot==1; theseData = cond3; thisTitle = 'Spatial Move'; thisColor = [0,204,0];thisX=1:2;
    elseif iPlot==2; theseData = cond4; thisTitle = 'Color Move'; thisColor = [204,0,204]; thisX=3:4;
    end
    
    % compute mean and SEM in pixels
    thisMean = mean(theseData);
    thisSEM = std(theseData,0,1)./sqrt(size(theseData,1));
    
    % convert to visual angles
    thisMean = visAngleCalculate_pix(thisMean,120);
    thisSEM = visAngleCalculate_pix(thisSEM,120);
    
    %subplot(1,2,iPlot);
    bar(thisX,thisMean,'FaceColor',thisColor./255); hold on
    
    errorbar(thisX,thisMean,thisSEM,...
        'LineStyle','none',...
        'LineWidth',2,...
        'color','k',...
        'CapSize',0)
    
    set(gca,...
        'box','off',...
        'linewidth',1.5,...
        'xlim',[.5,4.5],...
        'fontsize',22,...
        'xTick',[])
    
    pbaspect([1,1,1]);
    
    %title(thisTitle);
        
    
end


% generate plots for time (sample)
subplot(1,2,2);
for iPlot=1:2
    
    if iPlot==1; theseData = cond3_idx; thisTitle = 'Spatial Move'; thisColor = [0,204,0];thisX=1:2;
    elseif iPlot==2; theseData = cond4_idx; thisTitle = 'Color Move'; thisColor = [204,0,204]; thisX=3:4;
    end
    
    % compute mean and SEM in pixels
    thisMean = mean(theseData);
    thisSEM = std(theseData,0,1)./sqrt(size(theseData,1));
    
    % convert from samples (500 Hz) to milliseconds
    thisMean = thisMean*2;
    
%     % convert to visual angles
%     thisMean = visAngleCalculate_pix(thisMean,120);
%     thisSEM = visAngleCalculate_pix(thisSEM,120);
    
    %subplot(1,2,iPlot);
    bar(thisX,thisMean,'FaceColor',thisColor./255); hold on
    
    errorbar(thisX,thisMean,thisSEM,...
        'LineStyle','none',...
        'LineWidth',2,...
        'color','k',...
        'CapSize',0)
    
    set(gca,...
        'box','off',...
        'linewidth',1.5,...
        'xlim',[.5,4.5],...
        'ylim',[200,350],...
        'fontsize',22,...
        'xTick',[])
    
    pbaspect([1,1,1]);
    
    %title(thisTitle);
        
    
end





