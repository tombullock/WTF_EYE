%{
Self_Report_Stats_Resample

Author: Tom Bullock, UCSB Attention Lab
Date: 12.21.19

%}

clear
close all

%% set dirs
sourceDir = '/home/bullock/BOSS/CPT_Adaptation/Data_Compiled';

%% choose subjects
[~,mySubjects] = CPT_SUBJECTS;

%% remove sj102
mySubjects(1)=[];

%% load self-report data
load([sourceDir '/' 'Self_Report.mat'])

%% isolate subjects
allPain = allPain(mySubjects-100,:,:)

%% null loop (compare Tx and Ct)
for i=1:1000 % iteration loop
    
    for k=1:5 % exposure loop
        
        observedData = squeeze(allPain(:,:,k));
        
        for m=1:length(observedData)
            thisPerm = randperm(size(observedData,2));
            nullDataMat(m,:) = observedData(m,thisPerm);
        end
        
        [H,P,CI,STATS] = ttest(nullDataMat(:,1),nullDataMat(:,2));
        tValsNull(i,k) = STATS.tstat;
        clear STATS nullDataMat observedData
        
    end
end

%% observed data loop
for k=1:5 % exposure loop
    
    observedData = squeeze(allPain(:,:,k));   
    [H,P,CI,STATS] = ttest(observedData(:,1),observedData(:,2));
    tValsObs(k) = STATS.tstat;
    clear STATS nullDataMat observedData
    
end


%% compare observed t-test at each trial and timepoint with corresponding null distribution of t-values

% trial (exposure) loop
for k=1:size(tValsNull,2)
    
    thisNullDist = sort(tValsNull(:,k),1,'descend');
    [~,idx] = min(abs(thisNullDist - tValsObs(k)));
    tStatIdx(k) = idx;
    clear thisNullDist  idx
    
end


% convert idx value to probability
%pValuesPairwise = tStatIdx./1000;

% convert idx value to a vector of 0's and 1's (representing significance
% at p<.05)
clear sigVec
for k=1:size(tStatIdx,2)
    if tStatIdx(k)<25 || tStatIdx(k)>975
        sigVec(k) = 1;
    else
        sigVec(k)=0;
    end
end

%% quickly get ANOVA results
addpath(genpath('/home/bullock/BOSS/CPT_Adaptation/resampling'))
% name variables
var1_name = 'tx/ct';
var1_levels = 2;
var2_name = 'trial';
var2_levels = 5;

observedData = [squeeze(allPain(:,1,:)),squeeze(allPain(:,2,:))];

statOutput = teg_repeated_measures_ANOVA(observedData,[var1_levels var2_levels],{var1_name, var2_name});

%% get descriptive stats
descriptives.condMeans = nanmean(allPain);
descriptives.condSEMs = nanstd(allPain,0,1)./sqrt(size(allPain,1));


% save data
save([sourceDir '/' 'STATS_Resampled_Self_Report.mat'],'sigVec','tValsObs','tValsNull','mySubjects','statOutput','descriptives')














% 
% % iteration loop (null data)
% for i=1:1000
%     
%     disp(['Null Iteration ' num2str(i)])
%     
%     % sample (timepoint) loop
%     for j=1:size(paMatAll,5)
%         
%         % trial (CPT/WPT exposure) loop
%         for k=1:size(paMatAll,3)
%             
%             observedData = squeeze(mean(paMatAll(:,:,k,:,j),4)); % subs x conds (average across R/L eyes)
%             
%             for m=1:length(observedData)              
%                 thisPerm = randperm(size(observedData,2));              
%                 nullDataMat(m,:) = observedData(m,thisPerm);       
%             end
%             
%             [H,P,CI,STATS] = ttest(nullDataMat(:,1),nullDataMat(:,2));
%             tValsNull(i,j,k) = STATS.tstat;
%             clear STATS nullDataMat observedData
%             
%         end
%     end
% end
% 
% %% generate matrix of real t-test results
% 
% % sample (timepoint) loop
% for j=1:size(paMatAll,5)
%     % trial (CPT/WPT loop) 
%     for k=1:size(paMatAll,3)
%         
%         observedData = squeeze(mean(paMatAll(:,:,k,:,j),4)); % subs x conds (average across R/L eyes)
%         [H,P,CI,STATS] = ttest(observedData(:,1),observedData(:,2));
%         tValsObs(j,k) = STATS.tstat;
%         clear STATS observedData 
%  
%     end
% end
% 
%     
% %% compare observed t-test at each trial and timepoint with corresponding null distribution of t-values
% 
% % sample (timepoint) loop
% for j=1:size(tValsNull,2)
%     
%     % trial (exposure) loop
%     for k=1:size(tValsNull,3)
% 
%         thisNullDist = sort(tValsNull(:,j,k),1,'descend');
%         [~,idx] = min(abs(thisNullDist - tValsObs(j,k)));
%         tStatIdx(j,k) = idx;
%         clear thisNullDist  idx
%         
%     end
% end
% 
% % convert idx value to probability
% %pValuesPairwise = tStatIdx./1000;
% 
% % convert idx value to a vector of 0's and 1's (representing significance
% % at p<.05)
% clear sigVec
% for j=1:size(tStatIdx,1)
%     for k=1:size(tStatIdx,2)
%         if tStatIdx(j,k)<25 || tStatIdx(j,k)>975
%             sigVec(j,k) = 1;
%         else
%             sigVec(j,k)=0;
%         end
%     end
% end
% 
% % save data
% save([sourceDir '/' 'STATS_Resampled_Self_Report.mat'],'sigVec','tValsObs','tValsNull','subjects','badSubs')



