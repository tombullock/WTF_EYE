%==========================================================================
%{
Beh_WM_Modelling
Purpose: model WM data to get precision and resoution.  Script also 
 generates an "error" vector for each subject, which should go into a
 "Beh_errors" folder (this is useful for other analyses e.g. correlating
 error magnitude with EEG measures)

Author: Tom Bullock, UCSB Attention Lab

Date last modified: 10.30.19

To-Do: add double session compatibility

%}
%==========================================================================

clear 
close all

% add dependencies
%addpath(genpath('\\home\\bullock\\WTF\\Behavior\\visionlab-MemToolbox-fea8609'))
addpath(genpath('D:\\WTF_EYE\\Dependancies\\visionlab-MemToolbox-fea8609'))

% set dirs
rDir='/home/waldrop/Desktop/WTF_EYE';
%sourceDir = '\\home\\bullock\\WTF\\Behavior\\Beh_Merged';
sourceDir = [rDir '/' 'Beh_Data_Processed'];
%destDirCompiledData = '\\home\\bullock\\WTF\\Data_Compiled';
destDirCompiledData = [rDir '/' 'Data_Compiled'];
%destDirErrorsData = '\\home\\bullock\\WTF\\Behavior\\Beh_Errors';
destDirErrorsData = [rDir '/' 'Beh_Errors'];

%set verbosity (0=just get values, 2=get plots for each individual\\cond)
setVerbosity = 0;

% set subject numbers
subjectNumbers = [1:7,9:14,16:20,22:27,31];

%subjectNumbers = 5; % somthing up with sj05

%matlabpool open 72  % this just boosts MEMFIT speed (not running subs in parallel as per usual)
for subjLoop=1:size(subjectNumbers,2)
        
   sjNum= subjectNumbers(subjLoop);
   
   load([sourceDir '/' sprintf('sj%02d_allBeh.mat',sjNum)])
   
   % merge  data across 2 sessions
   for i=cbOrder
       masterStructMerged(i).allTrialData = [masterStruct(1,i).allTrialData, masterStruct(2,i).allTrialData];
   end
   masterStruct = [];
   masterStruct = masterStructMerged;

   clear allFit
   
   for condLoop=1:4
       
       disp(['Processing sj ' num2str(sjNum) ' condition ' num2str(condLoop)])
       
       clear fit x thisMouseAngleDegs clear 
       
       % remove broken trials
       cnt=0;
       for i=1:length(masterStruct(condLoop).allTrialData)
           if masterStruct(condLoop).allTrialData(i).brokenTrial==0
           cnt=cnt+1;
           unbrokenStruct(cnt)=masterStruct(condLoop).allTrialData(i);
           end
       end
       
       % isolate color distances
       if ismember(condLoop,[2,4])
           for i=1:length(unbrokenStruct)
               x(i,1)=unbrokenStruct(i).col_responsedistance;
           end
       end
       
       % isolate spatial distances
       if ismember(condLoop,[1,3])
           
           for i=1:length(unbrokenStruct)            
               thisMouseAngleDegs = 360-(180 - 180*(unbrokenStruct(i).thisMouseAngle/pi)); % from WTF script             
               x(i,1)= round(unbrokenStruct(i).stimLocAngle - thisMouseAngleDegs);              
           end
               
           for i = 1:length(x)
               if x(i,:) > 180
                   x(i,:) = x(i,:)-360;
               elseif x(i,:) < -180
                   x(i,:) = x(i,:)+360;
               end
           end

       end
       
      % fit mixture model with bias
       model=WithBias(StandardMixtureModel);
       fit = MemFit(x,model,'Verbosity',setVerbosity);
       mu=fit.maxPosterior(1);
       g=fit.maxPosterior(2);
       sd=fit.maxPosterior(3);
       
       modelMu(subjLoop,condLoop) = mu;
       modelGuess(subjLoop,condLoop) = g;
       modelSD(subjLoop,condLoop) = sd;
       
       %save error vector for other analyses
       save([destDirErrorsData '\\' sprintf('sj%d%02d_beh_errors.mat',sjNum,condLoop)],'x','model','fit')
   end
   
end
 
% save summary fit data stats
save([destDirCompiledData '\\' 'Modelling_Data.mat'],'modelGuess','modelSD','modelMu');








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old but potentially useful stuff...



% 
% 
%    for condLoop= 1:4 
%        
%        data = [];
% 
%        load([sourceDir '\\' sprintf('sj%d%02d_newBeh.mat',sjNum,condLoop)]);
% 
%        disp(['Processing Subject ' num2str(sjNum) '0' num2str(condLoop)])
% 
%        %loop through trialInfoUnbroken structure and extract two columns of
%        %data (actual location in degs and response location in degs)
%        for j = 1:length(trialInfoUnbroken)    
%           data(j,:) = [trialInfoUnbroken(j).stimLocAngle, trialInfoUnbroken(j).thisMouseAngleDegs];     
%        end
% 
%        x = round(data(:,2) - data(:,1));
% 
%        for i = 1:length(x)
%            if x(i,:) > 180
%                x(i,:) = x(i,:)-360;
%            elseif x(i,:) < -180
%                x(i,:) = x(i,:)+360;
%            end
%        end
%        
%        % do model with bias
%        model=WithBias(StandardMixtureModel);
%        
%        fit = [];
%        fit = MemFit(x,model,'Verbosity',setVerbosity);
%        
% % %        % if fitting standard mixture model       
% % %        g = fit.maxPosterior(1,1);
% % %        sd = fit.maxPosterior(1,2);
%        
%        % if fitting mixture model with bias
%        mu=fit.maxPosterior(1);
%        g=fit.maxPosterior(2);
%        sd=fit.maxPosterior(3);
%        
%        modelMu(subjLoop,condLoop) = mu;
%        modelGuess(subjLoop,condLoop) = g;
%        modelSD(subjLoop,condLoop) = sd;
%        
%        %save error vector for other analyses
%        save([destDirErrorsData '\\' sprintf('sj%d%02d_beh_errors.mat',sjNum,condLoop)],'x','model','fit')
%       
% %        % to plot the fit at later stage after running all subjects we can
% %        % use the following plotting functions (from MemFit.m, line 200 ish)
% %        % model
% %        % fit
% %        % data.errors = x
% %        data.errors = x;
% %        
% %        posteriorSamples = fit.posteriorSamples;
% %        
% %        % plots model fit
% %        PlotModelFitInteractive(model, fit.maxPosterior, data);
% %        
% %        % posterior etc.
% %        h = PlotModelParametersAndData(model, posteriorSamples, data);
% %        
% %        % Posterior predictive plot (i.e. how good the model fit is)
% %        h = PlotPosteriorPredictiveData(model, posteriorSamples, data);
%       
%    end
%    
%    %saveas(sprintf('sj%d_memFit.fig',sjNum),'fig')
%    
% end
% 
% save([destDirCompiledData '\\' 'Modelling_Data.mat'],'modelGuess','modelSD');
% 
% % create simple plot
% h=figure;
% errorbar(mean(modelSD,1),std(modelSD,1)\\sqrt(size(modelSD,1)),'o')
% 
% % save plot
% %saveas(h,'ALL_SUBS_BEH_SD.fig','fig')
% 
% matlabpool close