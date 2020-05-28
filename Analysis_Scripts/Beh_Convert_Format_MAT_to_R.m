%==========================================================================
%{
Beh_Convert_Format_MAT_to_R
Purpose: Convert behavioral data from MATLAB (wide) to R (long)format
ready for statistcal analysis
Author: Tom Bullock, UCSB Attention Lab
Date last modified: 03.04.18
%}
%==========================================================================

clear all
close all

% define behavioral source
sourceFolder = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% load behavioral data file (contains SD and Guess rate)
load([sourceFolder '/' 'Modelling_Data.mat'])

%convert data to "long" format for analysis in R (72 x 4)
for i=1:2

    thisData = [];
    
    if i==1
        thisData = modelSD;
    else
        thisData = modelGuess;
    end
    
    % convert into R "long" format
    dumMem = [zeros(25,1), ones(25,1), zeros(25,1), ones(25,1)];
    dumEyes = [zeros(25,1), zeros(25,1), ones(25,1), ones(25,1)];
    
    % create column vector of subjects
    sjNums = 1:size(thisData,1);
    sjNums = [sjNums, sjNums, sjNums, sjNums]';
    
    % use colon operator function to get into long format
    if i==1
        modelSD_LF = [thisData(:), sjNums, dumMem(:), dumEyes(:)];
    else
        modelGuess_LF = [thisData(:), sjNums, dumMem(:), dumEyes(:)];
    end
    
end

save([sourceFolder '/' 'Modelling_Data_LF.mat'],'modelSD_LF','modelGuess_LF','-v7.3')

clear all
close all