%{
WRAPPER
Author: Tom Bullock, UCSB Attention Lab
Date: 11.18.19

Wrapper function for WTF_EYE study

Notes:

Place raw behavior data in Beh_Data
Place raw EEG data in EEG_Raw

%}

% ADD SJ NUMS HERE AND INCLUDE INPUS TO FUNCTIONS (MAKE EASY TO RUN IN ONE
% GO!)

%% Trial Data 
Beh_Rename_Files
Beh_Merge_Blocks

%% EEG
EEG_Preprocessing1
EEG_Preprocessing2
IEM
IEM_Cross_TT

