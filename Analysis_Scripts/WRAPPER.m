%{
WRAPPER
Author: Tom Bullock, UCSB Attention Lab
Date: 11.18.19

Wrapper function for WTF_EYE study

Notes:

Place raw behavior data in Beh_Data dir
Place raw EEG data in EEG_Raw dir
Make sure all subfolders exist
%}

%% which subs to process?
subjects = [];
%subjects = [7,9];


%% Trial Data 
%Beh_Merge_Blocks(subjects)

%% EEG
%EEG_Preprocessing1(subjects)
EEG_Preprocessing2(subjects)

IEM(subjects)
IEM_Compute_Compile_Slopes_Within(subjects)
IEM_Cross_TT(subjects) 
IEM_Cross_TT_Compile_Slopes


