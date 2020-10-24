# Tracking Disruptions of Spatial Representations in Working Memory by Eye-Movements (SCAM!)

## Abstract: 

Oscillations in the alpha frequency band (~8-12 Hz) play an important role in supporting selective attention to visual items and maintenance of their spatial locations in human working memory (WM).  Recent findings suggest that the maintenance of spatial information in alpha is modulated by different forms of interruption to continuous visual input, such that attention shifts, eye-closure and backward masking of the encoded item cause reconstructed representations of remembered locations to become degraded.  Here, we investigated how another common visual disruption - eye-movement - modulates reconstructions of stored representations of items held in WM.  Participants completed trials of a delayed estimation task, where they encoded and recalled either the location (spatial memory) or color (non-spatial memory) of an object after a brief retention period.  During retention, participants either fixated at center or executed a sequence of guided eye-movements.  Electroencephalography (EEG) was recorded at the scalp and gaze was monitored using gaze-contingent eye-tracking.  The inverted encoding modeling (IEM) technique was applied to reconstruct location-selective responses from alpha-band activity during encoding and retention periods.  Eye movements disrupted the maintenance of stored location representations in alpha, and there was no evidence to suggest a distinct pattern of alpha activation emerged to support the location-specific representation during eye-movements.  However, location representations could be recovered when eye-gaze stabilized both between and after eye-movements.  These results suggest that eye-movements disrupt maintenance of spatial information in alpha in a manner that is consistent with other acute interruptions to continuous visual input. 

## Data:

Preprocessed EEG data files (ready for IEM) can be downloaded <b>[here](https://ucsb.box.com/s/17pfhwl4a8i1gnzyuwhz4rh2hplszhaq)</b>.  

## Scripts


### Meta Scripts

Use `WRAPPER.m` function to automatically process Trial and EEG data


### Trial Data

`Beh_Rename_Files` Remove unique date codes from trial data files (makes them easier to work with)

`Beh_Merge_Blocks` Rename and merge trial info blocks

`Beh_WM_Modelling` Run MEMFIT mixture model (w/bias) on behavioral data

`Beh_Convert_Format_MAT_to_R` Convert behavioral data from wide (MATLAB)to long (R)format


### EEG Preprocessing

`EEG_Preprocessing1` Import, filter, sync EEG with Trial Info
 
`EEG_Preprocessing2` Reject artifacts, split data into conditions

`badChannelInfo` Identify bad electrodes for each subject


### IEM Individual Conditions w/Blocked Approach

`IEM_Settings` Contains general IEM settings [BLOCKED MODEL ONLY]. Very important!

`IEM` Run all IEMs, save TFs for high temporal resolution "within" IEMs, save weights for lower temporal resolution "cross" IEMs
 
`IEM_Compute_Compile_Slopes_Within` Generate slopes for "within" IEM results

`IEM_Cross_TT` Run IEM on the training/testing weights for all conditions and timepoints

`IEM_All_Freqs` +job Run IEM across multiple frequencies (e.g. 4-30 Hz)

`IEM_Compute_Compile_Slopes_All_Freqs` Generate slopes for "within" IEM results computed across multiple freqs (e.g. 4-30 Hz)

`IEM_Cross_TT` +job Run cross-condition IEM training/testing

`IEM_Cross_TT_Calculate_slope` Calculate slopes for cross-condition IEMs

`IEM_Cross_TT_Compile_Slopes` Compile cross-condition slopes from all subjects into a single master file

`IEM_Eye` Run fixation locked IEM (doesn't really work, abandoned)

`IEM_Plot_Cross_TT_With_Threshold` Plots cross condition IEM results with threshold (Bayes or Frequentist)

`IEM_Cross_Session_Gen`+job Run cross session (i.e. days 1 and 2) IEM and get training/testing sets for each session/condition (control analysis)

`IEM_Cross_Session_TT`+job Run cross-session (i.e. days 1 and 2) training and testing (control analysis)

`IEM_Cross_Session_TT_Compile_Slopes` Get slopes from cross-session training/testing and compile into a single "master" file (control analysis)

`IEM_Plot_Cross_Session_TT.m` Plots cross session analysis (control analysis)

`SNR_Analysis_Compute`+job Calculate SNR for each condition

`SNR_Analysis_Get_Non_Alpha`+job Get non-alpha bandpass data for SNR analysis

`EEG_Compile_Bandpassed_Data` Grab bandpassed data from individual subjects and compile into a single file

`Compute_Alpha_Lateralization` Does exactly that (just for stim and ret periods)

`Compute_Alpha_Lateralization_Multiple_EMs` Does exactly that (for stim and multiple chunks of EM data)

`EEG_Extract_Preprocessed_Data_For_ERP_Analysis` +job Gets raw signal data for ERP analyses (control analysis)

`EEG_Compute_ERPs` Generate ERPs (control analysis)


### IEM Fixed Model w/ Leave One Absent Approach

`IEM_Fixed_Model_LOA.m`+job Run fixed model with leave one out cross-validation.  This doesn't seem to be workign well.

`IEM_Fixed_Model_LOA_Settings.m` IEM settings for ^

`Plot_IEM_Within_Surf_Fixed_Split_Conds.m` Create surface plots for fixed model LOA results (very rough script)


### Eye Data Analysis

`Process_Eye_Data1` Epoch eye data around stimulus trigger 102 and create matrix of all trials

`Process_Eye_Data2` Sync epoched eye data (all trials) with trial data

`Process_Eye_Data3` Remove broken trials and epoch around targets

`Eye_Movement_Analysis_Euclid` Finds euclidian distance between eye movements and fixation dot cues in "move" conditions

`Eye_Movement_Analysis_Euclid_Plot` Plots euclidian distance (error) results for ^^

`Process_Eye_Data4`Manually select stable eye time segments post EM1, EM2 and EM3 cues and save indices (this is for the fixation-locked IEM analysis that we didn't include in paper)


### Multimodal Scripts

`Correlate_Beh_With_CTFs` Look at relationship between CTF slope and behavior


### Plotting scripts

`Plot_Alpha_Topographies` Plot alpha topos (circular layout)

`Plot_Alpha_Topographies_Static_2_Periods` Plot alpha topos for stimulus and retention periods

`Plot_Alpha_Topographies_Static_4_Periods` Plot alpha topos for stimulus and EM1, EM2 and EM3 chunks of retention period

`Plot_Behavior_Modeling_Results` Plots mixture model SD and Guess rate results

`Plot_Alpha_Topographies_GIFs` Generate gifs showing alpha lateralization shifts over time on trial

`Plot_Trials_Schematics_With_Alpha_Topos` Generate schematics showing samples of alpha lateralization and corresponding target position on screen

`Plot_IEM_Cross_TT` Plot results of IEM cross training/testing analyses

`Plot_IEM_Within_Surf` Generate surface plots for CTFs

`Plot_SNR_Results` Compare SNR across conditions (control)


### Miscellaneous

`badBlock_beh` Fix messed up behavioral data

`badCB_beh` Fix messed up counterbalancing data for certain subjects

`EEG_Fix_Bad_Config_Files` Fix EEG files that were recorded with incorrect sample rate and channel lables in ActiView

`EEG_Merge_Broken_Files` Merge EEG files where continuous recording (in ActiView) was interrupted mid experiment e.g. subject bathroom break

`Eye_Diagnosis` Troubleshoot eye/trial data sync issues

`Eye_Location_Plot` Quickly plot a single subject's eye-movements

`Plot_EYE_Data` Plot eye data for single subs (old, remove)

`Plot_EYE_Data_Color` Kamryn? )

`Plot_EYE_Nwe` Kamryn?


### Dependencies (third party, should all really be placed in separate folder)

EEGLAB toolbox (v.2019_1)

`gif` Creates gifs

`hline` Draws horizontal lines in plots

`vline` Draws vertical lines in plots

`shadedErrorBar` Generates shaded error bars

`visAngleCalculatePix` Convert pixels to degrees visual angle

`wshift` shifts a vector or matrix


## Notes for GitHub/Coding Style Workshop (remove)

Give an overview of how to set up a github folder, how to add/commit/push etc. 

Show how I deal with having stuff on cluster and local machines (directory naming v.important here, needs to be unified)

How to deal with local scripts that are not run on cluster? e.g. R scripts?

Integrating plots with adobe illustrator (show auto updating)

Talk about MATLAB GitHub interaction e.g. useful for renaming files, deleting etc. But can also mysteriously not work for some stuff.

Script titles (add a "project" line, also "date created, date last updated")

Function names are capitalized, underscored