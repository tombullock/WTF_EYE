# SCAM (Space, Color, Attention and Memory) 

## Project Goals: 

In our previous "WTF" study we demonstrated that Process 2 emerges when participants are required to remember a spatial location in the absence of continuous sensory input (closed eyes).
In the current study we test whether Process 2 emerges for remembered spatial locations when sensory input is disrupted (by eye movements in the retention period) rather than cut off.
If true, then this would suggest that Process 2 is not just specific to eye closure.

In additon to eye movement, we also manipulate whether the task is to remember location (spatial) or color (non-spatial).


## Scripts

Use `WRAPPER.m` function to automatically process Trial and EEG data

### Process Trial Data (includes behavioral data) 

`Beh_Merge_Blocks` Rename and merge trial info blocks


### Process EEG Data

`EEG_Preprocessing1` Import, filter, sync EEG with Trial Info
 
`EEG_Preprocessing2` Reject artifacts, split data into conditions

`IEM` Run all IEMs, save TFs for high temporal resolution "within" IEMs, save weights for lower temporal resolution "cross" IEMs
 
`IEM_Compute_Compile_Slopes_Within` Generate slopes for "within" IEM results

`IEM_Cross_TT` Run IEM on the training/testing weights for all conditions and timepoints


### Plot/Analyze Data

`Beh_WM_Modelling` Run mixture model (w/bias) on data

`Plot_IEM_Within_Surf` Generate IEM 3D surface plots (training/testing on the same timepoint within participants)

`Plot_IEM_Cross_TT` Generate IEM heatmaps (training/testing across all timespoints and across all conditions)

