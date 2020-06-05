%{
Purpose: run IEM script in parallel
Author: Tom Bullock
Date: 05.08.20

%}

% select subjects
subjects = [1:7,9:14,16:20,22:27,31];


% run in serial for debugging (0) or parallel (1)
runInParallel=1;

if runInParallel
    s = parcluster();
    job = createJob(s);
end

for iSub = 1:length(subjects)
    sjNum = subjects(iSub);
    if runInParallel
        createTask(job,@EEG_Extract_Preprocessed_Data_For_ERP_Analysis,0,{sjNum});
    else
        EEG_Extract_Preprocessed_Data_For_ERP_Analysis(sjNum) 
    end
end

if runInParallel
    submit(job)
end