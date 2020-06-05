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
    s.ResourceTemplate = '--ntasks-per-node=6 --mem=65536';
    job = createJob(s);
end


if runInParallel
    createTask(job,@SNR_Analysis_Compute,0,{subjects});
else
    SNR_Analysis_Compute(subjects)
end


if runInParallel
    submit(job)
end