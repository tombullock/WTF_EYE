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
    disp('PARALLEL PROCESSING ACTIVATED!'); pause(3);
    s = parcluster();
    job = createJob(s);
else
    disp('SERIAL PROCESSING...'); pause(3);
end

for iSub = 1:length(subjects)
    sjNum = subjects(iSub);
    if runInParallel
        createTask(job,@IEM_Cross_Session_TT,0,{sjNum});
    else
        IEM_Cross_Session_TT(sjNum) 
    end
end

if runInParallel
    submit(job)
end