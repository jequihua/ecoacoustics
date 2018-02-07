function jobID = extractJobId(cmdOut)
% Extracts the job ID from the qsub command output for PBS

% Copyright 2010-2015 The MathWorks, Inc.

% jobId should be in the following format:
% [Moab.]1234

[notneeded,jobID] = strtok(cmdOut,'.');
if isempty(jobID)==true
    % Didn't find a '.', so it wasn't of the
    % format Moab.####, just ####
    jobID = strtrim(notneeded);
else
    % Found a '.', so remove the leading '.'
    jobID = jobID(2:end);
end
dctSchedulerMessage(0, '%s: Job ID %s was extracted from qstat output %s.', mfilename, jobID, cmdOut);
