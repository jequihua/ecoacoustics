function jobID = extractJobId(cmdOut)
% Extracts the job ID from the qsub command output for PBS

% Copyright 2010-2015 The MathWorks, Inc.

% jobId should be in the following format:
% 123.server-name

% The second piece of the regexp must pick out valid (fully-qualified) server names
% % jobID = regexp(cmdOut, '[0-9\[\]]+\.[a-zA-Z0-9-\.]*', 'match', 'once');

% As of Nov 2014, the output is a simple number job ID.
% % jobID = strsplit(cmdOut,'.');
% % jobID = jobID{2};

jobID = strtrim(cmdOut);
dctSchedulerMessage(0, '%s: Job ID %s was extracted from qstat output %s.', mfilename, jobID, cmdOut);
