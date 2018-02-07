function submitString = getSubmitString(jobName, quotedLogFile, quotedCommand, ...
    varsToForward, props)
%GETSUBMITSTRING Gets the correct msub command for a PBS cluster

% Copyright 2010-2015 The MathWorks, Inc.

envString = strjoin(varsToForward', ',');

ns = ClusterInfo.getNameSpace();
getAdditionalSubmitArguments = str2func([ns '.getAdditionalSubmitArguments']);
additionalSubmitArgs = getAdditionalSubmitArguments(props) %#ok<NOPRT>

% Submit to PBS using msub. Note the following:
% "-N Job#" - specifies the job name
% "-j oe" joins together output and error streams
% "-o ..." specifies where standard output goes to
% envString has the "-v 'NAME,NAME2'" piece.
submitString = sprintf('msub -N %s -j oe -o %s -v %s %s %s', ...
    jobName, quotedLogFile, envString, additionalSubmitArgs, quotedCommand);
