function showStart(job)
%
% If the following is displayed
%
%   INFO:  cannot determine start time for job ########
%
% then there's a good chance the job has already run.

% Copyright 2015 The MathWorks, Inc.
% Raymond S. Norris (raymond.norris@mathworks.com)

narginchk(1,1)
if numel(job)>1
    error('Must only supply one job.')
end

if ~isa(job,'parallel.job.CJSIndependentJob') ...
        && ~isa(job,'parallel.job.CJSCommunicatingJob')
    error('Must provide Independent or Communicating Job')
end

ns = ClusterInfo.getNameSpace();
getRemoteConnection = str2func([ns '.getRemoteConnection']);

% Get remote connection
%
%  An alternative to this could just be:
%
%   remoteConnection = job.Parent.UserData.RemoteConnection;
%
cluster = job.Parent;
data = cluster.getJobClusterData(job);

try
    clusterHost = data.RemoteHost;
    remoteJobStorageLocation = data.RemoteJobStorageLocation;
catch err
    ex = MException('FailedToRetrieveRemoteParameters', ...
        'Failed to retrieve remote parameters from the job cluster data.');
    ex = ex.addCause(err);
    throw(ex);
end

remoteConnection = getRemoteConnection(cluster, clusterHost, remoteJobStorageLocation);

jid = schedID(job);
commandToRun = sprintf('showstart -e all %s',jid);

try
    [~, cmdOut] = remoteConnection.runCommand(commandToRun);
catch err
    ex = MException('FailedToGetJobStart', ...
        'Failed to get job start from cluster.');
    ex = ex.addCause(err);
    throw(ex);
end

disp(strtrim(cmdOut))

end
