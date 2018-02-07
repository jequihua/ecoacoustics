function cleanJobStorageLocation(c)

% Copyright 2014-2015 The MathWorks, Inc.
% Raymond S. Norris (raymond.norris@mathworks.com)

narginchk(1,1)

if isa(c,'parallel.cluster.Generic')==false
    error('Must supply cluster object.')
end

% Store the current filename for the errors, warnings and dctSchedulerMessages
currFilename = mfilename;

% Delete local job storage location
jsl = c.JobStorageLocation;
if exist(jsl,'dir')==7
    [success,emsg,eid] = rmdir(jsl,'s');
    if success==false
        error(eid,emsg)
    end
end

% Create local job storage location
[success,emsg,eid] = mkdir(jsl);
if success==false
    error(eid,emsg)
end

if numel(c.IndependentSubmitFcn)<3
    % Running on the cluster, so we don't need to worry about remote.
    % Return early.
    return
end

% Get name of cluster and remote job storage location to get handle to
% remote connection
clusterHost = c.IndependentSubmitFcn{2};
remoteJobStorageLocation = c.IndependentSubmitFcn{3};

ns = ClusterInfo.getNameSpace();
getRemoteConnection = str2func([ns '.getRemoteConnection']);

% Get remote connection
remoteConnection = getRemoteConnection(c, clusterHost, remoteJobStorageLocation);

% Delete old remote job storage location
try
    commandToRun = sprintf('rm -rf %s',remoteJobStorageLocation);
    % Execute the command on the remote host.
    [cmdFailed, cmdOut] = remoteConnection.runCommand(commandToRun);
catch err
    cmdFailed = true;
    cmdOut = err.message;
end
if cmdFailed
    dctSchedulerMessage(1, '%s: Failed to delete remote job storage location on cluster.  Reason:\n\t%s', currFilename, cmdOut);
end

% We could recreate the remote job storage location, but the next time we
% submit a job, it'll automatically get created for us.  Additionally, we
% need to have the storage metadata file recreated anyway.  The next job we
% create will generate this warning:

% Warning: The storage metadata file did not exist. Recreating it.

% % try
% %     commandToRun = sprintf('mkdir -p %s',remoteJobStorageLocation);
% %     % Execute the command on the remote host.
% %     [cmdFailed, cmdOut] = remoteConnection.runCommand(commandToRun);
% % catch err
% %     cmdFailed = true;
% %     cmdOut = err.message;
% % end
% % if cmdFailed
% %     dctSchedulerMessage(1, '%s: Failed to create remote job storage location on cluster.  Reason:\n\t%s', currFilename, cmdOut);
% % end
