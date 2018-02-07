function state = getJobStateFcn(cluster, job, state)
%GETJOBSTATEFCN Gets the state of a job from PBS
%
% Set your cluster's GetJobStateFcn to this function using the following
% command:
%     set(cluster, 'GetJobStateFcn', @getJobStateFcn);

% Copyright 2010-2015 The MathWorks, Inc.

% Store the current filename for the errors, warnings and dctSchedulerMessages
currFilename = mfilename;
if ~isa(cluster, 'parallel.Cluster')
    error('parallelexamples:GenericPBS:SubmitFcnError', ...
        'The function %s is for use with clusters created using the parcluster command.', currFilename)
end
if cluster.HasSharedFilesystem
    error('parallelexamples:GenericPBS:SubmitFcnError', ...
        'The submit function %s is for use with nonshared filesystems.', currFilename)
end

ns = ClusterInfo.getNameSpace();
getRemoteConnection = str2func([ns '.getRemoteConnection']);

% Get the information about the actual cluster used
data = cluster.getJobClusterData(job);
if isempty(data)
    % This indicates that the job has not been submitted, so just return
    dctSchedulerMessage(1, '%s: Job cluster data was empty for job with ID %d.', currFilename, job.ID);
    return
end

try
    hasDoneLastMirror = data.HasDoneLastMirror;
catch err
    ex = MException('parallelexamples:GenericPBS:FailedToRetrieveRemoteParameters', ...
        'Failed to retrieve remote parameters from the job cluster data.');
    ex = ex.addCause(err);
    throw(ex);
end
% Shortcut if the job state is already finished or failed
jobInTerminalState = strcmp(state, 'finished') || strcmp(state, 'failed');
% and we have already done the last mirror
if jobInTerminalState && hasDoneLastMirror
    return;
end
try
    clusterHost = data.RemoteHost;
    remoteJobStorageLocation = data.RemoteJobStorageLocation;
catch err
    ex = MException('parallelexamples:GenericPBS:FailedToRetrieveRemoteParameters', ...
        'Failed to retrieve remote parameters from the job cluster data.');
    ex = ex.addCause(err);
    throw(ex);
end
remoteConnection = getRemoteConnection(cluster, clusterHost, remoteJobStorageLocation);
try
    jobIDs = data.ClusterJobIDs;
catch err
    ex = MException('parallelexamples:GenericPBS:FailedToRetrieveJobID', ...
        'Failed to retrieve clusters''s job IDs from the job cluster data.');
    ex = ex.addCause(err);
    throw(ex);
end

% Get the full display from checkjob
commandToRun = sprintf('checkjob --xml %s', sprintf('%s ', jobIDs{:}));
dctSchedulerMessage(4, '%s: Querying cluster for job state using command:\n\t%s', currFilename, commandToRun);

try
    % We will ignore the status returned from the state command because
    % a non-zero status is returned if the job no longer exists
    % Execute the command on the remote host.
    [~, cmdOut] = remoteConnection.runCommand(commandToRun);
catch err
    ex = MException('parallelexamples:GenericPBS:FailedToGetJobState', ...
        'Failed to get job state from cluster.');
    ex = ex.addCause(err);
    throw(ex);
end

[clusterState, br] = iExtractJobState(cmdOut, numel(jobIDs));
dctSchedulerMessage(6, '%s: State %s was extracted from cluster output.\n', currFilename, clusterState);
if isempty(br)==false
    for bidx = 1:length(br)
        dctSchedulerMessage(6, '%s: BlockReason is:\n', currFilename, br{bidx});
    end
end

% If we could determine the cluster's state, we'll use that, otherwise
% stick with MATLAB's job state.
if ~strcmp(clusterState, 'unknown')
    state = clusterState;
end
% Decide what to do with mirroring based on the cluster's version of job state and whether or not
% the job is currently being mirrored:
% If job is not being mirrored, and job is not finished, resume the mirror
% If job is not being mirrored, and job is finished, do the last mirror
% If the job is being mirrored, and job is finished, do the last mirror.
% Otherwise (if job is not finished, and we are mirroring), do nothing
isBeingMirrored = remoteConnection.isJobUsingConnection(job.ID);
isJobFinished = strcmp(state, 'finished') || strcmp(state, 'failed');
if ~isBeingMirrored && ~isJobFinished
    % resume the mirror
    dctSchedulerMessage(4, '%s: Resuming mirror for job %d.', currFilename, job.ID);
    try
        remoteConnection.resumeMirrorForJob(job);
    catch err
        warning('parallelexamples:GenericPBS:FailedToResumeMirrorForJob', ...
            'Failed to resume mirror for job %d.  Your local job files may not be up-to-date.\nReason: %s', ...
            err.getReport);
    end
elseif isJobFinished
    dctSchedulerMessage(4, '%s: Doing last mirror for job %d.', currFilename, job.ID);
    try
        remoteConnection.doLastMirrorForJob(job);
        % Store the fact that we have done the last mirror so we can shortcut in the future
        data.HasDoneLastMirror = true;
        cluster.setJobClusterData(job, data);
    catch err
        warning('parallelexamples:GenericPBS:FailedToDoFinalMirrorForJob', ...
            'Failed to do last mirror for job %d.  Your local job files may not be up-to-date.\nReason: %s', ...
            err.getReport);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state, br] = iExtractJobState(checkjobOut, numJobs)
% Function to extract the job state from the output of checkjob

[lstate, br] = iGenerateXmlTree(checkjobOut,numJobs);

% For PBSPro, the pending states are HQSTUW, the running states are BRE
% For Torque, the pending states are HQW, the running states are RE
% Based on https://computing.llnl.gov/tutorials/moab/#JobStates
numPending = numel(regexp(lstate, 'BatchHold|SystemHold|UserHold|Deferred|Idle|Staging|Suspended|Migrated'));
numRunning = numel(regexp(lstate, 'Running|Starting'));
% Canceling from what?  Could either be from Idle or Running.
numFinished = numel(regexp(lstate, 'Completed|Canceling|Removed|Vacated|NotQueued'));

% If all of the jobs that we asked about have finished, then we know the job has finished.
if numFinished == numJobs
    state = 'finished';
    return;
end

% Any running indicates that the job is running
if numRunning > 0
    state = 'running';
    return;
end

% We know numRunning == 0 so if there are some still pending then the
% job must be queued again, even if there are some finished
if numPending > 0
    state = 'queued';
    return
end

state = 'unknown';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state, br] = iGenerateXmlTree(q,numJobs)

[~,tn] = fileparts(tempname);
file = [tn '_job.xml'];
fid = fopen(file,'wt');
if fid<0
    emsg = ['Failed to open file ''' file ''' for writing to.'];
    dctSchedulerMessage(4, '%s: ', emsg);
    error(emsg)
end
cu = onCleanup(@()iCleanup(fid));

% Write out XML file
q = strtrim(q);
fprintf(fid,'%s',q);

% Read XML file
xDoc = xmlread(file);

% Search for State and BlockReason
job = xDoc.getElementsByTagName('job');
state = [];
br = {};
for jidx = 0:numJobs-1
    item = job.item(jidx);
    try
        the_state = char(item.getAttribute('State'));
    catch
        % Let's check if Moab threw a (specific) error code.  I can't find
        % a list of potential codes, but one to look for is
        %
        %   <Error code="700">server rejected request - invalid job specified: ######</Error>
        %
        % Since we've asked for a specific job ID, it must have existed at
        % one point.  We can get this error code if we wait too long to
        % query for the state, after it's finished.
        %
        % We can't assume all error codes are completed, but in this case,
        % we will.  We're also assuming that there's only one job.  If
        % there's more than one, we wouldn't know which job was invalid.
        % The problem becomes that we may never get out of this "queued"
        % state.  In all likelihood, if we've specified an invalid job ID,
        % there's only one job.
        if numJobs==1 && isempty(strfind(q,'invalid job specified'))==false
            the_state = 'Completed';
        else
            the_state = 'unknown';
        end
    end
    state = [state the_state]; %#ok<AGROW>
    try
        the_block_reason = char(item.getAttribute('BlockReason'));
    catch
        the_block_reason = '';
    end
    br{jidx+1,1} = the_block_reason; %#ok<AGROW>
end

function iCleanup(fid)
file = fopen(fid);
if isempty(file)==false
    fclose(fid);
    delete(file)
end
