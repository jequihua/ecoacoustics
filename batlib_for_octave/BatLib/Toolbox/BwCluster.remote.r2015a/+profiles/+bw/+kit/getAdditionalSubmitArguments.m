function asa = getAdditionalSubmitArguments(props)

% http://www.bwhpc-c5.de/wiki/index.php/BwUniCluster_Batch_Jobs

% Copyright 2010-2015 The MathWorks, Inc.

asa = '';

[~,flag] = strtok(props.MatlabArguments,'p');
if strcmp(flag,'parallel')
    w = props.NumberOfTasks;
else
    w = 1;
end


%% REQUIRED


%% OPTIONAL

% Walltime
wt = ClusterInfo.getWallTime();
if isempty(wt)==false
    asa = [asa ' -l walltime=' wt];
end

% Physical Memory used by a single process of the job
mu = ClusterInfo.getMemUsage();
if isempty(mu)==false
    asa = [asa ' -l pmem=' mu];
end

% Number of procs per node
ppn = ClusterInfo.getProcsPerNode();
if isempty(ppn)==false
    asa = [asa ' -l tpn=' num2str(ppn)];
    requireMultiNodeFlag = true;
else
    requireMultiNodeFlag = false;
end

asa = sprintf('%s -l procs=%d', asa, w);

% Queue name
if requireMultiNodeFlag==true
    % If the user specifies the procs per node, then it's safe to assume
    % they want more than one node.  Therefore, automatically assign
    % queue name to 'multinode'
    qn = 'multinode';
else
    % http://www.bwhpc-c5.de/wiki/index.php/Batch_Jobs_-_bwUniCluster_Features#msub_-q_queues
    qn = ClusterInfo.getQueueName();
    % %     switch qn
    % %         case 'develop'
    % %             % max: nodes=1:ppn=16, 30 minutes
    % %             % default: nodes=1:ppn=1, 30+ min, 4 GB
    % %             ppn = 16;
    % %         case 'singlenode'
    % %             % max: nodes=1:ppn=16, 3 days
    % %             % default: nodes=1:ppn=1, 30+ min, 4 GB
    % %             ppn = 16;
    % %         case 'multinode'
    % %             % max: nodes=16:ppn=16, 2 days
    % %             % default: nodes=2:ppn=1, 10 min, 4 GB
    % %             singlenode = false;
    % %             ppn = 16;
    % %         case 'verylong'
    % %             % max: nodes=1:ppn=16, 6 days
    % %             % default: nodes=1:ppn=1, 3+ days, 4 GB
    % %             ppn = 16;
    % %         case 'fat'
    % %             % max: nodes=1:ppn=32, 1 day
    % %             % default: nodes=1:ppn=1, 10 min, 32 GB
    % %             ppn = 32;
    % %     end
end
if isempty(qn)==false
    asa = [asa ' -q ' qn];
end

% Reservation
res = ClusterInfo.getReservation();
if isempty(res)==false
    asa = [asa ' -l advres=' res];
end

% Email notification
ea = ClusterInfo.getEmailAddress();
if isempty(ea)==false
    % User wants to be emailed the job status
    asa = [asa ' -M ' ea ' -m abe'];
end

% Every job is going to require a certain number of MDCS licenses.
% The Native RM job submission language provides a direct method of
% license specification.
%
% Moab:
%  1. Copy
%          $MOABHOMEDIR/tools/flexlm/license.mon.flexLM.pl
%     to
%          $MOABHOMEDIR/tools/flexlm/license.mon.matlab.flexLM.pl
%
%  2. Modify license.mon.matlab.flexLM.pl
%
%      my @FLEXlmCmds = ("/opt/bwhpc/common/math/matlab/R2014b/etc/glnxa64/lmutil lmstat -a -c /opt/bwhpc/common/math/matlab/R2014b/licenses/network.lic");
%
%  3. On uc1n996, modify $MOABHOMEDIR/etc/moab.cfg
%
%      RMCFG[FLEXlm1]     TYPE=NATIVE RESOURCETYPE=LICENSE
%      RMCFG[FLEXlm1]     CLUSTERQUERYURL=exec://$TOOLSDIR/flexlm/license.mon.matlab.flexLM.pl
%
%{
% TODO
asa = sprintf('%s -l software=MATLAB_Distrib_Comp_Engine+%d',asa,w);
%}

% Catch-all
udo = ClusterInfo.getUserDefinedOptions();
if isempty(udo)==false
    asa = [asa ' ' udo];
end

asa = strtrim(asa);
