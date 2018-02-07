function asa = getAdditionalSubmitArguments(props)

% https://www.tuegrid-doc.uni-tuebingen.de/dokuwiki/lib/exe/fetch.php?id=hpc-bw%3Ahpc-bw&cache=cache&media=hpc-bw:batch_doc.pdf

% Copyright 2010-2015 The MathWorks, Inc.

asa = '';
currFilename = mfilename;

[~,flag] = strtok(props.MatlabArguments,'p');
if strcmp(flag,'parallel')
    w = props.NumberOfTasks;
else
    w = 1;
end

% Default PPN
ppn = 8;


%% REQUIRED


%% OPTIONAL

% Walltime
wt = ClusterInfo.getWallTime();
if isempty(wt)==false
    asa = [asa ' -l walltime=' wt];
end

% Physical Memory used by the entire job
mu = ClusterInfo.getMemUsage();
if isempty(mu)==false
    asa = [asa ' -l mem=' mu];
end

qn = ClusterInfo.getQueueName();
if isempty(qn)==false
    asa = [asa ' -q ' qn];
    switch qn
        % Add geotest, teu-extralong, tue-long, tue-test, tue-short?
        case 'cfc'
            ppn = 20;
        case 'chem01'
            ppn = 64;
        case 'cpt'
            ppn = 24;
        case 'hpc-uni'
            % MW: This can either be 4 or 8
            ppn = 4;
        case 'linguist'
            ppn = 8;
        case 'unicore'
            % MW: Told 16, should this be 8?
            ppn = 8;
        case 'user'
            % Routing Q?
            ppn = 8;
        case 'workshop'
            ppn = 8;
    end
end

% Number of nodes & procs per node
u_ppn = ClusterInfo.getProcsPerNode();
if isempty(u_ppn)==false
    % Override ppn
    ppn = u_ppn;
end
ppn = min(w,ppn);

% Nodes
numberOfNodes = ceil(w/ppn);

asa = sprintf('%s -l nodes=%d:ppn=%d', asa, numberOfNodes, ppn);
dctSchedulerMessage(4, '%s: Requesting %d nodes with %d processors per node', currFilename, ...
    numberOfNodes, ppn);

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
%      my @FLEXlmCmds = ("/opt/bwhpc/common/math/matlab/R2015a/etc/glnxa64/lmutil lmstat -a -c /opt/bwhpc/common/math/matlab/R2014b/licenses/network.lic");
%
%  3. On u-003-scfe03, modify $MOABHOMEDIR/.../moab.cfg
%
%      RMCFG[FLEXlm1]     TYPE=NATIVE RESOURCETYPE=LICENSE
%      RMCFG[FLEXlm1]     CLUSTERQUERYURL=exec://$TOOLSDIR/flexlm/license.mon.matlab.flexLM.pl
%
asa = sprintf('%s -l software=MATLAB_Distrib_Comp_Engine:%d',asa,w);

% Catch-all
udo = ClusterInfo.getUserDefinedOptions();
if isempty(udo)==false
    asa = [asa ' ' udo];
end

asa = strtrim(asa);
