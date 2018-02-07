function configCluster
% Configure MATLAB for users running remotely to the cluster.
%   1. Imports partial profile
%   2. Determines and creates (if needed) local job storage location
%   3. Determines remote job storage location
%   4. Updates profile and sets it as the default

% Copyright 2013-2015 The MathWorks, Inc.

% The location of MATLAB as compared to the cluster
type = 'remote';

% The version of MATLAB being supported
release = ['R' version('-release')];

% Baden-Wurttemberg
site = 'bw';

% The location of the cluster setting file(s)
pfile_loc = fileparts(mfilename('fullpath'));

% Listing of setting file(s).  Derive the specific one to use.
pfiles = fullfile(pfile_loc,'*.settings');
pfile_listing = dir(pfiles);
len = length(pfile_listing);
if len==0
    error('Failed to find profiles.  Contact your System Administrator.')
elseif len==1
    pfile = pfile_listing.name;
else
    pfile = lExtractPfile(pfile_listing);
end

% From the settings file, derive the name of the cluster
cluster = strtok(pfile,'_');

% We want to ensure that the user is running the correct version of MATLAB.
% We'll use the release string on the setting files (which matches the
% version of MDCS on the cluster) with the release of MATLAB.  If they're
% not the same, we'll throw a warning to the user.
expression = 'r20[0-9]+[ab]';
mdcs_release = regexp(pfile,expression,'match');
if strcmpi(mdcs_release,release)==false
    warning('This version of MATLAB does not match the supported MDCS release on %s.',upper(cluster))
end

% Create the user's local Job Storage Location folder
rootd = lGetLocalRoot();
jfolder = fullfile(rootd,'MdcsDataLocation',site,cluster,release,type);
if exist(jfolder,'dir')==false
    [status,err,eid] = mkdir(jfolder);
    if status==false
        error(eid,err)
    end
end

% From the settings file, derive the name of the Profile
[~,profile] = fileparts(pfile);

% Delete the old profile (if it exists)
profiles = parallel.clusterProfiles();
idx = strcmp(profiles,profile);
ps = parallel.Settings;
ws = warning;
warning off %#ok<WNOFF>
ps.Profiles(idx).delete
warning(ws)

% Import the new profile
p = parallel.importProfile(fullfile(pfile_loc,pfile));

% User's remote Job Storage Location folder
rootd = lGetRemoteRoot(cluster);
rjsl = [rootd '/MdcsDataLocation/' site '/' cluster '/' release '/' type];

% Get a handle to the profile
c = parcluster(p);
c.JobStorageLocation = jfolder;
c.IndependentSubmitFcn{3} = rjsl;
c.CommunicatingSubmitFcn{3} = rjsl;
c.saveProfile

% Set as default profile
parallel.defaultClusterProfile(p);

% Since we've imported a new profile, in all likelihood, any of the current
% ClusterInfo settings are void
disp('Clearing all ClusterInfo settings.')
ClusterInfo.clear

ns = ['profiles.' site '.' cluster];
ClusterInfo.setNameSpace(ns)

lNotifyUserOfCluster(upper(cluster))

% % Validate if you want to
% profiles = parallel.clusterProfiles();
% idx = strcmp(profiles,profile);
% ps = parallel.Settings;
% ps.Profiles(idx).validate

end


function pfile = lExtractPfile(pl)
% Display profile listing to user to select from
len = length(pl);
for pidx = 1:len
    name = pl(pidx).name;
    names{pidx,1} = strtok(name,'_'); %#ok<AGROW>
end

selected = false;
while selected==false
    for pidx = 1:len
        fprintf('\t[%d] %s\n',pidx,names{pidx});
    end
    idx = input(sprintf('Select a cluster [1-%d]: ',len));
    selected = idx>=1 && idx<=len;
end
pfile = pl(idx).name;

end


function r = lGetLocalRoot()

if isunix
    % Some Mac user's have noticed that the [/private]/tmp
    % directory gets cleared when the system is reboot, so for UNIX
    % in general, let's just use the user's local home directory.
    uh = java.lang.System.getProperty('user.home');
    uh = uh.toLowerCase;
    
    r = char(uh);
else
    % If this returns an empty string (some how user name is not defined),
    % it's a no-op for FULLFILE, so there's no strong need to error out.
    un = java.lang.System.getProperty('user.name');
    un = un.toLowerCase;
    un = char(un);
    
    r = fullfile(tempdir,un);
end

end


function r = lGetRemoteRoot(cluster)

r = input(['Home directory on ' upper(cluster) ' (e.g. /home/joe): '],'s');
if isempty(r)
    error(['Failed to configure cluster: ' cluster])
end

end


function lNotifyUserOfCluster(cluster)

switch lower(cluster)
    case 'freiburg'
        fprintf('\nBefore submitting a job to %s, you must specify the wall time.\n',cluster);
        fprintf('\n\t\t>> %% E.g. set wall time to 1 hour\n\t\t>> ClusterInfo.setWallTime(''01:00:00'')\n\n');
    case {'kit','tuebingen'}
        fprintf('\nBefore submitting a job to %s, you can specify the wall time.\n',cluster);
        fprintf('\n\t\t>> %% E.g. set wall time to 1 hour\n\t\t>> ClusterInfo.setWallTime(''01:00:00'')\n\n');
    otherwise
        error('Unsupported cluster %s',cluster)
end

end
