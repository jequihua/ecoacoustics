classdef (Sealed) ClusterInfo < handle
    %CLUSTERINFO ClusterInfo class.
    %
    %   Provides cluster information about:
    %      Arch
    %      ClusterHost
    %      DataParallelism
    %      DebugMessagesTurnedOn
    %      DiskSpace
    %      EmailAddress
    %      GpusPerNode
    %      MemUsage
    %      PrivateKeyFile
    %      PrivateKeyFileHasPassPhrase
    %      ProcsPerNode
    %      ProjectName
    %      QueueName
    %      RequireExclusiveNode
    %      Reservation
    %      SshPort
    %      UseGpu
    %      UserDefinedOptions
    %      UserNameOnCluster
    %      WallTime
    %
    %   To add properties, both a get and set function need to be coded.
    %   For example (error handling has been removed for simplicity):
    %
    %      function setPropertyName(val)
    %      setpref(ClusterInfo.Group,PropertyName,val)
    %
    %      function val = getPropertyName()
    %      val = getpref(ClusterInfo.Group,PropertyName);
    %
    %   The reference to PropertyName in the function name would be
    %   replaced with the actual name of the property.  The reference to
    %   PropertyName in the calls to setpref and getpref would be replaced
    %   with a quoted reference of the actual property name.
    %
    %   Example:
    %     % Prior to calling a submit script, set a property, as such:
    %     ClusterInfo.setWallTime('02:30:00');
    %
    %     % Within the submit script, access the property, as such:
    %     wt = ClusterInfo.getWallTime();

    %   Copyright 2009-2015 The MathWorks, Inc.
    %   Raymond S. Norris (raymond.norris@mathworks.com)

    methods (Static,Access='private')
        function gp = Group()
            gp = 'ClusterInfo';
        end
        function va = ValidArch()
            va = {'32-bit'; '64-bit'};
        end
        function vdp = ValidDataParallelism()
            % Needs to be verified
            vdp = {'gigE'; 'ib'; 'myri'};
        end
        function ac = ArchComplex()
            ac = {'lx24-x86'; 'lx24-amd64'};
        end
    end

    methods (Access='private')
        function obj = ClusterInfo()
        end
    end

    methods (Static)

        function state()
            % Display the values of each of the fields.  The values could
            % be empty, which means the user has not set them.  We don't
            % assume a default value.

            % Pick the longest field
            len = length('PrivateKeyFileHasPassPhrase');
            spaces = 3;
            sHeaderFormat = ['%' num2str(len+spaces) 's : %s\n'];
            dHeaderFormat = ['%' num2str(len+spaces) 's : %d\n'];

            fprintf('\n');
            fprintf(sHeaderFormat, 'Arch', ClusterInfo.getArch());
            fprintf(sHeaderFormat, 'ClusterHost',  ...
                ClusterInfo.getClusterHost());
            fprintf(sHeaderFormat, 'DataParallelism',  ...
                ClusterInfo.getDataParallelism());
            fprintf(sHeaderFormat, 'DiskSpace', ...
                ClusterInfo.getDiskSpace());
            fprintf(sHeaderFormat, 'EmailAddress', ...
                ClusterInfo.getEmailAddress());
            fprintf(dHeaderFormat, 'GpusPerNode', ...
                ClusterInfo.getGpusPerNode());
            fprintf(sHeaderFormat, 'MemUsage', ...
                ClusterInfo.getMemUsage());
            fprintf(sHeaderFormat, 'PrivateKeyFile', ...
                ClusterInfo.getPrivateKeyFile());
            fprintf(sHeaderFormat, 'PrivateKeyFileHasPassPhrase', ...
                num2str(ClusterInfo.getPrivateKeyFileHasPassPhrase()));
            fprintf(dHeaderFormat, 'ProcsPerNode', ...
                ClusterInfo.getProcsPerNode());
            fprintf(sHeaderFormat, 'ProjectName', ...
                ClusterInfo.getProjectName());
            fprintf(sHeaderFormat, 'QueueName', ClusterInfo.getQueueName());
            fprintf(dHeaderFormat, 'RequireExclusiveNode', ...
                ClusterInfo.getRequireExclusiveNode());
            fprintf(sHeaderFormat, 'Reservation', ...
                ClusterInfo.getReservation());
            fprintf(dHeaderFormat, 'SshPort', ClusterInfo.getSshPort());
            fprintf(dHeaderFormat, 'UseGpu', ClusterInfo.getUseGpu());
            fprintf(sHeaderFormat, 'UserDefinedOptions', ...
                ClusterInfo.getUserDefinedOptions());
            fprintf(sHeaderFormat, 'UserNameOnCluster', ...
                ClusterInfo.getUserNameOnCluster());
            fprintf(sHeaderFormat, 'WallTime', ClusterInfo.getWallTime());
        end

        function clear()
            ns = ClusterInfo.getNameSpace();
            ClusterInfo.delete()
            ClusterInfo.setNameSpace(ns);
        end
        
        function delete()
            try
                rmpref(ClusterInfo.Group)
            catch E %#ok<NASGU>
            end
        end

        function setArch(a)
            if nargin==0 || ischar(a)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Arch must be a character string.')
            end
            va = ClusterInfo.ValidArch();
            idx = strcmp(a,va);
            if ~any(idx)
                lst = '';
                for lidx = 1:length(va)-1
                    lst = [lst va{lidx} ' ']; %#ok<AGROW>
                end
                lst = [lst va{end}];
                error('Valid arches: %s', lst)
            end
            ac = ClusterInfo.ArchComplex;
            setpref(ClusterInfo.Group,'Arch',ac{idx})
        end
        function a = getArch()
            try
                val = getpref(ClusterInfo.Group,'Arch');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            a = val;
        end

        function setClusterHost(ch)
            if nargin==0 || ischar(ch)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Cluster host must be a character string.')
            end
            setpref(ClusterInfo.Group,'ClusterHost',ch)
        end
        function ch = getClusterHost()
            try
                val = getpref(ClusterInfo.Group,'ClusterHost');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            ch = val;
        end

        function setDataParallelism(dp)
            if nargin==0 || ischar(dp)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Data parallelism must be a character string.')
            end
            vdp = ClusterInfo.ValidDataParallelism();
            idx = strcmp(dp,vdp);
            if ~any(idx)
                lst = '';
                for lidx = 1:length(vdp)-1
                    lst = [lst vdp{lidx} ' ']; %#ok<AGROW>
                end
                lst = [lst vdp{end}];
                error('Valid data parallelism: %s', lst)
            end
            setpref(ClusterInfo.Group,'DataParallelism',dp)
        end
        function dp = getDataParallelism()
            try
                val = getpref(ClusterInfo.Group,'DataParallelism');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            dp = val;
        end
        function ddp = getDefaultDataParallelism()
            vdp = ClusterInfo.ValidDataParallelism();
            ddp = vdp{1};
        end

        function setDebugMessagesTurnedOn(dmto)
            if nargin==0 || islogical(dmto)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Flag must be boolean.')
            end
            setpref(ClusterInfo.Group,'DebugMessagesTurnedOn',dmto)
        end
        function dmto = getDebugMessagesTurnedOn()
            try
                val = getpref(ClusterInfo.Group,'DebugMessagesTurnedOn');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = false;
            end
            dmto = val;
        end

        function setDiskSpace(ds)
            if nargin==0 || ischar(ds)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Disk space must be a character string.')
            end
            setpref(ClusterInfo.Group,'DiskSpace',ds)
        end
        function ds = getDiskSpace()
            try
                val = getpref(ClusterInfo.Group,'DiskSpace');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            ds = val;
        end

        function setEmailAddress(ea)
            if nargin==0 || ischar(ea)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Email address must be a character string.')
            end
            exp = '[A-Z0-9._%+-]+@(?:[A-Z0-9-]+\.)+[A-Z]{2,4}';
            if isempty(regexpi(ea,exp)) && isempty(ea)==false
                warning('distcomp:clusterinfo:EmailAddress', ...
                    ['''' ea ''' doesn''t appear to be a valid email address.'])
            end
            setpref(ClusterInfo.Group,'EmailAddress',ea)
        end
        function ea = getEmailAddress()
            try
                val = getpref(ClusterInfo.Group,'EmailAddress');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            ea = val;
        end

        function setGpusPerNode(gpn)
            if nargin==0 || isnumeric(gpn)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'GPUs per node must be an integer.')
            end
            setpref(ClusterInfo.Group,'GpusPerNode',gpn)
        end
        function gpn = getGpusPerNode()
            try
                val = getpref(ClusterInfo.Group,'GpusPerNode');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            gpn = val;
        end

        function setMemUsage(mu)
            if nargin==0 || ischar(mu)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Memory usage must be a character string.')
            end
            setpref(ClusterInfo.Group,'MemUsage',mu)
        end
        function mu = getMemUsage()
            try
                val = getpref(ClusterInfo.Group,'MemUsage');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            mu = val;
        end

        function setNameSpace(ns)
            if nargin==0 || ischar(ns)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Name space must be a character string.')
            end
            setpref(ClusterInfo.Group,'NameSpace',ns)
        end
        function ns = getNameSpace()
            try
                val = getpref(ClusterInfo.Group,'NameSpace');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            ns = val;
        end

        function setPrivateKeyFile(pkf)
            if nargin==0 || ischar(pkf)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Private key file must be a character string.')
            end
            setpref(ClusterInfo.Group,'PrivateKeyFile',pkf)
        end
        function pkf = getPrivateKeyFile()
            try
                val = getpref(ClusterInfo.Group,'PrivateKeyFile');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            pkf = val;
        end

        function setPrivateKeyFileHasPassPhrase(pkhpp)
            if nargin==0 || islogical(pkhpp)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Flag must be boolean.')
            end
            setpref(ClusterInfo.Group,'PrivateKeyFileHasPassPhrase',pkhpp)
        end
        function pkhpp = getPrivateKeyFileHasPassPhrase()
            try
                val = getpref(ClusterInfo.Group,'PrivateKeyFileHasPassPhrase');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = true;
            end
            pkhpp = val;
        end

        function setProcsPerNode(ppn)
            if nargin==0 || isnumeric(ppn)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Procs per node must be an integer.')
            end
            setpref(ClusterInfo.Group,'ProcsPerNode',ppn)
        end
        function ppn = getProcsPerNode()
            try
                val = getpref(ClusterInfo.Group,'ProcsPerNode');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            ppn = val;
        end

        function setProjectName(pn)
            if nargin==0 || ischar(pn)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Project name must be a character string.')
            end
            setpref(ClusterInfo.Group,'ProjectName',pn)
        end
        function pn = getProjectName()
            try
                val = getpref(ClusterInfo.Group,'ProjectName');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            pn = val;
        end

        function setQueueName(qn)
            if nargin==0 || ischar(qn)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Queue name must be a character string.')
            end
            setpref(ClusterInfo.Group,'QueueName',qn)
        end
        function ch = getQueueName()
            try
                val = getpref(ClusterInfo.Group,'QueueName');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            ch = val;
        end

        function setRequireExclusiveNode(ren)
            if nargin==0 || islogical(ren)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Flag must be boolean.')
            end
            setpref(ClusterInfo.Group,'RequireExclusiveNode',ren)
        end
        function ren = getRequireExclusiveNode()
            try
                val = getpref(ClusterInfo.Group,'RequireExclusiveNode');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = false;
            end
            ren = val;
        end

        function setReservation(r)
            if nargin==0 || ischar(r)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Reservation must be a character string.')
            end
            setpref(ClusterInfo.Group,'Reservation',r)
        end
        function r = getReservation()
            try
                val = getpref(ClusterInfo.Group,'Reservation');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            r = val;
        end

        function setSshPort(sp)
            if nargin==0 || isnumeric(sp)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'SSH port must be an integer.')
            end
            setpref(ClusterInfo.Group,'SshPort',sp)
        end
        function sp = getSshPort()
            try
                val = getpref(ClusterInfo.Group,'SshPort');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = [];
            end
            sp = uint32(val);
        end

        function setUseGpu(ug)
            if nargin==0 || islogical(ug)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Flag must be boolean.')
            end
            if ug==true
                % MW: Is this correct?  Should we necessarily clear out the
                % queue?
                ClusterInfo.setQueueName('')
            end
            setpref(ClusterInfo.Group,'UseGpu',ug)
        end
        function ug = getUseGpu()
            try
                val = getpref(ClusterInfo.Group,'UseGpu');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = false;
            end
            ug = val;
        end

        function setUserDefinedOptions(udo)
            if nargin==0 || ischar(udo)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Options must be a character string.')
            end
            setpref(ClusterInfo.Group,'UserDefinedOptions',udo)
        end
        function udo = getUserDefinedOptions()
            try
                val = getpref(ClusterInfo.Group,'UserDefinedOptions');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            udo = val;
        end

        function setUserNameOnCluster(unoc)
            if nargin==0 || ischar(unoc)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Username must be a character string.')
            end
            setpref(ClusterInfo.Group,'UserName',unoc)
        end
        function unoc = getUserNameOnCluster()
            try
                val = getpref(ClusterInfo.Group,'UserName');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            unoc = val;
        end

        function setWallTime(wt)
            if nargin==0 || ischar(wt)==false
                error('distcomp:clusterinfo:InvalidType', ...
                    'Wall time must be a character string.')
            end
            setpref(ClusterInfo.Group,'WallTime',wt)
        end
        function wt = getWallTime()
            try
                val = getpref(ClusterInfo.Group,'WallTime');
            catch E %#ok<NASGU>
                % TODO: Should this throw an error
                val = '';
            end
            wt = val;
        end

    end

end
