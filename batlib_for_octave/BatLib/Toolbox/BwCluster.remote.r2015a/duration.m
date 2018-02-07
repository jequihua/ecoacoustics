function [sec, tsec] = duration(obj)
%DURATION Calculate the number of seconds from start time to finish time.
%   SEC = DURATION(OBJ) calculates total time (in seconds) to run a
%   distributed job or task.
%
%   [SEC, TSEC] = DURATION(OBJ) also calculates total time (in seconds) to
%   run a distributed job or task (TSEC), from submit time to finish time.
%
%   If no return argument is assigned, the duration is only displayed.  If
%   the job or task has not finished, use the current time to show duration
%   to this point.
%
%   DURATION only calculates time for tasks or jobs completed within 30
%   days.  Assumes submission host and execution host are in the same time
%   zone.
%
%   Examples
%   ========
%   time = duration(job);
%   duration(job.Task(1))
%    515 seconds

%   Copyright 2007-2015 The MathWorks, Inc.
%   Raymond S. Norris (raymond.norris@mathworks.com)

if isempty(strfind(class(obj),'job'))==true
    typ = 'Task';
else
    typ = 'Job';
end

% Get the create and finished times.  (Not sure if this is when create
% started or create finished).
ct = get(obj,'CreateTime');
st = get(obj,'StartTime');
ft = get(obj,'FinishTime');

if isempty(ct)
    % Job/Task is in a bad state if it doesn't have a creation time
    disp([typ ' has no creation time.'])
    sec = NaN;
    tsec = NaN;
    return
end

if isempty(st)
    % Job/Task has not started running yet.
    disp([typ ' has not started.'])
    sec = NaN;
    tsec = NaN;
    return
end

if isempty(ft)
    disp([typ ' is still running.  Using current time as finished time.'])
    % datestr() doesn't support a time zone, so just put in EST.  It
    % doesn't matter for our calculations.
    ft = [datestr(now,'ddd mmm dd HH:MM:SS') ' EST ' datestr(now,'yyyy')];
end

try
    % Remove the standard time marking (e.g. EST) and the day (e.g. Thu).
    sidx = strfind(ct,' ');
    ct([1:4 sidx(4):sidx(5)-1]) = [];
    st([1:4 sidx(4):sidx(5)-1]) = [];
    ft([1:4 sidx(4):sidx(5)-1]) = [];

    % Convert to a date number
    fmt = 'mmm dd HH:MM:SS yyyy';
    ct = datenum(ct,fmt);
    st = datenum(st,fmt);
    ft = datenum(ft,fmt);

    % Get time difference.  Let's assume that it doesn't go past 30 days.
    [~,~,d,h,mi,s] = datevec(ft-st);
    s = s + mi*60 + h*60*60 + d*24*60*60;

    % Get time difference.  Let's assume that it doesn't go past 30 days.
    [~,~,d,h,mi,ts] = datevec(ft-ct);
    ts = ts + mi*60 + h*60*60 + d*24*60*60;

catch le
    disp(le.message)
    s = NaN;
    ts = NaN;
end

if nargout==0
    if isnan(s), s=0; end
    disp([char(9) num2str(s) ' seconds'])
    %     if isnan(ts), ts=0; end
else
    sec = s;
    tsec = ts;
end

end
