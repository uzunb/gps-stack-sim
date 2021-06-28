%% GPSDATETIME
%   [T, D, W, R, LS] = GPSDATETIME(...,'LeapSecond',LS) returns GPS time, in terms of:
%    - Time Of Week T (optional),
%    - Day number D (optional), 
%    - Week W (optional),
%    - Rollover R (optional),
%    - Leap second LS (optional),
%   from:
%    - DATETIME function arguments,
%    - Leap second LS (optional, default value based on history if not specified).
%
%  Example: display of current GPS time (no outputs)
%   GPSDATETIME('now')
%
%  Example: elements from current GPS time
%   [T, D, W, R, LS] = GPSDATETIME('now')
%
%  Example: GPS time at 06/04/2019 00:00:00 - leap second (second rollover)
%   [T, D, W, R] = GPSDATETIME('06/04/2019 23:59:42',...
%                              'format',     'dd/MM/uuuu HH:mm:ss',...
%                              'leapsecond', 18)
%
%  See also DATETIME.

function varargout = GPSdatetime(varargin)

% GPS reference date
ReferenceDate = ...
    datetime('06/01/1980',...
             'InputFormat', 'dd/MM/yyyy',...
             'TimeZone',    'UTC');
         
% Leap seconds         
ls = [];
i = find(strcmpi(varargin,'leapsecond'));
if ~isempty(i)
    ls = varargin{i+1};
    varargin(i:i+1) = [];
end

% Current time
T = feval(@datetime,varargin{:});
T.TimeZone = 'UTC';

% Default leap second
if isempty(ls)    
    
    % Leap seconds dates (01/01/1980: LS = 0)
    LeapSecondsDates = ...
        datetime({'01/01/1980',...
                  '01/07/1981',...
                  '01/07/1982',...
                  '01/07/1983',...
                  '01/07/1985',...
                  '01/01/1988',...
                  '01/01/1990',...
                  '01/01/1991',...
                  '01/07/1992',...
                  '01/07/1993',...
                  '01/07/1994',...
                  '01/01/1996',...
                  '01/07/1997',...
                  '01/01/1999',...
                  '01/01/2006',...
                  '01/01/2009',...
                  '01/07/2012',...
                  '01/07/2015',...
                  '01/01/2017'},...
                  'InputFormat', 'dd/MM/yyyy',...
                  'TimeZone',    'UTC');
    
    % Leap seconds
    ls = find(ge(T,LeapSecondsDates),1,'last')-1;
    
    % Control
    dt = T-LeapSecondsDates(end);
    if isempty(ls)
        error('gpsdatetime: invalid GPS date (GPS date shall be posterior to 01/01/1980).');
    elseif gt(dt,diff(LeapSecondsDates(end-1:end)))
        y = years(dt);
        d = mod(days(dt),365);
        S = {'','s'};
        s = @(n)S{gt(n,1)+1};
        warning('gpsdatetime: possibly missed leap seconds (the last leap second occurred %.0f year%c %.0f day%c ago).',y,s(y),d,s(d));
    end

end

% Adding of leap seconds
T = T+ls/86400;

% Current date
DV = datevec(T);

% Day number
d = day(T,'dayofweek');

% GPS time of week [s]
t = 86400*(d-1)+[3600 60 1]*DV(4:6)';

% Elapsed time from GPS reference date [day]
et = days(T-ReferenceDate);

% Rollover number
r = floor(et/7/1024);

% Week number
w = floor(et/7-r*1024);

% Outputs
switch nargout
    case 0
        fprintf(1,' TOW[s] | DAY |  WN  | RO | LS \n'); 
        fprintf(1,'-------------------------------\n');
        fprintf(1,' %06.0f |  %u  | %04u |  %u | %02u\n',t,d,w,r,ls);
    otherwise
        varargout = {t,d,w,r,ls};
        varargout = varargout(1:nargout);
end

end
