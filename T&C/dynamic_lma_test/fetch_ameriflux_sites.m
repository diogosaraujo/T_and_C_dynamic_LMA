%% fetch_ameriflux_sites.m
% -------------------------------------------------------------------------
% Pulls AmeriFlux site metadata (coordinates, IGBP, tower years, URL)
% directly from the AmeriFlux data API and writes AmeriFlux_site_list.csv,
% formatted for pair_ecoregion_ameriflux.m.
% No R, no Rtools, no manual download needed -- just internet access.
%
% NOTE ON YEARS: DATA_START/DATA_END written here are the tower OPERATIONAL
% years (TOWER_BEGAN / TOWER_END) -- a first-pass proxy for record length.
% They are NOT necessarily the years with submitted BASE flux data. Good
% enough to get the pairing running; refine with true BASE data-availability
% later if you need exact coverage for the >=10-year rule.
% -------------------------------------------------------------------------
clear; clc;

url = 'https://amfcdn.lbl.gov/api/v1/site_display/AmeriFlux';
D   = webread(url, weboptions('ContentType','json','Timeout',90));

% API returns a JSON array of site objects -> struct array OR cell of structs
if isstruct(D),      N = numel(D); getEl = @(i) D(i);
elseif iscell(D),    N = numel(D); getEl = @(i) D{i};
else, error('Unexpected JSON structure returned by the API.');
end

SITE_ID=strings(N,1); SITE_NAME=strings(N,1); IGBP=strings(N,1);
LAT=nan(N,1); LON=nan(N,1); ELEV=nan(N,1);
DSTART=strings(N,1); DEND=strings(N,1); URLA=strings(N,1);

for i = 1:N
    s = getEl(i);
    SITE_ID(i)   = gstr(s,'SITE_ID');
    SITE_NAME(i) = gstr(s,'SITE_NAME');
    IGBP(i)      = gstr(s,'IGBP');
    DSTART(i)    = gstr(s,'TOWER_BEGAN');
    DEND(i)      = gstr(s,'TOWER_END');
    URLA(i)      = gstr(s,'URL_AMERIFLUX');

    % coordinates live in the nested GRP_LOCATION group
    g = [];
    if isfield(s,'GRP_LOCATION') && ~isempty(s.GRP_LOCATION)
        g = s.GRP_LOCATION;
        if iscell(g), g = g{1}; elseif isstruct(g) && numel(g)>1, g = g(1); end
    end
    if ~isempty(g)
        LAT(i)  = gnum(g,'LOCATION_LAT');
        LON(i)  = gnum(g,'LOCATION_LONG');
        ELEV(i) = gnum(g,'LOCATION_ELEV');
    end
end

T = table(SITE_ID,SITE_NAME,LAT,LON,ELEV,IGBP,DSTART,DEND,URLA, ...
    'VariableNames',{'SITE_ID','SITE_NAME','LOCATION_LAT','LOCATION_LONG', ...
    'LOCATION_ELEV','IGBP','DATA_START','DATA_END','URL_AMERIFLUX'});
T = T(T.SITE_ID~="",:);
writetable(T,'AmeriFlux_site_list.csv');
fprintf('Wrote AmeriFlux_site_list.csv: %d sites (%d with coordinates)\n', ...
    height(T), sum(~isnan(T.LOCATION_LAT)));

%% helpers
function v = gstr(s,f)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
        v = string(s.(f)); v = v(1);
    else, v = "";
    end
end
function x = gnum(s,f)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
        x = str2double(string(s.(f)));
    else, x = NaN;
    end
end
