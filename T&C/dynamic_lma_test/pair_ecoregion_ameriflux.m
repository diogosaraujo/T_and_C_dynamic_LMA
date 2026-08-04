%% pair_ecoregion_ameriflux.m
% -------------------------------------------------------------------------
% Lists EVERY matching AmeriFlux tower for each valid PLSR ecoregion x
% forest-type cell (deciduous / evergreen) and writes a CSV.
%
% For each valid (ecoregion, forest type):
%   - find ALL towers of the matching IGBP class located inside that EPA L3
%     ecoregion;
%   - output one row per tower, ranked by record length (Rank 1 = longest);
%   - flag each as 'ge10yr' or 'under_10yr' (relative to MIN_YEARS);
%   - if there is NO matching tower, write a single 'not_available' row.
%
% REQUIREMENTS: MATLAB Mapping Toolbox (shaperead, projcrs, projfwd, inpolygon).
%
% INPUTS (edit USER SETTINGS):
%   shpFile : EPA L3 ecoregion shapefile with ECO_IDX, US_L3CODE, US_L3NAME,
%             and skill fields TQ2_DEC / TQ2_EVG (Albers CRS).
%   amfFile : AmeriFlux site list CSV (run fetch_ameriflux_sites.m first).
%   outFile : output CSV path.
%
% OUTPUT columns:
%   ECO_IDX, US_L3CODE, US_L3NAME, ForestType, PLSR_TemporalQ2,
%   nStations, Rank, StationID, StationName, Lat, Lon, IGBP,
%   DataStart, DataEnd, nYears, Status, StationURL
%   Status in {'ge10yr','under_10yr','not_available'}
% -------------------------------------------------------------------------
clear; clc;

%% ===================== USER SETTINGS =====================
shpFile = 'C:\Users\dalca\Documents\Droughts\New Analysis\ecoregions\PLSR_skill_metric_maps\us_eco_l3_LMA_plsr_skill_metrics.shp';
amfFile = 'AmeriFlux_site_list.csv';
outFile = 'ecoregion_ameriflux_pairing.csv';

MIN_YEARS  = 10;         % threshold used only to flag ge10yr vs under_10yr
igbpDec    = {'DBF'};    % IGBP classes treated as "deciduous"  (add 'DNF','MF','WSA' to broaden)
igbpEvg    = {'ENF'};    % IGBP classes treated as "evergreen"  (add 'EBF','MF','WSA' to broaden)
ONLY_VALID = true;       % restrict to cells with a fitted PLSR model (TQ2 ~= -9999)
%% =========================================================

%% 1) Ecoregion shapefile (USA Contiguous Albers Equal Area Conic USGS)
S = shaperead(shpFile);
try
    p = projcrs(102039,'Authority','ESRI');
catch
    p = projcrs(fileread(strrep(shpFile,'.shp','.prj')));
end

%% 2) AmeriFlux site list (flexible column detection)
A   = readtable(amfFile,'TextType','string','VariableNamingRule','preserve');
sid = string(pickcol(A,{'SITE_ID','Site ID','SiteID'}));
snm = string(pickcol(A,{'SITE_NAME','Name','Site Name','SiteName'}));
igbp= upper(strtrim(string(pickcol(A,{'IGBP','Vegetation Abbreviation (IGBP)','Vegetation Abbreviation'}))));
lat = double(pickcol(A,{'LOCATION_LAT','Latitude (degrees)','Latitude','Lat'}));
lon = double(pickcol(A,{'LOCATION_LONG','Longitude (degrees)','Longitude','Long','Lon'}));
ds  = pickcol_opt(A,{'DATA_START','Data Start','AmeriFlux Data Start'});
de  = pickcol_opt(A,{'DATA_END','Data End','AmeriFlux Data End'});
if isempty(ds) || isempty(de)
    avail = pickcol_opt(A,{'Data Availability','AmeriFlux BASE Data Availability','Years','DATA_AVAILABILITY'});
    if isempty(avail)
        warning('No data-year columns found; nYears will be NaN.');
        ds = nan(numel(sid),1); de = nan(numel(sid),1);
    else
        [ds,de] = parseAvail(avail);
    end
else
    ds = toYear(ds); de = toYear(de);
end
nyr = de - ds + 1;
url = "https://ameriflux.lbl.gov/sites/siteinfo/" + sid;

%% 3) Assign each tower to an ecoregion (ECO_IDX) by point-in-polygon
[tx,ty] = projfwd(p, lat, lon);
ecoOfTower = nan(numel(sid),1);
for i = 1:numel(sid)
    if isnan(tx(i)); continue; end
    for k = 1:numel(S)
        if inpolygon(tx(i),ty(i), S(k).X, S(k).Y)
            ecoOfTower(i) = S(k).ECO_IDX;  break
        end
    end
end

%% 4) One row per ecoregion: code, name, per-type skill
[uEco,ia] = unique([S.ECO_IDX]');
L3code = string({S(ia).US_L3CODE})';
L3name = string({S(ia).US_L3NAME})';
tq2dec = arrayfun(@(s) s.TQ2_DEC, S(ia))';
tq2evg = arrayfun(@(s) s.TQ2_EVG, S(ia))';

%% 5) List ALL matching towers per valid (ecoregion, type) cell
rows  = cell(0,17);
types = {'deciduous','evergreen'};
for t = 1:2
    ftype = types{t};
    if t==1, igbpSet = string(igbpDec); tq2 = tq2dec;
    else,    igbpSet = string(igbpEvg); tq2 = tq2evg; end
    for e = 1:numel(uEco)
        eco = uEco(e);  q2 = tq2(e);
        isValid = (q2 ~= -9999) && ~isnan(q2);
        if ONLY_VALID && ~isValid; continue; end

        cand = find(ecoOfTower==eco & ismember(igbp,igbpSet));
        if isempty(cand)
            rows(end+1,:) = {eco,L3code(e),L3name(e),ftype,q2, ...
                0,NaN,"","",NaN,NaN,"",NaN,NaN,NaN,"not_available",""}; %#ok<AGROW>
            continue
        end
        ny = nyr(cand); ny(isnan(ny)) = -1;      % rank longest-first
        [~,ord] = sort(ny,'descend');
        cand = cand(ord);  nc = numel(cand);
        for r = 1:nc
            k = cand(r);  yrs = nyr(k);
            if ~isnan(yrs) && yrs >= MIN_YEARS, st = "ge10yr"; else, st = "under_10yr"; end
            rows(end+1,:) = {eco,L3code(e),L3name(e),ftype,q2, ...
                nc,r,sid(k),snm(k),lat(k),lon(k),igbp(k), ...
                ds(k),de(k),yrs,st,url(k)}; %#ok<AGROW>
        end
    end
end

%% 6) Write CSV + summary
T = cell2table(rows,'VariableNames',{'ECO_IDX','US_L3CODE','US_L3NAME', ...
    'ForestType','PLSR_TemporalQ2','nStations','Rank','StationID','StationName', ...
    'Lat','Lon','IGBP','DataStart','DataEnd','nYears','Status','StationURL'});
T = sortrows(T,{'ForestType','ECO_IDX','Rank'});
writetable(T,outFile);

covered = T(T.Status~="not_available",:);
key = string(covered.ECO_IDX) + "_" + covered.ForestType;
fprintf('Wrote %s\n', outFile);
fprintf('  cells covered: %d | cells not_available: %d | total station-matches: %d (%d with >=%d yr)\n', ...
    numel(unique(key)), sum(T.Status=="not_available"), height(covered), ...
    sum(covered.Status=="ge10yr"), MIN_YEARS);

%% ===================== helpers =====================
function col = pickcol(T,names)
    col = pickcol_opt(T,names);
    if isempty(col); error('Required column not found. Tried: %s', strjoin(names,', ')); end
end
function col = pickcol_opt(T,names)
    vn = string(T.Properties.VariableNames); col = [];
    for n = string(names)
        h = find(strcmpi(vn,n),1);
        if isempty(h); h = find(contains(lower(vn),lower(n)),1); end
        if ~isempty(h); col = T.(char(vn(h))); return; end
    end
end
function y = toYear(v)
    if isnumeric(v); y = double(v); return; end
    s = string(v); y = nan(numel(s),1);
    for i = 1:numel(s)
        m = regexp(s(i),'\d{4}','match','once');
        if strlength(m)>0; y(i) = double(m); end
    end
end
function [ds,de] = parseAvail(v)
    s = string(v); n = numel(s); ds = nan(n,1); de = nan(n,1);
    for i = 1:n
        yy = double(regexp(s(i),'\d{4}','match'));
        if ~isempty(yy); ds(i) = min(yy); de(i) = max(yy); end
    end
end
