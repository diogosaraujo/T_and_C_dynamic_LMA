function finish_meteo(raw_dir, out_dir, partition_dir, year_tag, t_bef_in, t_aft_in)
%FINISH_METEO  Second stage of the T&C forcing build.
%
%   finish_meteo(raw_dir, out_dir, partition_dir, year_tag)
%
% DESTINATION. Each finished file goes to the 'dest_dir' stamped inside its own
% raw file by the Python builder -- model_run/<ST>/<scenario>/<GCM>/ for the GCM
% path, model_run/<ST>/era5_land/ for ERA5-Land -- so the forcing sits directly
% above the fixed_lma/dyn_lma pair that reads it and model_run needs nothing from
% input_data. out_dir is only the fallback for raw files written before that field
% existed. raw_dir stays a staging area: raw files are intermediates, consumed
% here and safe to delete afterwards.
%
% build_meteo_input.py writes Meteo_<ST>_raw.mat with everything that is a plain
% unit conversion. This adds the fields that are not: SAB1/SAB2/SAD1/SAD2,
% PARB/PARD and the cloud fraction N, by calling the existing, tested
% C_Automatic_Radiation_Partition. That routine is several hundred lines of solar
% geometry, Gueymard clear-sky and Slingo cloud physics -- reused rather than
% ported, because a Python port of comparable code shipped two bugs that only a
% verification harness caught.
%
% t_bef/t_aft are DERIVED from the product definition, not optimised -- but ONLY
% from the radiation variable, because ERA5-Land does not treat its variables
% alike:
%
%   ACCUMULATED from 00 UTC, reset daily   ssrd, tp
%       de-accumulated here, so the value at hour H is the mean over (H-1, H]
%   INSTANTANEOUS at the timestamp         t2m, d2m, sp, u10, v10
%       the value AT hour H, not an interval mean
%
% t_bef/t_aft is the window over which solar altitude is averaged to match the
% RADIATION timestamp, so it follows ssrd and nothing else:
%
%       t_bef = 1,  t_aft = 0
%
% Of the partition's own inputs, Pr comes from tp and shares that convention
% (it only enters through the N = 1 when Pr > 0 rule), while Tdew comes from d2m
% and is instantaneous. Tdew feeds the clear-sky water-vapour attenuation, where a
% half-hour offset is second order -- worth stating, not worth correcting.
%
% The same half-hour offset applies to Ta, ea and Ws in the forcing at large:
% instantaneous values at H are paired with interval-mean Rsw and Pr for (H-1, H].
% Averaging consecutive instantaneous values onto the interval would remove it; it
% is left uncorrected and stated, because it is small and because introducing a
% smoothing step would change the forcing in ways that are harder to audit.
%
% The optimiser inside C_Automatic_Radiation_Partition exists for datasets whose
% timestamp convention is undocumented. ERA5-Land's is documented, so fitting it
% would estimate a quantity we already know -- and it would return a slightly
% different answer per station, making the network inconsistent for no physical
% reason. (The shipped Meteo_US_xRM file carries 0.75/0.25, the optimiser's
% answer; that is what is being replaced here, deliberately.)
%
% This is a per-PRODUCT constant. The GCM path, whose radiation is disaggregated
% from daily to hourly with an imposed convention, takes 0/1 instead.
%
% ERA5-Land is UTC, so DeltaGMT must be 0 and Lon must be the true site longitude,
% or the solar geometry is wrong by hours.

if nargin < 4 || isempty(year_tag),      year_tag = '1985_2020';                  end
if nargin < 3 || isempty(partition_dir), partition_dir = fullfile('..','T&C','Diogo'); end
if nargin < 2 || isempty(out_dir),       out_dir = raw_dir;                       end

% Defaults are the ERA5-Land convention: de-accumulated hour H covers (H-1, H].
% The GCM path passes 0/1 instead, because build_gcm_meteo.py CONSTRUCTS the
% hourly series and places the value at hour H over (H, H+1]. A constructed
% series has no convention of its own to discover, so it is imposed at both ends
% and must agree -- getting this backwards shifts the solar geometry by an hour
% and shows up as a systematic Rn bias, not as an error.
if nargin < 6 || isempty(t_aft_in), t_aft_in = 0.0; end
if nargin < 5 || isempty(t_bef_in), t_bef_in = 1.0; end
T_BEF = t_bef_in;
T_AFT = t_aft_in;
why = {'imposed (constructed hourly series)', 'ERA5-Land, de-accumulated'};
fprintf('t_bef/t_aft = %g/%g -- %s\n', T_BEF, T_AFT, ...
    why{1 + double(T_BEF == 1.0 && T_AFT == 0.0)});

addpath(partition_dir);
if ~exist('C_Automatic_Radiation_Partition', 'file')
    error('finish_meteo:noPartition', ...
        'C_Automatic_Radiation_Partition not found on the path (looked in %s)', ...
        partition_dir);
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

files = dir(fullfile(raw_dir, 'Meteo_*_raw.mat'));
if isempty(files)
    error('finish_meteo:noInput', 'no Meteo_*_raw.mat in %s', raw_dir);
end
fprintf('%d raw file(s) in %s\n\n', numel(files), raw_dir);

nok = 0;
dests = {};
for k = 1:numel(files)
    f = files(k).name;
    site = regexprep(f, '^Meteo_|_raw\.mat$', '');
    S = load(fullfile(raw_dir, f));

    need = {'Date','Lat','Lon','Zbas','DeltaGMT','Pr','Tdew','Rsw'};
    miss = need(~isfield(S, need));
    if ~isempty(miss)
        fprintf('  ! %-10s missing %s -- skipped\n', site, strjoin(miss, ', '));
        continue
    end
    if S.DeltaGMT ~= 0
        fprintf('  ! %-10s DeltaGMT = %g, expected 0 for ERA5-Land (UTC)\n', ...
            site, S.DeltaGMT);
    end

    % Force COLUMN vectors before the partition sees them. This is the root cause
    % of the 742 GB request in job 35673: scipy.io.savemat writes a 1-D array as a
    % 1xN row, and the partition does LWP0 = zeros(size(Tdew)), so a row Tdew makes
    % LWP0 1xN, which then meets the Nx1 solar arrays and implicitly expands to
    % NxN. Cheap to enforce here as well as on the Python side, and it costs
    % nothing when the file already has columns.
    for fn = {'Date','Pr','Tdew','Rsw','Ta','Ws','ea','esat','Pre','Ca','Ds','U'}
        if isfield(S, fn{1})
            S.(fn{1}) = S.(fn{1})(:);
        end
    end

    % Validate BEFORE calling, rather than wrapping the call in try/catch. A
    % station whose ERA5 retrieval came back empty is a known-bad input and there
    % is nothing to learn from watching the partition fail on it -- job 35685 lost
    % 83 stations that way, dying inside Estimate_CloudCover where interp1 got an
    % empty X because every Nsim was NaN. Genuine errors, dimension mismatches
    % included, still abort the run and are still read.
    if ~any(isfinite(S.Rsw)) || ~any(S.Rsw > 0)
        fprintf('  ! %-10s Rsw has no usable values -- ERA5 retrieval empty; skipped\n', site);
        continue
    end
    unusable = {};
    for fn = {'Tdew','Pr','Date'}
        if ~any(isfinite(S.(fn{1}))), unusable{end+1} = fn{1}; end %#ok<AGROW>
    end
    if ~isempty(unusable)
        fprintf('  ! %-10s all-NaN in %s -- skipped\n', site, strjoin(unusable, ', '));
        continue
    end

    % One call on the whole series, exactly as the shipped Meteo_US_xRM.mat was
    % produced. No chunking and no per-chunk fallback: a dimension error here
    % should stop the run and be read, not be papered over by a shorter window
    % that quietly changes what the t_bef/t_aft calibration block sees.
    %
    % The earlier 742 GB request (job 35673) was not a size problem at all. It
    % came from passing ROW vectors: the partition does LWP0 = zeros(size(Tdew)),
    % so a 1xN Tdew makes LWP0 1xN, which then meets the Nx1 solar-geometry arrays
    % and implicitly expands to NxN. Hence the column coercion above -- and note
    % the routine column-ifies Date itself (line 21) and Rsw/Pr (line 171), but
    % only AFTER the calibration block at line 94, which is the gap we fell into.
    [~,~,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,~,t_bef,t_aft] = ...
        C_Automatic_Radiation_Partition(S.Date, S.Lat, S.Lon, S.Zbas, ...
            S.DeltaGMT, S.Pr, S.Tdew, S.Rsw, 0, T_BEF, T_AFT);

    S.SAD1 = SAD1(:); S.SAD2 = SAD2(:);
    S.SAB1 = SAB1(:); S.SAB2 = SAB2(:);
    S.PARB = PARB(:); S.PARD = PARD(:);
    S.N    = N(:);
    S.t_bef = t_bef;  S.t_aft = t_aft;
    if abs(t_bef - T_BEF) > 1e-9 || abs(t_aft - T_AFT) > 1e-9
        fprintf('  ! %-10s partition returned t_bef/t_aft = %g/%g, not the forced %g/%g\n', ...
            site, t_bef, t_aft, T_BEF, T_AFT);
    end

    % Band closure. These never close exactly, and the benchmark is the shipped
    % Meteo_US_xRM.mat, which this routine produced: max residual 77.2 W/m2, 3.27%
    % of hours above 1 W/m2, and every residual over 10 W/m2 occurring where Rsw is
    % only 10-77 W/m2. So the mismatch lives at low light, where the partition
    % cannot reconcile a small measured flux with its clear-sky decomposition.
    %
    % Demanding exact closure therefore flags every run. What would actually be
    % wrong is a large residual at HIGH irradiance, so that is what is tested.
    band_sum = S.SAB1(:) + S.SAB2(:) + S.SAD1(:) + S.SAD2(:);
    r = abs(band_sum - S.Rsw(:));
    rswc = S.Rsw(:);
    resid = max(r);
    frac = 100 * mean(r > 1);
    bad = r > 10 & rswc > 200;
    if any(bad)
        fprintf('  ! %-10s %d hour(s) with |SAB+SAD - Rsw| > 10 W/m2 at Rsw > 200: max %.1f\n', ...
            site, sum(bad), max(r(bad)));
    end

    % Restore the orientation the shipped Meteo_US_xRM.mat uses: the
    % meteorological inputs as 1xN rows, Date, Ca and the partition outputs as Nx1
    % columns. The partition needed columns; T&C should see what it always has.
    for fn = {'Pr','Tdew','Rsw','Ta','Ws','ea','esat','Pre','Ds'}
        if isfield(S, fn{1}), S.(fn{1}) = reshape(S.(fn{1}), 1, []); end
    end
    for fn = {'Date','Ca','N','PARB','PARD','SAB1','SAB2','SAD1','SAD2','U'}
        if isfield(S, fn{1}), S.(fn{1}) = reshape(S.(fn{1}), [], 1); end
    end

    % WHERE THE FINISHED FILE GOES. The forcing belongs with the runs that read
    % it, one copy per (station, scenario, GCM) directly above the fixed_lma /
    % dyn_lma pair that shares it -- so model_run is self-contained and the two
    % arms read one file instead of holding a symlink each into input_data.
    %
    % The destination is carried IN the raw file, stamped by whichever Python
    % builder wrote it, rather than parsed out of the filename here. Parsing
    % cannot work: 'Meteo_US_Wrc_GFDL_ESM4_historical_raw.mat' splits on
    % underscores that belong to the station, the GCM and the scenario alike, and
    % no rule separates them. The builder already knows all three.
    %
    % A raw file with no dest_dir is an ERROR, not a reason to fall back to
    % out_dir. Silently writing the forcing somewhere other than the run tree
    % produces a file nothing reads, while the log says the station succeeded.
    % Every builder stamps this field; a file without one predates the move and
    % must be rebuilt, not guessed at.
    if ~isfield(S, 'dest_dir') || isempty(S.dest_dir)
        error('finish_meteo:noDestination', ...
              ['%s has no dest_dir. It was written before the forcing moved ' ...
               'into model_run; re-run stage 1 rather than partitioning it ' ...
               'into %s, where nothing will read it.'], f, out_dir);
    end
    dest = S.dest_dir;
    if ~ischar(dest), dest = char(dest); end
    dest = strtrim(reshape(dest, 1, []));
    S = rmfield(S, 'dest_dir');         % routing metadata, not forcing
    if ~exist(dest, 'dir'), mkdir(dest); end

    outfile = fullfile(dest, sprintf('Meteo_%s_%s.mat', site, year_tag));
    save(outfile, '-struct', 'S', '-v7.3');
    nok = nok + 1;
    dests{end+1} = dest;   %#ok<AGROW> -- one per file, numel(files) at most
    fprintf(['  %-10s %7d h  N %.2f-%.2f  PARBmax %6.1f  Rswmax %6.1f  ' ...
             'band resid max %5.1f W/m2 (%.2f%% of hours > 1; reference ' ...
             'US_xRM: 77.2, 3.27%%)\n'], ...
        site, numel(S.Date), min(N), max(N), max(PARB), max(S.Rsw), resid, frac);
end

% Destinations are per file now, so name them rather than echoing one out_dir --
% "written to <the one place>" would be a lie the moment routing is in use.
u = unique(dests);
if numel(u) == 1
    fprintf('\n%d/%d written to %s\n', nok, numel(files), u{1});
else
    fprintf('\n%d/%d written across %d destination folder(s), e.g. %s\n', ...
        nok, numel(files), numel(u), u{1});
end
rmpath(partition_dir);
end
