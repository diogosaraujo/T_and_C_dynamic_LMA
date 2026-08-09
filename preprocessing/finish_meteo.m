function finish_meteo(raw_dir, out_dir, partition_dir, year_tag)
%FINISH_METEO  Second stage of the T&C forcing build.
%
%   finish_meteo(raw_dir, out_dir, partition_dir, year_tag)
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

T_BEF = 1.0;    % de-accumulated ERA5-Land hour H covers (H-1, H]
T_AFT = 0.0;

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

    % Run one CALENDAR MONTH at a time. Something inside the partition scales as
    % N^2: at N = 315576 (36 years hourly) MATLAB refused a 315576 x 315576 array
    % outright (742 GB, job 35673), and at N = 8784 (one year) each such matrix is
    % still 617 MB, which OOM-killed a 24 GB job once MATLAB was willing to try
    % (35676, exit 137). A month is at most 744 hours, so the same term is ~4 MB.
    %
    % Chunking is safe because the calculation is per-timestep throughout: solar
    % geometry, Gueymard clear sky, the clearness index and the N = 1 when Pr > 0
    % rule all act on one hour at a time. There is no cumsum, filter or smoothing
    % across the series, and t_bef/t_aft are forced, so the calibration loop that
    % would need a long record (NI = min(8760,NT)) never runs.
    dts = datetime(S.Date, 'ConvertFrom', 'datenum');
    ym  = year(dts)*100 + month(dts);
    uy  = unique(ym(:)).';
    n   = numel(S.Date);
    SAD1 = zeros(n,1); SAD2 = zeros(n,1); SAB1 = zeros(n,1); SAB2 = zeros(n,1);
    PARB = zeros(n,1); PARD = zeros(n,1); N = zeros(n,1);
    failed = false;
    nchunk = 0;
    for y = uy
        ix = find(ym == y);
        try
            [~,~,d1,d2,b1,b2,pb,pd,nn,~,t_bef,t_aft] = ...
                C_Automatic_Radiation_Partition(S.Date(ix), S.Lat, S.Lon, S.Zbas, ...
                    S.DeltaGMT, S.Pr(ix), S.Tdew(ix), S.Rsw(ix), 0, T_BEF, T_AFT);
        catch ME
            fprintf('  ! %-10s %06d: partition failed: %s\n', site, y, ME.message);
            failed = true;
            break
        end
        SAD1(ix)=d1(:); SAD2(ix)=d2(:); SAB1(ix)=b1(:); SAB2(ix)=b2(:);
        PARB(ix)=pb(:); PARD(ix)=pd(:); N(ix)=nn(:);
        nchunk = nchunk + 1;
        if mod(nchunk, 60) == 0
            fprintf('      %s: %d/%d months\n', site, nchunk, numel(uy));
        end
        clear d1 d2 b1 b2 pb pd nn
    end
    if failed
        continue
    end

    S.SAD1 = SAD1(:); S.SAD2 = SAD2(:);
    S.SAB1 = SAB1(:); S.SAB2 = SAB2(:);
    S.PARB = PARB(:); S.PARD = PARD(:);
    S.N    = N(:);
    S.t_bef = t_bef;  S.t_aft = t_aft;
    if abs(t_bef - T_BEF) > 1e-9 || abs(t_aft - T_AFT) > 1e-9
        fprintf('  ! %-10s partition returned t_bef/t_aft = %g/%g, not the forced %g/%g\n', ...
            site, t_bef, t_aft, T_BEF, T_AFT);
    end

    % The bands must sum to the total shortwave, or light is being invented or
    % lost before it reaches the canopy.
    band_sum = S.SAB1 + S.SAB2 + S.SAD1 + S.SAD2;
    resid = max(abs(band_sum - S.Rsw(:)));
    if resid > 1e-6 * max(1, max(S.Rsw))
        fprintf('  ! %-10s bands do not close: max |SAB+SAD - Rsw| = %.3g W/m2\n', ...
            site, resid);
    end

    outfile = fullfile(out_dir, sprintf('Meteo_%s_%s.mat', site, year_tag));
    save(outfile, '-struct', 'S', '-v7.3');
    nok = nok + 1;
    fprintf('  %-10s %7d h  N %.2f-%.2f  PARBmax %6.1f  bands close to %.1e  -> %s\n', ...
        site, numel(S.Date), min(N), max(N), max(PARB), resid, ...
        sprintf('Meteo_%s_%s.mat', site, year_tag));
end

fprintf('\n%d/%d written to %s\n', nok, numel(files), out_dir);
rmpath(partition_dir);
end
