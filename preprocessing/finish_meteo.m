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
% t_bef and t_aft are FORCED to the ERA5-Land convention rather than optimised per
% site. ERA5-Land accumulations run from 00 UTC and are de-accumulated to the
% preceding hour, so the documented offsets are t_bef = 1, t_aft = 0; the
% optimiser lands near 0.75/0.25, which is what the shipped US_xRM file carries.
% One value per PRODUCT, reused across sites -- an optimiser run per station would
% make the network inconsistent for no physical reason.
%
% ERA5-Land is UTC, so DeltaGMT must be 0 and Lon must be the true site longitude,
% or the solar geometry is wrong by hours.

if nargin < 4 || isempty(year_tag),      year_tag = '1985_2020';                  end
if nargin < 3 || isempty(partition_dir), partition_dir = fullfile('..','T&C','Diogo'); end
if nargin < 2 || isempty(out_dir),       out_dir = raw_dir;                       end

T_BEF = 0.75;   % see the note above
T_AFT = 0.25;

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

    try
        [~,~,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,~,t_bef,t_aft] = ...
            C_Automatic_Radiation_Partition(S.Date, S.Lat, S.Lon, S.Zbas, ...
                S.DeltaGMT, S.Pr, S.Tdew, S.Rsw, 0, T_BEF, T_AFT);
    catch ME
        fprintf('  ! %-10s partition failed: %s\n', site, ME.message);
        continue
    end

    S.SAD1 = SAD1(:); S.SAD2 = SAD2(:);
    S.SAB1 = SAB1(:); S.SAB2 = SAB2(:);
    S.PARB = PARB(:); S.PARD = PARD(:);
    S.N    = N(:);
    S.t_bef = t_bef;  S.t_aft = t_aft;

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
