function status = GRAPH_COMPARE(station_dir, out_dir, dpi)
%GRAPH_COMPARE  Fixed-LMA vs dynamic-LMA comparison figures for one station.
%
%   GRAPH_COMPARE('/.../model_run/US-HBK')                 -> <station>/figures_compare
%   GRAPH_COMPARE(station_dir, out_dir, dpi)
%
% Reads BOTH arms' RES files and redraws the GRAPH_MOD figures that carry the
% fixed-vs-dynamic signal, with the two arms overlaid or set side by side.
% GRAPH_MOD itself is untouched -- it plots one run at a time by design, and the
% treatment is only visible in the difference.
%
% Which figure becomes what, following the original layout as closely as possible:
%   table   every parameter in the generated MOD_PARAM except Sl_H/LMA (which is
%           the treatment and differs between arms by construction)
%   19c     total weighted LAI, both arms overlaid
%   21c     annual GPP/NPP/ANPP -- ORIGINAL COLOURS kept (k/r/g), arms separated
%           by line style (fixed solid, dynamic dashed)
%   22c     the six water/rain-use efficiencies -- ORIGINAL LINE STYLE kept
%           (dashed), arms separated by colour (fixed black, dynamic red)
%   101c    assimilation and foliage respiration
%   103c    LAI, GPP/NPP and leaf age
%   104c    phenology state
%   105c    carbon pools B(1:4)
%           ...101c/103c/104c/105c replace the original High-vs-Low vegetation
%           contrast with fixed-vs-dynamic. Low vegetation is off at these forest
%           sites (Ccrown = 1, High only), so those panels were empty anyway.
%   106c    soil temperature and beta factor, 1x2 panels (fixed | dynamic)
%   1001c-1004c  monthly water balance, 1x2 panels -- these are stacked/filled
%           area plots that cannot be legibly overlaid, hence side by side.
%
% Returns 0 if every figure was produced, 1 otherwise.

if nargin < 2 || isempty(out_dir)
    out_dir = fullfile(station_dir, 'figures_compare');
end
if nargin < 3 || isempty(dpi), dpi = 150; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

FIXED = fullfile(station_dir, 'era5_land', 'fixed_lma');
DYN   = fullfile(station_dir, 'era5_land', 'dyn_lma');
[~, station] = fileparts(station_dir);

fprintf('station : %s\n  fixed : %s\n  dyn   : %s\n  out   : %s\n', ...
    station, FIXED, DYN, out_dir);

F = load_arm(FIXED);
D = load_arm(DYN);
fprintf('  loaded both arms (%d daily, %d hourly steps)\n', ...
    numel(F.LAI_H(:,1)), numel(F.An_H(:,1)));

% Colours: fixed is the reference, dynamic is the treatment.
CF = [0 0 0]; CD = [0.85 0.1 0.1];
set(0, 'DefaultFigureVisible', 'off');
nfail = 0; nfig = 0;

%% ---------------------------------------------------- parameter table
try
    nfig = nfig + 1;
    draw_table(fullfile(FIXED, ['MOD_PARAM_' matlab_name(station) '.m']), station);
    save_fig(out_dir, station, 'table_parameters', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'parameter table');
end

%% ---------------------------------------------------- 19c  total LAI
try
    nfig = nfig + 1;
    newfig();
    plot(F.NNd, F.LAItot, '-',  'Color', CF, 'LineWidth', 1.2); hold on; grid on;
    plot(D.NNd, D.LAItot, '--', 'Color', CD, 'LineWidth', 1.2);
    ylabel('[LAI]'); xlabel('Day');
    title(sprintf('%s  LAI - Leaf Area Index (total weighted)', station));
    legend('fixed LMA', 'dynamic LMA', 'Location', 'best');
    save_fig(out_dir, station, 'fig19c_LAI', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig19c');
end

%% ---------------------------------------------------- 21c  productivities
% Original colours retained; the arms differ by line style.
try
    nfig = nfig + 1;
    newfig();
    hold on; grid on;
    plot(F.Yrs_yr, F.GPP_yr,  '-k', 'LineWidth', 1.5);
    plot(F.Yrs_yr, F.NPP_yr,  '-r', 'LineWidth', 1.5);
    plot(F.Yrs_yr, F.ANPP_yr, '-g', 'LineWidth', 1.5);
    plot(D.Yrs_yr, D.GPP_yr,  '--k', 'LineWidth', 1.5);
    plot(D.Yrs_yr, D.NPP_yr,  '--r', 'LineWidth', 1.5);
    plot(D.Yrs_yr, D.ANPP_yr, '--g', 'LineWidth', 1.5);
    xlabel('Year'); ylabel('[gC/m^2]');
    title(sprintf('%s  Vegetation Productivities   (solid = fixed, dashed = dynamic)', station));
    legend('GPP fixed','NPP fixed','ANPP fixed','GPP dyn','NPP dyn','ANPP dyn', ...
        'Location','best','NumColumns',2);
    save_fig(out_dir, station, 'fig21c_productivity', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig21c');
end

%% ---------------------------------------------------- 22c  efficiencies
% Original dashed line style retained; the arms differ by colour.
try
    nfig = nfig + 1;
    newfig();
    E = {  @(a) a.GPP_yr./a.ET_yr,                'EWUE = GPP/ET',        '[gC/m^2 mm]'
           @(a) a.GPP_yr./a.T_yr,                 'WUE_T = GPP/T',        '[gC/m^2 mm]'
           @(a) 0.001*a.GPP_yr.*a.VPD_yr./a.ET_yr,'IWUE = GPP*VPD/ET',    '[gC/m^2 mm]'
           @(a) a.GPP_yr./a.T_yr,                 'WUE_{leaf} = Ag/T',    '[gC/m^2 mm]'
           @(a) a.ANPP_yr./a.Pr_yr,               'RAIN USE EFFICIENCY',  '[gC/m^2 mm]'
           @(a) a.ET_yr./(a.ET_yr+a.Lk_yr),       'HORTON INDEX',         '[-]' };
    for k = 1:6
        subplot(3,2,k); hold on; grid on;
        plot(F.Yrs_yr, E{k,1}(F), '--', 'Color', CF, 'LineWidth', 1.2);
        plot(D.Yrs_yr, E{k,1}(D), '--', 'Color', CD, 'LineWidth', 1.2);
        xlabel('Year'); ylabel(E{k,3}); title(E{k,2});
        if k == 1, legend('fixed','dynamic','Location','best'); end
    end
    save_fig(out_dir, station, 'fig22c_efficiency', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig22c');
end

%% ------------------------------- 101c assimilation (was High vs Low)
try
    nfig = nfig + 1;
    newfig();
    subplot(2,1,1); hold on; grid on;
    plot(F.NN, F.An_H, 'b', 'LineWidth', 1.0);
    plot(F.NN, F.Rdark_H, 'r', 'LineWidth', 1.0);
    ylabel('[\mu mol CO_2 / m^2 s]'); title('Assimilation Rate -- FIXED LMA');
    legend('Net Assimilation','Foliage Respiration','Location','best');
    subplot(2,1,2); hold on; grid on;
    plot(D.NN, D.An_H, 'b', 'LineWidth', 1.0);
    plot(D.NN, D.Rdark_H, 'r', 'LineWidth', 1.0);
    ylabel('[\mu mol CO_2 / m^2 s]'); xlabel('Hour');
    title('Assimilation Rate -- DYNAMIC LMA');
    save_fig(out_dir, station, 'fig101c_assimilation', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig101c');
end

%% ------------------------------- 103c LAI / GPP-NPP / leaf age
try
    nfig = nfig + 1;
    newfig();
    subplot(3,1,1); hold on; grid on;
    plot(F.NNd, F.LAI_H, '-', 'Color', CF, 'LineWidth', 1.2);
    plot(D.NNd, D.LAI_H, '--','Color', CD, 'LineWidth', 1.2);
    plot(F.NNd, F.LAIdead_H, ':', 'Color', CF, 'LineWidth', 1.0);
    plot(D.NNd, D.LAIdead_H, ':', 'Color', CD, 'LineWidth', 1.0);
    ylabel('[LAI]'); title('LAI (dotted = dead LAI)');
    legend('fixed','dynamic','Location','best');
    subplot(3,1,2); hold on; grid on;
    plot(F.NNd, F.NPP_H, '-', 'Color', CF, 'LineWidth', 1.2);
    plot(D.NNd, D.NPP_H, '--','Color', CD, 'LineWidth', 1.2);
    plot(F.NNd, F.NPP_H+F.RA_H, ':', 'Color', CF, 'LineWidth', 1.0);
    plot(D.NNd, D.NPP_H+D.RA_H, ':', 'Color', CD, 'LineWidth', 1.0);
    ylabel('[gC / m^2 d]'); title('GPP/NPP (dotted = GPP)');
    subplot(3,1,3); hold on; grid on;
    plot(F.NNd, F.AgeL_H, '-', 'Color', CF, 'LineWidth', 1.2);
    plot(D.NNd, D.AgeL_H, '--','Color', CD, 'LineWidth', 1.2);
    ylabel('[days]'); xlabel('Day'); title('Average Leaf Age');
    save_fig(out_dir, station, 'fig103c_LAI_NPP_age', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig103c');
end

%% ------------------------------- 104c phenology state
try
    nfig = nfig + 1;
    newfig();
    subplot(2,1,1); plot(F.NNd, F.PHE_S_H, 'Color', CF, 'LineWidth', 1.2); grid on;
    ylabel('[#]'); title('PHENOLOGY STATE -- FIXED LMA   (1 dormant, 2 max growth, 3 normal, 4 senescence)');
    ylim([0 5]);
    subplot(2,1,2); plot(D.NNd, D.PHE_S_H, 'Color', CD, 'LineWidth', 1.2); grid on;
    ylabel('[#]'); xlabel('Day'); title('PHENOLOGY STATE -- DYNAMIC LMA'); ylim([0 5]);
    save_fig(out_dir, station, 'fig104c_phenology', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig104c');
end

%% ------------------------------- 105c carbon pools
try
    nfig = nfig + 1;
    newfig();
    nm = {'Foliage','Sapwood','Fine Roots','Carbohydrate Reserve'};
    subplot(2,1,1); hold on; grid on;
    for k = 1:4, plot(F.NNd, F.B_H(:,k), 'LineWidth', 1.2); end
    title('Carbon Pool H_{VEG} -- FIXED LMA'); ylabel('gC/m^2');
    legend(nm, 'Location','best','NumColumns',2);
    subplot(2,1,2); hold on; grid on;
    for k = 1:4, plot(D.NNd, D.B_H(:,k), 'LineWidth', 1.2); end
    title('Carbon Pool H_{VEG} -- DYNAMIC LMA'); ylabel('gC/m^2'); xlabel('Days');
    save_fig(out_dir, station, 'fig105c_carbon_pools', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig105c');
end

%% ------------------------------- 106c soil T and beta, side by side
try
    nfig = nfig + 1;
    newfig();
    A = {F, D}; ttl = {'FIXED LMA','DYNAMIC LMA'};
    for k = 1:2
        subplot(1,2,k); hold on; grid on;
        yyaxis left;  plot(A{k}.NNd, A{k}.TdpI_H, 'LineWidth', 1.0); ylabel('Soil T 30-day mean [\circC]');
        yyaxis right; plot(A{k}.NNd, A{k}.Bfac_dayH, 'LineWidth', 1.0); ylabel('\beta factor [-]');
        xlabel('Day'); title(ttl{k});
    end
    save_fig(out_dir, station, 'fig106c_soilT_beta', dpi);
catch ME
    nfail = nfail + 1; warn(ME, 'fig106c');
end

%% ------------------------------- 1001c-1004c monthly water balance
mb = { '1001c_monthly_fluxes',   1
       '1002c_monthly_balance',  2
       '1003c_fraction_fluxes',  3
       '1004c_fraction_balance', 4 };
for m = 1:4
    try
        nfig = nfig + 1;
        newfig();
        A = {F, D}; ttl = {'FIXED LMA','DYNAMIC LMA'};
        for k = 1:2
            subplot(1,2,k);
            draw_monthly(A{k}, mb{m,2});
            title(ttl{k});
        end
        save_fig(out_dir, station, ['fig' mb{m,1}], dpi);
    catch ME
        nfail = nfail + 1; warn(ME, mb{m,1});
    end
end

close all
fprintf('\n%d of %d figure(s) written to %s\n', nfig-nfail, nfig, out_dir);
status = double(nfail > 0);
end

%% ===================================================================== helpers

function A = load_arm(d)
% Load only what the comparison needs. A full workspace is several GB and two
% arms at once would not fit comfortably; naming the variables keeps it small.
hits = dir(fullfile(d, 'RES_*.mat'));
if numel(hits) ~= 1
    error('GRAPH_COMPARE:res', 'expected exactly one RES_*.mat in %s, found %d', d, numel(hits));
end
v = {'LAI_H','LAI_L','LAIdead_H','NPP_H','NPP_L','RA_H','RA_L','ANPP_H','ANPP_L', ...
     'AgeL_H','PHE_S_H','B_H','TdpI_H','Bfac_dayH','An_H','Rdark_H', ...
     'T_H','T_L','EG','EIn_H','EIn_L','EIn_urb','EIn_rock','ELitter','ESN','ESN_In', ...
     'EICE','Lk','Rh','Rd','Qi_in','Qi_out','Pr','Pr_liq','Pr_sno','Ds', ...
     'Date','Datam','Ccrown','dth'};
A = load(fullfile(hits(1).folder, hits(1).name), v{:});
A.NNd = (1:size(A.LAI_H,1))';
A.NN  = (1:size(A.An_H,1))';
A.LAItot = (A.LAI_H + A.LAI_L) * A.Ccrown';
A.B_H = squeeze(A.B_H(:,1,:));            % one crown at these sites
A = aggregate(A);
end

function A = aggregate(A)
% Annual and monthly aggregates, replicated from GRAPH_MOD's switch_summary block
% so the comparison figures show exactly the same quantities.
Yrs = year(A.Date); r = 0;
for i = min(Yrs):max(Yrs)
    k = find(Yrs == i);
    if numel(k) > 350*24
        r = r + 1;
        A.Pr_yr(r) = sum(A.Pr(k));
        A.ET_yr(r) = sum(sum(A.T_H(k,:),2) + sum(A.T_L(k,:),2) + A.EG(k) + ...
            sum(A.EIn_H(k,:),2) + sum(A.EIn_L(k,:),2) + A.EIn_urb(k) + A.EIn_rock(k) + ...
            A.ESN(k) + A.ESN_In(k));
        A.Lk_yr(r) = sum(A.Lk(k));
        A.T_yr(r)  = sum(sum(A.T_H(k,:),2) + sum(A.T_L(k,:),2));
        A.VPD_yr(r)= mean(A.Ds(k));
    end
end
Yrs = year(A.Date(1):1:A.Date(end)+1);
GPP_H = A.NPP_H + A.RA_H; GPP_L = A.NPP_L + A.RA_L;
Yrs = Yrs(1:size(GPP_H,1)); r = 0;
for i = min(Yrs):max(Yrs)
    k = find(Yrs == i);
    if numel(k) > 350
        r = r + 1;
        A.GPP_yr(r)  = sum((GPP_H(k,:)   + GPP_L(k,:))   * A.Ccrown');
        A.NPP_yr(r)  = sum((A.NPP_H(k,:) + A.NPP_L(k,:)) * A.Ccrown');
        A.ANPP_yr(r) = sum((A.ANPP_H(k,:)+ A.ANPP_L(k,:))* A.Ccrown');
        A.Yrs_yr(r)  = i;
    end
end
A.Nyr = numel(A.Date)/8766;
PRECIP = (A.Pr_liq + A.Pr_sno)*A.dth;
TRASP  = (A.T_H + A.T_L)*A.dth;
EINT   = (A.ELitter + A.EIn_H + A.EIn_L + A.EIn_urb + A.EIn_rock)*A.dth;
QPER   = sum(A.Qi_out - A.Qi_in, 2);
ESNOW  = (A.EICE + A.ESN + A.ESN_In)*A.dth;
for j = 1:12
    s = A.Datam(:,2) == j;
    A.Tm(j)=sum(TRASP(s)); A.ESNm(j)=sum(ESNOW(s)); A.EGm(j)=sum(A.EG(s));
    A.Einm(j)=sum(EINT(s)); A.Qim(j)=sum(QPER(s));  A.Rhm(j)=sum(A.Rh(s));
    A.Rdm(j)=sum(A.Rd(s));  A.Prm(j)=sum(PRECIP(s)); A.Lkm(j)=sum(A.Lk(s));
    A.TT(j)=A.Tm(j)+A.ESNm(j)+A.EGm(j)+A.Einm(j)+A.Qim(j)+A.Lkm(j)+A.Rhm(j)+A.Rdm(j);
end
end

function draw_monthly(A, which)
n = A.Nyr; m = 1:12; mm = [0, m, 13];
hold on; grid on;
switch which
    case 1   % cumulative monthly fluxes [mm]
        fill(mm, [0 (A.EGm+A.Tm+A.Lkm+A.Qim+A.Einm+A.Rdm+A.Rhm)/n 0], 'y');
        fill(mm, [0 (A.EGm+A.Tm+A.Lkm+A.Qim+A.Einm)/n 0], 'c');
        fill(mm, [0 (A.EGm+A.Tm+A.Lkm+A.Qim)/n 0], 'b');
        fill(mm, [0 (A.EGm+A.Tm)/n 0], 'g');
        fill(mm, [0 A.EGm/n 0], 'r');
        plot(m, A.Prm/n, 'o--k', 'LineWidth', 1.5);
        legend('Runoff','Interc. Evap.','Rec. + Lat. Flow','Transp.','Soil Evap.','Precip.', ...
            'Location','best'); ylabel('[mm]'); xlabel('Month'); xlim([1 12]);
    case 2   % cumulative balance [mm]
        fill(mm, [0 (A.Tm+A.ESNm+A.EGm+A.Einm+A.Qim+A.Lkm+A.Rhm+A.Rdm)/n 0], 'y');
        fill(mm, [0 (A.Rhm+A.Rdm+A.Lkm+A.Qim)/n 0], 'r');
        fill(mm, [0 (A.Rhm+A.Rdm+A.Qim)/n 0], 'b');
        fill(mm, [0 (A.Rhm+A.Rdm)/n 0], 'g');
        fill(mm, [0 A.Rhm/n 0], 'c');
        plot(m, A.Prm/n, 'o--k', 'LineWidth', 1.5);
        legend('Evapotransp.','Rec.','Lat. flow','Sat. Excess','Infilt. Excess','Precip.', ...
            'Location','best'); ylabel('[mm]'); xlabel('Month'); xlim([1 12]);
    case 3   % fractions
        T = A.TT;
        fill(mm, [0 (A.EGm+A.Tm+A.Lkm+A.Qim+A.Einm+A.Rdm+A.Rhm)./T 0], 'y');
        fill(mm, [0 (A.EGm+A.Tm+A.Lkm+A.Qim+A.Einm)./T 0], 'c');
        fill(mm, [0 (A.EGm+A.Tm+A.Lkm+A.Qim)./T 0], 'b');
        fill(mm, [0 (A.EGm+A.Tm)./T 0], 'g');
        fill(mm, [0 A.EGm./T 0], 'r');
        legend('Runoff','Interc. Evap.','Rec. + Lat. Flow','Transp.','Soil Evap.', ...
            'Location','best'); ylabel('Fraction'); xlabel('Month'); axis([1 12 0 1]);
    case 4
        T = A.TT;
        fill(mm, [0 (A.Rhm+A.Rdm+A.Qim+A.ESNm+A.Einm+A.Lkm+A.Tm+A.EGm)./T 0], 'y');
        fill(mm, [0 (A.Rhm+A.Rdm+A.Lkm+A.Qim)./T 0], 'r');
        fill(mm, [0 (A.Rhm+A.Rdm+A.Qim)./T 0], 'b');
        fill(mm, [0 (A.Rhm+A.Rdm)./T 0], 'g');
        fill(mm, [0 A.Rhm./T 0], 'c');
        legend('Evapotransp.','Rec.','Lat. flow','Sat. Excess','Infilt. Excess', ...
            'Location','best'); ylabel('Fraction'); xlabel('Month'); axis([1 12 0 1]);
end
end

function draw_table(mod_param_file, station)
% Every scalar/short-array assignment in the generated MOD_PARAM, laid out in
% three columns. Sl_H and anything LMA-related is dropped: it IS the treatment,
% so it differs between arms and does not belong in a shared parameter table.
txt = fileread(mod_param_file);
lines = regexp(txt, '\r?\n', 'split');
names = {}; vals = {};
for i = 1:numel(lines)
    s = strtrim(regexprep(lines{i}, '%.*$', ''));
    if isempty(s), continue; end
    tok = regexp(s, '^([A-Za-z_]\w*)\s*(\([^)]*\))?\s*=\s*(.+?);\s*$', 'tokens', 'once');
    if isempty(tok), continue; end
    nm = tok{1}; vl = strtrim(tok{3});
    if any(strcmpi(nm, {'Sl_H','Sl_L','mSl_H','mSl_L'})), continue; end   % the treatment
    if numel(vl) > 26, vl = [vl(1:23) '...']; end
    names{end+1} = [nm tok{2}]; %#ok<AGROW>
    vals{end+1}  = vl;          %#ok<AGROW>
end
fh = newfig();
axis off;
n = numel(names); ncol = 3; nrow = ceil(n/ncol);
fs = max(4, min(7, 380/nrow));
for c = 1:ncol
    x0 = 0.02 + (c-1)*0.335;
    for r = 1:nrow
        k = (c-1)*nrow + r;
        if k > n, break; end
        y = 0.95 - (r-1)*(0.92/nrow);
        text(x0, y, names{k}, 'Units','normalized','FontSize',fs, ...
            'FontName','FixedWidth','Interpreter','none');
        text(x0+0.20, y, vals{k}, 'Units','normalized','FontSize',fs, ...
            'FontName','FixedWidth','Interpreter','none','Color',[0 0 0.7]);
    end
end
annotation(fh,'textbox',[0 0.96 1 0.04],'String', ...
    sprintf('%s -- parameters used (LMA/Sl excluded: it is the treatment)  [%d entries]', ...
    station, n), 'EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',9,'Interpreter','none');
end

function fh = newfig()
% 11 x 6.5 inches, the maximum requested, applied to every figure for consistency.
fh = figure('Units','inches','Position',[0 0 11 6.5], ...
            'PaperUnits','inches','PaperSize',[11 6.5], ...
            'PaperPosition',[0 0 11 6.5],'Color','w');
set(gca,'FontSize',9);
end

function save_fig(out_dir, station, name, dpi)
f = fullfile(out_dir, sprintf('%s_%s.png', matlab_name(station), name));
print(gcf, f, '-dpng', sprintf('-r%d', dpi));
close(gcf);
fprintf('  wrote %s\n', f);
end

function s = matlab_name(x), s = strrep(x, '-', '_'); end

function warn(ME, what)
fprintf(2, '  FAILED %s: %s\n', what, ME.message);
for k = 1:min(3, numel(ME.stack))
    fprintf(2, '      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
end
end
