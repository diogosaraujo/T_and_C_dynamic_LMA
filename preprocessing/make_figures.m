function status = make_figures(rundir, force)
%MAKE_FIGURES  Draw GRAPH_MOD figures from a finished run, without re-running it.
%
%   status = make_figures(rundir)          skip if figures/ already has PNGs
%   status = make_figures(rundir, true)    redraw regardless
%
% GO_<ST>.m saves the WHOLE workspace into RES_<ST>.mat, so loading that file
% restores everything GRAPH_MOD needs. A script called from a function runs in
% that function's workspace, so GRAPH_MOD sees the loaded variables directly --
% no need to push them to the base workspace.
%
% This exists because GRAPH_MOD's inputdlg calls threw under 'matlab -batch' and
% aborted before the first figure, so 182 completed runs saved results but no
% figures. Those runs are perfectly good; only the plotting has to be redone, and
% at ~2 min per arm that is far cheaper than repeating ~0.5 h simulations.
%
% Returns 0 figures written, 2 nothing to do (no RES, or already drawn), 1 failed.

if nargin < 2 || isempty(force)
    force = false;
end
if ~exist(rundir, 'dir')
    error('make_figures:noDir', 'not a directory: %s', rundir);
end

d = dir(fullfile(rundir, 'RES_*.mat'));
if isempty(d)
    fprintf('  no RES_*.mat in %s -- run not finished, nothing to plot\n', rundir);
    status = 2;
    return
end
if numel(d) > 1
    % Ambiguous input is worth stopping for: picking one silently could pair a
    % station's figures with another station's results.
    error('make_figures:manyRES', '%d RES_*.mat files in %s', numel(d), rundir);
end

fig_dir = fullfile(rundir, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
elseif ~force
    % Skip only if the figures are NEWER than the results. Keying the skip on
    % existence alone silently left stale PNGs in place after a re-run, so the
    % figures showed one model configuration while RES held another -- and the
    % only way to notice was to remember which came first.
    png = dir(fullfile(fig_dir, '*.png'));
    if ~isempty(png) && max([png.datenum]) >= d(1).datenum
        fprintf('  figures in %s are newer than %s -- nothing to do\n', ...
            fig_dir, d(1).name);
        status = 2;
        return
    elseif ~isempty(png)
        fprintf('  figures are OLDER than %s -- redrawing\n', d(1).name);
    end
end

% GRAPH_MOD lives at the model_run root and the T&C source in root/Code, the same
% two directories GO_<ST>.m adds. Derived from rundir rather than assumed
% relative, so this works from any working directory:
%   <root>/<STATION>/era5_land/<arm>  ->  <root>
root = fileparts(fileparts(fileparts(rundir)));
addpath(fullfile(root, 'Code'));
addpath(root);

resfile = fullfile(d(1).folder, d(1).name);
fprintf('  loading %s (%.0f MB)\n', d(1).name, d(1).bytes / 1e6);
load(resfile);                                          %#ok<LOAD> whole workspace

status = 1;
old = cd(rundir);
restore = onCleanup(@() cd(old));
try
    set(0, 'DefaultFigureVisible', 'off');
    GRAPH_MOD;
    h = findobj('Type', 'figure');
    for kf = 1:numel(h)
        nm = get(h(kf), 'Name');
        if isempty(nm)
            nm = sprintf('fig%02d', get(h(kf), 'Number'));
        end
        nm = regexprep(nm, '[^0-9A-Za-z_-]', '_');
        saveas(h(kf), fullfile('figures', sprintf('%s_%s.png', id_location, nm)));
    end
    close all
    fprintf('  wrote %d figure(s) to %s\n', numel(h), fig_dir);
    status = 0;
catch ME
    close all force
    fprintf(2, '  GRAPH_MOD failed for %s: %s\n', rundir, ME.message);
    for k = 1:numel(ME.stack)
        fprintf(2, '      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
end
end
