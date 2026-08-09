function status = check_generated(root)
%CHECK_GENERATED  Parse and lint the MATLAB that build_model_run.py generates.
%
%   status = check_generated(root)      root defaults to $TC_INPUT_DATA/../model_run
%
% Two different checks, because they catch different things:
%
%   mtree     the real MATLAB parser. Catches syntax errors outright -- the
%             unterminated character vector in job 35691 was one of these, and it
%             cost a full MATLAB startup plus six minutes to discover.
%   checkcode the code analyser (mlint). Catches undefined or possibly-unset
%             variables, which is the class the Kbot_gen ordering bug in job 35696
%             belonged to. Less certain than mtree, since a variable assigned
%             later in the same script may not be flagged, so its output is
%             reported but only specific ids are treated as fatal.
%
% Returns 0 if every generated file parses and has no fatal message, 1 otherwise,
% so a SLURM wrapper can gate the model runs on it.

if nargin < 1 || isempty(root)
    root = fullfile(getenv('TC_INPUT_DATA'), '..', 'model_run');
end
if ~exist(root, 'dir')
    error('check_generated:noRoot', 'not a directory: %s', root);
end

files = [dir(fullfile(root, '*', 'era5_land', '*', 'GO_*.m'));
         dir(fullfile(root, '*', 'era5_land', '*', 'MOD_PARAM_*.m'))];
if isempty(files)
    fprintf('no generated .m files under %s\n', root);
    status = 1;
    return
end

% Undefined or possibly-unset variables. These are what a code generator gets
% wrong -- a name that only exists further down the file, or not at all.
FATAL_IDS = {'MATLAB:undefinedVarOrFunction', 'UDIM', 'NODEF', 'STRNU', 'NASGU'};

nerr = 0; nwarn = 0;
fprintf('checking %d generated file(s) under %s\n\n', numel(files), root);
for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);
    rel = strrep(p, [root filesep], '');

    % 1. Does it parse at all?
    try
        T = mtree(p, '-file');
        if T.anyerr
            [ln, ~, msg] = T.errmsg;
            fprintf('  PARSE ERROR  %s\n', rel);
            for j = 1:numel(ln)
                fprintf('      line %d: %s\n', ln(j), msg{j});
            end
            nerr = nerr + 1;
            continue
        end
    catch ME
        fprintf('  PARSE ERROR  %s : %s\n', rel, ME.message);
        nerr = nerr + 1;
        continue
    end

    % 2. What does the analyser say?
    m = checkcode(p, '-id', '-struct');
    if isempty(m)
        continue
    end
    fatal = false;
    for j = 1:numel(m)
        isfatal = any(strcmp(m(j).id, FATAL_IDS));
        fatal = fatal || isfatal;
        if isfatal
            fprintf('  UNDEFINED    %s line %d [%s] %s\n', rel, m(j).line, m(j).id, m(j).message);
        end
    end
    if fatal
        nerr = nerr + 1;
    else
        nwarn = nwarn + numel(m);
    end
end

fprintf('\n%d file(s) with a parse error or undefined variable\n', nerr);
fprintf('%d non-fatal analyser message(s) (style, unused values)\n', nwarn);
status = double(nerr > 0);
if status == 0
    fprintf('OK -- every generated file parses and defines what it uses\n');
end
end
