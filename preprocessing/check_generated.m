function status = check_generated(root)
%CHECK_GENERATED  Parse and lint the MATLAB that build_model_run.py generates.
%
%   status = check_generated(root)      root defaults to $TC_INPUT_DATA/../model_run
%
% checkcode is the code analyser behind the editor's warnings. It reports both
% syntax errors and suspicious code, so one pass covers the two failures that
% motivated this: the unterminated character vector in job 35691 (a syntax error)
% and Kbot_gen used before assignment in 35696 (an undefined name).
%
% mtree, the parser itself, would be the more direct syntax check, but its anyerr
% property is not available in MATLAB 2025a (job 35698 threw on all 8 files), so
% checkcode does both jobs here.
%
% Messages are classified rather than counted: a parse error or an undefined name
% means the file cannot run, while style and unused-value messages are inherited
% from the template and are not this script's business.
%
% Returns 0 if every generated file is runnable, 1 otherwise, so a SLURM wrapper
% can gate the model runs on it.

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

% What makes a generated file unrunnable, as opposed to merely untidy:
%   - a syntax error, which checkcode reports with 'Parse error' in the text
%   - a name that is used but never assigned, which is what a code generator gets
%     wrong when a block lands in the wrong order
% Style and unused-value messages are inherited from the template and are not
% this script's business.
FATAL_IDS   = {'UDIM', 'NODEF', 'STRNU', 'MATLAB:undefinedVarOrFunction'};
FATAL_TEXT  = {'parse error', 'not terminated', 'unexpected', 'unbalanced', ...
               'might be undefined', 'is undefined'};

nerr = 0; nwarn = 0;
fprintf('checking %d generated file(s) under %s\n\n', numel(files), root);
for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);
    rel = strrep(p, [root filesep], '');

    try
        m = checkcode(p, '-id', '-struct');
    catch ME
        fprintf('  CHECK FAILED %s : %s\n', rel, ME.message);
        nerr = nerr + 1;
        continue
    end
    if isempty(m)
        continue
    end

    fatal = false;
    for j = 1:numel(m)
        txt = lower(m(j).message);
        isfatal = any(strcmp(m(j).id, FATAL_IDS)) || ...
                  any(cellfun(@(w) contains(txt, w), FATAL_TEXT));
        if isfatal
            fatal = true;
            fprintf('  FATAL  %s line %d [%s] %s\n', ...
                rel, m(j).line, m(j).id, m(j).message);
        else
            nwarn = nwarn + 1;
        end
    end
    if fatal
        nerr = nerr + 1;
    end
end

fprintf('\n%d file(s) with a parse error or undefined variable\n', nerr);
fprintf('%d non-fatal analyser message(s) (style, unused values)\n', nwarn);
status = double(nerr > 0);
if status == 0
    fprintf('OK -- every generated file parses and defines what it uses\n');
end
end
