function data = loadAllComsolResults()
%LOADALLCOMSOLRESULTS Finds and combines every real 6D LHS results
%export in this folder (files matching comsol_lhs_*_results.txt, e.g.
%comsol_lhs_initial40_results.txt, comsol_lhs_batch2_results.txt), so
%the training set automatically grows as new adaptive batches are
%exported -- no need to manually re-combine each time a new batch
%comes in. Deliberately excludes comsol_sweep_*_results.txt (the
%single-axis characterization sweeps), which aren't genuinely-coupled
%6D points and would contaminate this dataset.

scriptDir = fileparts(mfilename('fullpath'));
files = dir(fullfile(scriptDir, 'comsol_lhs_*_results.txt'));

if isempty(files)
    error('No comsol_lhs_*_results.txt files found in %s', scriptDir);
end

data = table();
for i = 1:numel(files)
    f = fullfile(files(i).folder, files(i).name);
    batchData = parseComsolResults(f);
    fprintf('  %-40s  %d points\n', files(i).name, height(batchData));
    data = [data; batchData]; %#ok<AGROW>
end

fprintf('Total: %d points from %d file(s)\n', height(data), numel(files));

end
