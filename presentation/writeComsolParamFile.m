function writeComsolParamFile(outFile, paramNames, paramUnits, paramData)
%WRITECOMSOLPARAMFILE Write a COMSOL Parametric Sweep "Load from File"
%table: one row per parameter, name "v1 v2 ..." [unit]. Confirmed
%against a file COMSOL itself exported (a single-value row read back
%as "thx 0 [rad]") -- name/value/unit are space-separated, not
%comma-separated.
%
%   outFile     destination path
%   paramNames  cell array of parameter name strings
%   paramUnits  cell array of unit strings (COMSOL unit syntax, e.g. 'm', 'deg')
%   paramData   cell array of numeric vectors, one per parameter, all
%               the same length (one entry per combination)

fid = fopen(outFile, 'wt');   % text mode: writes CRLF row endings
for i = 1:numel(paramNames)
    valsStr = strtrim(sprintf('%.10g ', paramData{i}));
    fprintf(fid, '%s "%s" [%s]\n', paramNames{i}, valsStr, paramUnits{i});
end
fclose(fid);

end
