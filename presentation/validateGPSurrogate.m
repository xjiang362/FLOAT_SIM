%% validateGPSurrogate.m
% Validates the fitted GP surrogate against the real COMSOL data three
% ways: (1) leave-one-out cross-validation (predicted vs actual on
% held-out points -- the honest generalization check, since a GP
% trivially reproduces its own training points), (2) 1D marginal
% slices through the 6D space for a visual sanity check against
% physical intuition, and (3) a direct query section to feed in any
% pose and get predicted outputs back.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = loadAllComsolResults();

inputNames  = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
outputNames = {'Fx', 'Fy', 'Fz', 'Taux', 'Tauy', 'Tauz'};

X = table2array(data(:, inputNames));

%% Fit one GP per output (same setup as adaptiveGPNextBatch.m)
gprModels = struct();
yMean = struct(); yStd = struct();
for i = 1:numel(outputNames)
    y = data.(outputNames{i});
    yMean.(outputNames{i}) = mean(y);
    yStd.(outputNames{i})  = std(y);
    yNorm = (y - yMean.(outputNames{i})) / yStd.(outputNames{i});
    gprModels.(outputNames{i}) = fitrgp(X, yNorm, ...
        'KernelFunction', 'ardsquaredexponential', ...
        'Standardize', true, ...
        'BasisFunction', 'constant');
end

%% 1. Leave-one-out cross-validation: predicted vs actual
figure('Color', 'w', 'Name', 'Leave-one-out validation');
fprintf('%8s  %10s  %10s\n', 'output', 'RMSE', 'R^2');
for i = 1:numel(outputNames)
    name = outputNames{i};
    y = data.(name);
    yNorm = (y - yMean.(name)) / yStd.(name);

    cvMdl = crossval(gprModels.(name), 'Leaveout', 'on');
    yLOOnorm = kfoldPredict(cvMdl);
    yLOO = yLOOnorm * yStd.(name) + yMean.(name);   % back to real units

    rmse = sqrt(mean((yLOO - y).^2));
    ss_res = sum((y - yLOO).^2);
    ss_tot = sum((y - mean(y)).^2);
    r2 = 1 - ss_res / ss_tot;
    fprintf('%8s  %10.3g  %10.3f\n', name, rmse, r2);

    subplot(2, 3, i);
    scatter(y, yLOO, 30, 'filled'); hold on; grid on;
    lims = [min([y; yLOO]), max([y; yLOO])];
    plot(lims, lims, 'r--');   % perfect-prediction reference line
    xlabel('actual (COMSOL)'); ylabel('predicted (LOO)');
    title(sprintf('%s (R^2=%.2f)', name, r2));
end
exportgraphics(gcf, fullfile(scriptDir, 'gpValidation_leaveOneOut.png'));

%% 2. 1D marginal slices through the 6D space
bounds = max(abs(X), [], 1);   % per-input range from the training data itself
nGrid = 100;

figure('Color', 'w', 'Name', '1D marginal slices (all other inputs held at 0)');
for i = 1:6
    grid_i = linspace(-bounds(i), bounds(i), nGrid)';
    Xslice = zeros(nGrid, 6);
    Xslice(:, i) = grid_i;

    subplot(2, 3, i); hold on; grid on;
    for k = 1:numel(outputNames)
        name = outputNames{k};
        [muNorm, sdNorm] = predict(gprModels.(name), Xslice);
        mu = muNorm * yStd.(name) + yMean.(name);
        plot(grid_i, mu, 'DisplayName', name, 'LineWidth', 1.2);
    end
    xlabel(inputNames{i}); ylabel('predicted output');
    title(sprintf('slice along %s', inputNames{i}));
    if i == 1, legend('Location', 'best'); end
end
exportgraphics(gcf, fullfile(scriptDir, 'gpValidation_marginalSlices.png'));

%% 3. Direct query: feed in any pose, get predicted outputs back
% Edit testPoint below to whatever combination you want to check.
testPoint = [10e-6, 5e-6, -8e-6, 0.02, -0.01, 0.1];   % [rx ry rz thx thy thz]

fprintf('\nQuery point: rx=%.3g ry=%.3g rz=%.3g thx=%.4g thy=%.4g thz=%.4g\n', testPoint);
fprintf('%8s  %14s  %14s\n', 'output', 'predicted', 'std');
for i = 1:numel(outputNames)
    name = outputNames{i};
    [muNorm, sdNorm] = predict(gprModels.(name), testPoint);
    mu = muNorm * yStd.(name) + yMean.(name);
    sd = sdNorm * yStd.(name);   % std scales linearly, no mean shift needed
    fprintf('%8s  %14.4g  %14.4g\n', name, mu, sd);
end
