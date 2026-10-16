%% fitParametricModel.m
% Fits a full quadratic (2nd-order polynomial, including all pairwise
% interactions) response-surface model to each output, as an
% alternative to the GP surrogate -- much more data-efficient for
% smooth, moderate nonlinearity since it only estimates a small fixed
% number of coefficients (28 for 6 inputs: 1 intercept + 6 linear + 6
% squared + 15 pairwise interaction terms) via ordinary least squares,
% rather than fitting a flexible kernel that needs far more points to
% constrain reliably -- exactly what the GP's poor leave-one-out R^2
% at 40-55 points showed. Validated the same honest way: leave-one-out
% cross-validation, refitting with each point held out in turn.

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

data = loadAllComsolResults();

inputNames  = {'rx', 'ry', 'rz', 'thx', 'thy', 'thz'};
outputNames = {'Fx', 'Fy', 'Fz', 'Taux', 'Tauy', 'Tauz'};

X = table2array(data(:, inputNames));
n = height(data);

% Normalize inputs to roughly [-1,1] for a well-conditioned design
% matrix -- quadratic terms in raw units span wildly different scales
% (e.g. 1e-5 m vs 0.1 deg), which badly conditions the regression.
inputScale = max(abs(X), [], 1);
Xn = X ./ inputScale;

models = struct();
figure('Color', 'w', 'Name', 'Quadratic model: leave-one-out validation');
fprintf('%8s  %14s  %14s\n', 'output', 'LOO RMSE', 'LOO R^2');
for k = 1:numel(outputNames)
    name = outputNames{k};
    y = data.(name);

    % Full-data fit, kept for querying/reporting later
    mdl = fitlm(Xn, y, 'quadratic');
    models.(name) = mdl;

    % Leave-one-out cross-validation (cheap for OLS -- instant refits)
    yLOO = zeros(n, 1);
    for i = 1:n
        keep = true(n, 1);
        keep(i) = false;
        mdl_i = fitlm(Xn(keep, :), y(keep), 'quadratic');
        yLOO(i) = predict(mdl_i, Xn(i, :));
    end

    rmse = sqrt(mean((yLOO - y).^2));
    ss_res = sum((y - yLOO).^2);
    ss_tot = sum((y - mean(y)).^2);
    r2 = 1 - ss_res / ss_tot;
    fprintf('%8s  %14.4g  %14.3f\n', name, rmse, r2);

    subplot(2, 3, k);
    scatter(y, yLOO, 30, 'filled'); hold on; grid on;
    lims = [min([y; yLOO]), max([y; yLOO])];
    plot(lims, lims, 'r--');
    xlabel('actual (COMSOL)'); ylabel('predicted (LOO)');
    title(sprintf('%s (R^2=%.2f)', name, r2));
end
exportgraphics(gcf, fullfile(scriptDir, 'quadraticModel_leaveOneOut.png'));

%% Direct query: feed in any pose, get predicted outputs back
testPoint = [10e-6, 5e-6, -8e-6, 0.02, -0.01, 0.1];   % [rx ry rz thx thy thz]
testPointN = testPoint ./ inputScale;

fprintf('\nQuery point: rx=%.3g ry=%.3g rz=%.3g thx=%.4g thy=%.4g thz=%.4g\n', testPoint);
fprintf('%8s  %14s  %14s\n', 'output', 'predicted', '95%% CI half-width');
for k = 1:numel(outputNames)
    name = outputNames{k};
    [pred, ci] = predict(models.(name), testPointN);
    fprintf('%8s  %14.4g  %14.4g\n', name, pred, (ci(2) - ci(1)) / 2);
end
