%% plotRzZoom.m
% Heave only, first 20 ms, from the saved compareDesignModels run:
% exact sampled map vs double-integrator design model.
%
% The full 150 ms view is dominated by the slow uncontrolled lateral drift,
% which hides the part that actually distinguishes the two design models --
% the capture transient. This is that window.
%
% Regenerates the data by running compareDesignModels if the .mat is absent.

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'config'), fullfile(root,'dynamics'), ...
        fullfile(root,'control'), fullfile(root,'estimator'), ...
        fullfile(root,'presentation'), root);

matFile = fullfile(root,'presentation','CompareDesignModels.mat');
if ~isfile(matFile)
    fprintf('%s not found -- running compareDesignModels first.\n', matFile);
    compareDesignModels; clear;
    root = fileparts(fileparts(mfilename('fullpath')));
    matFile = fullfile(root,'presentation','CompareDesignModels.mat');
end
S = load(matFile);
res = S.res; nmv = S.nmv; col = S.col; x0 = S.x0;

tZoom = 20;     % [ms]

f = figure('Color','w','Position',[60 60 1250 520]);

% ---- (a) full zoom window ----
subplot(1,2,1); hold on; grid on; box on;
for i = 1:numel(res)
    L = res(i).L; if isempty(L), continue; end
    sel = L.t*1e3 <= tZoom;
    plot(L.t(sel)*1e3, L.x(sel,3)*1e6, '-', 'LineWidth',2.0, ...
         'Color',col(i,:), 'DisplayName',nmv{i});
    if ~isempty(L.contact) && L.contactTime*1e3 <= tZoom
        plot(L.contactTime*1e3, L.x(end,3)*1e6, 'x', 'Color',col(i,:), ...
             'MarkerSize',14,'LineWidth',2.5,'HandleVisibility','off');
        text(L.contactTime*1e3+0.4, L.x(end,3)*1e6, ...
             sprintf('contact %.1f ms', L.contactTime*1e3), ...
             'Color',col(i,:),'FontSize',8);
    end
end
yline(0,'k:','HandleVisibility','off');
yline( 15,'k--','LineWidth',1.1,'DisplayName','wall (\pm15 \mum)');
yline(-15,'k--','LineWidth',1.1,'HandleVisibility','off');
xlabel('t [ms]'); ylabel('r_z [\mum]'); xlim([0 tZoom]); ylim([-15.5 15.5]);
title('heave, first 20 ms'); legend('Location','southeast','FontSize',8);

% ---- (b) same window, tight on the settling band ----
subplot(1,2,2); hold on; grid on; box on;
for i = 1:numel(res)
    L = res(i).L; if isempty(L), continue; end
    sel = L.t*1e3 <= tZoom;
    plot(L.t(sel)*1e3, L.x(sel,3)*1e6, '-', 'LineWidth',2.0, 'Color',col(i,:));
end
yline(0,'k:');
xlabel('t [ms]'); ylabel('r_z [\mum]'); xlim([0 tZoom]); ylim([-2 2]);
title('same window, zoomed on the setpoint');

sgtitle(sprintf(['Heave capture: exact sampled map vs double integrator  \\bullet  ' ...
    'start (%.0f, %.0f, %.0f) \\mum'], x0(1)*1e6, x0(2)*1e6, x0(3)*1e6),'FontSize',11);
for ax = findobj(f,'Type','axes')', ax.Toolbar.Visible='off'; end
out = fullfile(root,'presentation','CompareDesignModels_rzZoom.png');
exportgraphics(f, out, 'Resolution',140);
fprintf('Saved %s\n', out);

%% ---- numbers over the zoom window ----
fprintf('\n%-20s %12s %12s %12s %12s\n','model','min rz[um]','max rz[um]','|rz|<0.1um','fate');
for i = 1:numel(res)
    L = res(i).L; if isempty(L), continue; end
    sel = L.t*1e3 <= tZoom;
    z = L.x(sel,3)*1e6; t = L.t(sel)*1e3;
    k = find(abs(z) > 0.1, 1, 'last');
    if isempty(k) || k >= numel(z), ts = 'not reached';
    else, ts = sprintf('%.1f ms', t(k+1)); end
    if isempty(L.contact), fate = 'survives';
    else, fate = sprintf('contact %.1f ms', L.contactTime*1e3); end
    fprintf('%-20s %12.3f %12.3f %12s %12s\n', nmv{i}, min(z), max(z), ts, fate);
end
