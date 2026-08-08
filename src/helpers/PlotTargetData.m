function [] = PlotTargetData(Target)
% Plots the target phase-velocity dispersion curve into the interactive dashboard
% axes (GUI) or the standalone figure 999.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

% Draw into the GUI dashboard axes when present; otherwise into standalone figure 999.
global WMI_RUN %#ok<GVMIS>
guiOn = ~isempty(WMI_RUN) && isfield(WMI_RUN,'dash') && isstruct(WMI_RUN.dash) ...
        && isfield(WMI_RUN.dash,'disp') && isgraphics(WMI_RUN.dash.disp);

if guiOn
    % Clear the whole dashboard at the start of a trial (the standalone path
    % gets the same fresh start from 'close all' between trials).
    axList = struct2cell(WMI_RUN.dash);
    for kk = 1:numel(axList)
        if isgraphics(axList{kk}); cla(axList{kk}); end
    end
else
    figure(999);
    figPos = GetFigurePosition(0.55, 0.4, 5, 3); set(gcf, 'Position', figPos);
end

% --- Frequency - phase velocity ---
[~, numCols] = size(Target.fv_file);
ax = PM_Inv_DashAxes('disp', 2,7, 1:2);
if numCols >= 3
    errorbar(ax, Target.fv_file(:,1), Target.fv_file(:,2), Target.fv_file(:,3), ...
        'k.-', 'LineWidth', 1, 'MarkerSize', 5, 'CapSize', 4);
else
    plot(ax, Target.fv_file(:,1), Target.fv_file(:,2), 'k.-', 'LineWidth', 1, 'MarkerSize', 10);
end
ylabel(ax,'Phase velocity (m/s)'); xlabel(ax,'Frequency (Hz)'); hold(ax,'off')
set(ax, 'XScale', 'log'); set(ax, 'YScale', 'log');

end
