function ax = PM_Inv_DashAxes(key, m, n, p)
% Resolves the target axes for a dashboard panel: the GUI uiaxes when the interactive
% app is active, otherwise the classic figure-999 subplot.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

    global WMI_RUN %#ok<GVMIS>
    if ~isempty(WMI_RUN) && isfield(WMI_RUN,'dash') && isstruct(WMI_RUN.dash) ...
            && isfield(WMI_RUN.dash, key) && isgraphics(WMI_RUN.dash.(key))
        ax = WMI_RUN.dash.(key);
    else
        figure(999);
        ax = subplot(m, n, p);
    end
end
