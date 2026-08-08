function PlotLayerLimits(lbinfo, ubinfo, Layer_param, d_max_limit)
% Draws the per-layer shear-velocity bound envelope (lower/upper limits) into the
% model dashboard axes (GUI) or figure 999.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

% --- Vs limits ---
if strcmpi(Layer_param.h, 'variable') && isfield(lbinfo,'vs')
    ax = PM_Inv_DashAxes('vs', 2,7, [3:4 10:11]); hold(ax,'on');
    lb_vs_plot = repelem(lbinfo.vs, 2);
    lb_d_plot  = [0 repelem(cumsum(ubinfo.h),2) d_max_limit];
    plot(ax, lb_vs_plot, lb_d_plot, '--', 'LineWidth', 1)

    ub_vs_plot = repelem(ubinfo.vs, 2);
    ub_d_plot  = [0 repelem(cumsum(lbinfo.h),2) d_max_limit];
    plot(ax, ub_vs_plot, ub_d_plot, '--', 'LineWidth', 1)
    axis(ax,'ij'); ylabel(ax,'Depth (m)'); xlabel(ax,'Shear wave velocity (m/s)')
end

end
