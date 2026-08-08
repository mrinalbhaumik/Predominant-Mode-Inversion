function [] = PM_Inv_PlotProfile (gen, fbest, h_vs_nu, N_layer ,Layer_param, initial_param, HTLM, Target)
% Plots the current best model (Vs-depth profile and modeled dispersion) each
% generation and records the iteration history in the global 'param'.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

[h, vs, nu, rho] = DE_UnpackVariables (h_vs_nu, N_layer, Layer_param, initial_param);

w = Target.fv_file(:,1);

[v, Domi_v] = PM_Inv_Forward(vs, nu, rho, h, w, HTLM, Target.Measured_component);

global param %#ok<GVMIS>

param{1,gen} = gen;
param{2,gen} = fbest;
param{3,gen} = h;
param{4,gen} = vs;
param{5,gen} = nu;
param{6,gen} = rho;
param{7,gen} = v(:,1:5);
param{8,gen} = Domi_v;


%% Plot Vs profiles

ax = PM_Inv_DashAxes('vs', 2,7, [3:4 10:11]); hold(ax,'on')
obj = findobj(ax,'Type','line');
if length(obj) >= 3
    delete(obj(1))
end
if isrow(h)
    h = h';
end
h(isnan(h)) = [];
vs_plot     = repelem(vs, 2);
% d_act       = repelem(h, 2);
% MaxDepth    = h(end);
d_act       = repelem(cumsum(h), 2);
MaxDepth    = sum(h);

if isrow(d_act)
    d_act   = d_act';
end
d_act       = [0; d_act; 1.2*MaxDepth];
plot(ax, vs_plot, d_act, '-', 'LineWidth', 1)
axis(ax,'ij')

ax = PM_Inv_DashAxes('vszoom', 2,7, 7); hold(ax,'on')
lines = findobj(ax, 'Type', 'line');
errs  = findobj(ax, 'Type', 'errorbar');
obj   = [lines; errs];
if length(obj) >= 1
    delete(obj(1))
end
plot(ax, vs_plot, d_act, '-', 'LineWidth', 1); axis(ax,'ij'); ylim(ax,[0 0.5*MaxDepth])

% plot dispersion curve

ax = PM_Inv_DashAxes('disp', 2,7, 1:2); hold(ax,'on')
axis(ax,'tight');
lines = findobj(ax, 'Type', 'line');
errs  = findobj(ax, 'Type', 'errorbar');
obj   = [lines; errs];
if length(obj)>1
    delete(obj(1:end-1))
end
plot(ax, w, Domi_v, 'or', 'MarkerSize', 3);
plot(ax, w, v(:,1), '-m', w, v(:,2), '-g', w, v(:,3), '-b'); % Plot theoretical

pause(0.1)
end
