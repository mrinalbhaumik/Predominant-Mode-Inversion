function [h, vs, nu, rho] = DE_UnpackVariables(u, N_layer, Layer_param, initial_param)
% Splits the packed DE model vector into per-layer thickness, shear velocity,
% Poisson's ratio and density.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

% Read the variable parameters from u (canonical order: h, Vs, nu, rho); fixed
% ones are taken from the initial model. h has N_layer values; Vs/nu/rho have N_layer+1.
idx = 1;

if strcmpi(Layer_param.h, 'variable')
    h = u(idx : idx + N_layer - 1);   idx = idx + N_layer;
else
    h = initial_param.h_i;
end

if strcmpi(Layer_param.Vs, 'variable')
    vs = u(idx : idx + N_layer);      idx = idx + N_layer + 1;
else
    vs = initial_param.Vs_i;
end

if strcmpi(Layer_param.nu, 'variable')
    nu = u(idx : idx + N_layer);      idx = idx + N_layer + 1;
else
    nu = initial_param.nu_i;
end

if strcmpi(Layer_param.rho, 'variable')
    rho = u(idx : idx + N_layer);
else
    rho = initial_param.rho_i;
end

end
