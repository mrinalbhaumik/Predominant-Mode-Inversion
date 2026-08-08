function u = DE_PackVariables(h, vs, nu, rho, Layer_param)
% Assembles per-layer thickness, shear velocity, Poisson's ratio and density into
% the packed DE model vector (elastic version — no damping).
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

% Keep only the variable parameters, in the canonical order h, Vs, nu, rho.
u = [];
if strcmpi(Layer_param.h,  'variable'); u = [u h];   end
if strcmpi(Layer_param.Vs, 'variable'); u = [u vs];  end
if strcmpi(Layer_param.nu, 'variable'); u = [u nu];  end
if strcmpi(Layer_param.rho,'variable'); u = [u rho]; end

end
