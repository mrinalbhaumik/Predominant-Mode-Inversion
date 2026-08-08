function [initPop] = GetInitialPopulation(lbinfo, ubinfo, Layer_param, NP)
% Generates the initial Differential Evolution population within the parameter
% bounds.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

% Build the population column-block by column-block, in the canonical order
% h, Vs, nu, rho, keeping only the variable parameters. Vs uses structured
% (increasing-with-depth) starting profiles; the rest are uniform-random in-bounds.
initPop = [];

if strcmpi(Layer_param.h, 'variable')
    initPop_h = lbinfo.h + rand(NP, length(lbinfo.h)) .* (ubinfo.h - lbinfo.h);
    initPop   = [initPop initPop_h];
end

if strcmpi(Layer_param.Vs, 'variable')
    initPop_vs = Get_Vs_pop(lbinfo, ubinfo, NP);
    initPop    = [initPop initPop_vs];
end

if strcmpi(Layer_param.nu, 'variable')
    initPop_nu = lbinfo.nu + rand(NP, length(lbinfo.nu)) .* (ubinfo.nu - lbinfo.nu);
    initPop    = [initPop initPop_nu];
end

if strcmpi(Layer_param.rho, 'variable')
    initPop_rho = lbinfo.rho + rand(NP, length(lbinfo.rho)) .* (ubinfo.rho - lbinfo.rho);
    initPop     = [initPop initPop_rho];
end

    function [initPop_vs] = Get_Vs_pop(lbinfo, ubinfo, NP)
        top_vs1 = lbinfo.vs(1); bottom_vs1 = lbinfo.vs(end);
        top_vs2 = ubinfo.vs(1); bottom_vs2 = ubinfo.vs(end);
        top_vs    = linspace(top_vs1, top_vs2, NP);
        bottom_vs = linspace(bottom_vs1, bottom_vs2, NP);
        initPop_vs = zeros(NP, length(lbinfo.vs));
        for ii = 1 : NP
            initPop_vs(ii, :) = linspace(top_vs(ii), bottom_vs(ii), length(lbinfo.vs));
        end
    end

end
