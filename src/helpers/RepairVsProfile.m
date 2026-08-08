function u = RepairVsProfile(u, lb, ub, Layer_param, initial_param, N_layer)
% Enforces structural constraints on a candidate model: forces the half-space to be
% the fastest layer and limits shear-velocity reversals below a set depth.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

[h, vs, nu, rho] = DE_UnpackVariables (u, N_layer, Layer_param, initial_param);

%% This section makes the half-space velocity maximum
if strcmpi(Layer_param.Vs_halfspace_max, 'yes')
    if vs(end) < max(vs(1:end-1))
        vs(end) = max(vs);
    end
end
%% This section ensure no velocity reversal after a given depth
rev_lim = Layer_param.Vs_rev_lim;
if rev_lim > 0
    Depth = cumsum(h);
    [~, ind] = find(Depth > rev_lim);
    if ~isempty(ind)
        for LL = 1 : length(ind)
            if ind(LL) > 1
                if vs(ind(LL)) < vs(ind(LL)-1)
                    vs(ind(LL)) = vs(ind(LL)-1);
                end
            end
        end
    end
end

%% Reform u from the updated parameters

[u] = DE_PackVariables (h, vs, nu, rho, Layer_param);

end
