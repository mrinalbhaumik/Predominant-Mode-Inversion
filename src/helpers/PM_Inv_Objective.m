function f_phase = PM_Inv_Objective(h_vs_nu, N_layer, HTLM, Layer_param, initial_param, Target)
% Objective (misfit) function: runs the forward model and returns the RMS error
% between the modeled and target dominant-mode phase velocity.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

%% Unpack the model vector ------------------------------------------------------------------------
[h, vs, nu, rho] = DE_UnpackVariables(h_vs_nu, N_layer, Layer_param, initial_param);

w  = Target.fv_file(:,1);

%% Phase-velocity misfit --------------------------------------------------------------------------
if any(diff(vs) < -Layer_param.Vs_reversal * vs(1:end-1))   % velocity reversal beyond allowed drop
    f_phase = 1e6;                 % penalize: reject this model

else
    [~, Domi_v] = PM_Inv_Forward(vs, nu, rho, h, w, HTLM, Target.Measured_component);
    D_R = Domi_v;

    if size(Target.fv_file,2) > 2
        v_f     = Target.fv_file(:,2);
        std_vel = Target.fv_file(:,3);
        elem_1  = (D_R(:)-v_f(:)).^2;
        f_phase = sqrt( sum(elem_1(:) ./ std_vel(:).^2) / length(v_f) );

    else
        v_f     = Target.fv_file(:,2);
        elem_1  = (D_R(:)-v_f(:)).^2;
        f_phase = sqrt( sum(elem_1) / length(v_f) );
    end
end

end
