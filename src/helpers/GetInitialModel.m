
function[initial_param, lbinfo, ubinfo, N_layer, d_max_limit] = GetInitialModel(ii, Target, Layer_param)
% Builds the initial layered earth model and the per-layer lower/upper parameter
% bounds (h, Vs, nu, rho) from the target dispersion curve.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

global WMI_RUN %#ok<GVMIS>
silent = ~isempty(WMI_RUN) && isfield(WMI_RUN,'silent') && WMI_RUN.silent;  % GUI-only: skip limit plots when just computing the table

w_v(:)  = Target.fv_file(:,1);
v_ph(:) = Target.fv_file(:,2);
lambda = v_ph./w_v;

min_val         = min(Layer_param.DCR_range);
max_val         = max(Layer_param.DCR_range);
random_number   = min_val + (max_val - min_val) * rand;
d_max_limit     = max(lambda) / random_number;

%% Upper and Lower limit

lbinfo.all = [];
ubinfo.all = [];

%  layer thickness limits __________________________

if strcmpi(Layer_param.h, 'variable')

    if strcmpi(Layer_param.h_range, 'LR_based')

        LR = Layer_param.h_count(ii);
        d_min       = zeros(1,1);
        d_max       = zeros(1,1);
        d_max_cal   = 0;
        i           = 1;
        while d_max_cal < d_max_limit

            % D_min
            if i==1
                d_min(i) = min(lambda) / 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                d_min(i) = d_max(i-1);
            end

            % D_max
            if i==1
                d_max(i)    = min(lambda)/2;
            elseif i==2
                d_max(i)    = d_min(i) + LR * min(lambda)/2; %%%%%%%%%%%%%%% was 2
            else
                d_max_cal   = d_min(i) + LR * (d_max(i-1) - d_min(i-1));
                % d_max(i) = min(d_max_cal, d_max_limit);
                if d_max_cal > d_max_limit
                    d_max(i-1)    = d_max_limit;
                    d_min(i) = [];
                else
                    d_max(i)    = d_max_cal;
                end
            end

            i = i + 1;
        end

        lb_h1 = [d_min(1) diff(d_min)];
        ub_h1 = [d_max(1) diff(d_max)];

        lb_h = min(lb_h1, ub_h1);
        ub_h = max(lb_h1, ub_h1);

        N_layer = length(lb_h);

        lbinfo.h   = lb_h;
        ubinfo.h   = ub_h;
        lbinfo.all = [lbinfo.all lb_h];
        ubinfo.all = [ubinfo.all ub_h];

    elseif strcmpi(Layer_param.h_range, 'LN_based')

        N_layer    = Layer_param.h_count(ii);
        lb_h       = repelem(min(lambda)/3, N_layer);%%%%%%%%%%
        ub_h       = repelem(max(d_max_limit ./ (N_layer-1), 10), N_layer); % Change it later

        d_min      = cumsum(lb_h);
        d_max      = cumsum(ub_h);

        lbinfo.h   = lb_h;
        ubinfo.h   = ub_h;
        lbinfo.all = [lbinfo.all lb_h];
        ubinfo.all = [ubinfo.all ub_h];

    else
        error('Correctly choose the h_range option, either LR_based or LN_based')
    end

elseif strcmpi(Layer_param.h, 'fixed')
    N_layer = Layer_param.h_count(ii);
else
    error('Correctly choose the h parameter, either variable or fixed')
end


%  layer Vs limits __________________________


if numel(Layer_param.Vs_range)==2 % if the Vs range is provided
    lb_vs = zeros(1, N_layer+1);
    ub_vs = zeros(1, N_layer+1);
    lb_vs(:) = min(Layer_param.Vs_range);
    ub_vs(:) = max(Layer_param.Vs_range);
else
    % wavelength dependent
    [lb_vs, ub_vs] = Get_WavelengthBased_Vs_Limits(N_layer,lambda, v_ph, d_min, d_max);
end


lbinfo.vs  = lb_vs;
ubinfo.vs  = ub_vs;
lbinfo.all = [lbinfo.all lb_vs];
ubinfo.all = [ubinfo.all ub_vs];

% Poisson's ratio limit _____________________________________
if strcmpi(Layer_param.nu, 'variable')
    if numel(Layer_param.nu_range) < 2
        error('Alart !!! Provide a range for nu');
    end
    lb_nu = min(Layer_param.nu_range) .* ones(1, N_layer+1);
    ub_nu = max(Layer_param.nu_range) .* ones(1, N_layer+1);

    lbinfo.nu  = lb_nu;
    ubinfo.nu  = ub_nu;
    
    lbinfo.all = [lbinfo.all lb_nu];
    ubinfo.all = [ubinfo.all ub_nu];
end

% Density limits _____________________________________
if strcmpi(Layer_param.rho, 'variable')
    if numel(Layer_param.rho_range) < 2
        error('Alart !!! Provide a range for density');
    end
    lb_rho = min(Layer_param.rho_range) .* ones(1, N_layer+1);
    ub_rho = max(Layer_param.rho_range) .* ones(1, N_layer+1);

    lbinfo.rho  = lb_rho;
    ubinfo.rho  = ub_rho;
    lbinfo.all = [lbinfo.all lb_rho];
    ubinfo.all = [ubinfo.all ub_rho];
end

%% Initial model
initial_param.all = [];
% Initial layer thickness __________________________

if strcmpi(Layer_param.h, 'fixed')

    if strcmpi(Layer_param.fixed, 'increasing')
        N_layer = Layer_param.h_count(ii);
        % Based on depth conversion ratio (DCR)
        LR           = unifrnd(1.25, 1.5, 1, 1);
        DCR          = (LR ^ (N_layer + 1) - LR) / ((LR-1) * LR);
        h1           = (1 / DCR) * d_max_limit;
        h_i          = [1, LR * ones(1, N_layer-1)];
        h_i          = h1.* cumprod(h_i);
    elseif strcmpi(Layer_param.fixed, 'equal')
        N_layer      = Layer_param.h_count(ii);
        h_i          = repelem(d_max_limit/N_layer, N_layer);
    else
        error('Correctly choose the type of layer for fixed h, either increasing or equal')
    end

elseif strcmpi(Layer_param.h, 'variable')

    if strcmpi(Layer_param.h_range, 'LR_based')

        % Based on depth conversion ratio (DCR)
        LR           = unifrnd(1.25, 2, 1, 1);
        DCR          = (LR ^ (N_layer + 1) - LR) / ((LR-1) * LR);
        h1           = (1 / DCR) * d_max_limit;
        h_i          = [1, LR * ones(1, N_layer-1)];
        h_i          = h1.* cumprod(h_i);
        initial_param.all = [initial_param.all h_i];
    elseif strcmpi(Layer_param.h_range, 'LN_based')

        h_i               = repelem(d_max_limit/N_layer, N_layer);
        initial_param.all = [initial_param.all h_i];
    else
        error('Correctly choose the h_range option, either LR_based or LN_based')
    end

else
    error('Correctly choose the h type, either fixed or variable')
end

d_i = cumsum(h_i);

% Initial Vs_____________________________________

    function v = make_unique(v)
        [s, j] = sort(v(:));
        s      = s + 0.1 * (0:length(s)-1)';
        v(j)   = s;
    end

lambda       = make_unique(lambda);
vs_i         = interp1(lambda, v_ph, d_i*3.5, 'linear', 'extrap'); % adjust the value to get a closer initial model
vs_i_HS      = 1.5 * max(vs_i) ;
vs_i         = [vs_i vs_i_HS];
initial_param.all = [initial_param.all vs_i];

% Initial Poisson's ratio _____________________________________
if strcmpi(Layer_param.nu, 'fixed')
    nu_i = repelem(Layer_param.nu_range(1), N_layer+1);
elseif strcmpi(Layer_param.nu, 'variable')
    nu_i  = unifrnd(Layer_param.nu_range(1), Layer_param.nu_range(2), 1, N_layer+1);
    initial_param.all = [initial_param.all nu_i];
else
    error('Correctly select the nu type, either fixed or variable')
end

% Initial density _____________________________________
if strcmpi(Layer_param.rho, 'fixed')
    rho_i = repelem(Layer_param.rho_range(1), N_layer+1);
elseif strcmpi(Layer_param.rho, 'variable')
    rho_i  = unifrnd(Layer_param.rho_range(1), Layer_param.rho_range(2), 1, N_layer+1);
    initial_param.all = [initial_param.all rho_i];
else
    error('Correctly select the density type, either fixed or variable')
end

initial_param.h_i   = h_i;
initial_param.Vs_i  = vs_i;
initial_param.nu_i  = nu_i;
initial_param.rho_i = rho_i;

%% Plot limits
% Standalone (non-GUI) path plots the computed limits here. In GUI mode the
% dashboard is present and the caller re-plots the (possibly edited) limits
% after applyLayerLimits, so we skip plotting here to avoid a stale envelope.
guiDash = ~isempty(WMI_RUN) && isfield(WMI_RUN,'dash') && isstruct(WMI_RUN.dash) ...
          && isfield(WMI_RUN.dash,'vs') && isgraphics(WMI_RUN.dash.vs);   % only when the GUI dashboard is live
if ~silent && ~guiDash
    PlotLayerLimits(lbinfo, ubinfo, Layer_param, d_max_limit);
end

%% This function calculate Vs limits based on wavelength

    function[lb_vs, ub_vs] = Get_WavelengthBased_Vs_Limits(N_layer,lambda, v_ph, d_min, d_max)

        vs_min = zeros(1,N_layer);
        vs_max = zeros(1,N_layer);

        % This section require when there is clear mode jump, i.e., same wavelength
        % may have two velocity.

        [lambda, ~, lambda_idx] = (unique(lambda));
        v_ph                     = accumarray(lambda_idx, v_ph, [], @mean);
        v_ph                     = flip(v_ph); 
        lambda                  = flip(lambda);
        %

        for jj = 1 : N_layer

            Lambda_min = d_min(jj) * 1.5;
            vs_min(jj) = interp1(lambda, v_ph, Lambda_min, 'linear', 'extrap');

            Lambda_max = d_max(jj) * 3.5;
            vs_max(jj) = interp1(lambda, v_ph, Lambda_max, 'linear', 'extrap');

        end
        vs_min_HS = max(vs_min);
        vs_max_HS = max(v_ph) * 2;
        lb_vs1 = [vs_min vs_min_HS];
        ub_vs1 = [vs_max vs_max_HS];
        lb_vs = min(lb_vs1,ub_vs1);
        ub_vs = max(lb_vs1,ub_vs1);

        lb_vs = 0.6 * lb_vs; % additional lower bound room.

    end

end