
% =====================================================================================================================================
% ============================= Differential Evolution (DE) for Predominant Mode Inversion (PMinv) ==============================
% =====================================================================================================================================
% It has both VCPM (Vertical component predominant mode) and RCPM (Radial
% component predominant mode)
% The GUI Version of this code is also available in the same directory.

%  Author : [Mrinal Bhaumik]
%  Date   : [Jun_2026]
%  Version: 1.0
% ------------------------------------------------------------
%  Description:
%  Main driver script for waveform/dispersion inversion of surface-wave
%  data. It imports a target dispersion curve, preprocesses it, builds an
%  initial layered earth model with parameter bounds, and inverts for the
%  layer properties (thickness, Vs, nu, rho) using a Differential
%  Evolution optimizer. Runs over the requested layering ratios and trials,
%  saving inputs and per-run results to an output folder named after the
%  target file.

%   This package will be **monitored, maintained, and upgraded** by the authors over time 
%   to ensure continued functionality, compatibility, and performance improvements.
%
%   The authors provide no warranty; use at your own risk.
%
% CONTACT / SUPPORT:
%   - Mrinal Bhaumik:  mrinal.bhaumik2012@gmail.com
%                   :  mrinal.bhaumik@usu.edu
%
%   For issues, feedback, collaboration opportunities, or any explanation or suggestions, 
%   please feel free to reach out via email. We welcome discussion and user input 
%   to further improve this package.
%
% IMPORTANT NOTE ON PERFORMANCE:
%   The performance and accuracy of the inversion heavily depend on the model parameterization. 
%   It is strongly recommended to experiment with different parameterizations and perform 
%   multiple trial runs to achieve stable and geophysical meaningful results.

% This software accompanies the publication: 
% [1]. Bhaumik, M., & Cox, B. R. (2026). Radial-Component Predominant-Mode Inversion of Rayleigh Waves:
% Application to DAS-based Site Characterization. arXiv preprint arXiv:2605.16717. https://doi.org/10.48550/arXiv.2605.16717

% [2]. Bhaumik, M., & Cox, B. R. (2026). Predominant-mode inversion of surface waves: Inherently addressing 
% inconspicuous low frequency mode jumps. Engineering Geology, 108834. https://doi.org/10.1016/j.enggeo.2026.108834

% [3]. Bhaumik, M., & Naskar, T. (2024). Computation of Surface Wave's Dominating Mode for Stratified Media.
% Indian Geotechnical Journal, 55(3), 1699-1714. https://doi.org/10.1007/s40098-024-01045-x

% ------------------------------------------------------------

%% Initialize session ===================================================================================================================
% clc; 
clearvars; close all       % clearvars (not 'clear all') keeps breakpoints/compiled fcns/parpool intact
disp('Running Inversion')
dbstop if error

% Add the helper functions folder (relative to this script) to the path
scriptDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(scriptDir, 'helpers')));

% Ensure a parallel pool exists (the DE solver uses parfor); print status
if license('test','Distrib_Computing_Toolbox')
    if isempty(gcp('nocreate'))
        parpool;
    end
    fprintf('Parallel pool active: %d workers.\n', gcp().NumWorkers);
else
    fprintf('Parallel Computing Toolbox not available - running serially.\n');
end

%% Target data — import & preprocess ====================================================================================================

% Import target curve: (.txt, .xlsx, .csv)
% Col-1:Frequency; Col-2:Phase velocity/slowness; Col-3: std_dev (optional)
% Path is resolved relative to this script so it runs from any working directory

Target.fv_file = fullfile(scriptDir, '50m_R0_dispersion_stats.txt');
Target.Measured_component = 'VCPM';   % measured component: 'VCPM' or 'RCPM' (used in PM_Inv_Forward > filtering)

% Preprocess data:
% Convert slowness to velocity, handle std, resample if required.
Include_min_COV = [];       % Include minimum COV(e.g, 0.05,0.1), if not keep it empty.
Resample        = "no";     % Set "yes" for log-scale resampling. Set to "no" if resampling not needed
Number          = [];       % if 'Resample' is 'yes', Set number of log samples
Stat_3rd_Col    = [];       % 3rd column information if available: 'COV', 'STD' etc

Target = PreprocessTargetData(Target, Include_min_COV, Resample, Number, Stat_3rd_Col);

%% Optimization (Differential Evolution) parameters =====================================================================================
DE.Trial_num = [1:2];     % trials to run per layering ratio
DE.Maxitr    = 100;       % generations
DE.TolStall  = 1e-7;      % relative stall tol on best f
DE.StallGens = 40;        % stop if no progress for this many gens
DE.NP        = 50;        % NP ≈ (15-20) * number_of_variables
DE.CR        = 0.9;       % crossover probability
DE.pbestFrac = 0.15;      % fraction of population in the p-best elite pool (JADE current-to-pbest)
DE.Verbose   = true;      % print early-stop messages

%% Forward-model parameters =============================================================================================================
HTLM.d   = 5;       % polynomial order
HTLM.EPW = 4;       % elements per wavelength
HTLM.PML = 20;      % PML layers

%% Layer parameter settings =============================================================================================================
Layer_param.h                = 'variable';     % 'variable' or 'fixed'
Layer_param.h_range          = 'LN_based';     % If 'h' is 'variable': choose 'LR_based' (Layer Ratio) or 'LN_based' (Layer by Number)
Layer_param.fixed            = 'increasing';   % If 'h' is 'fixed': choose layering type:'increasing' or 'equal'
Layer_param.h_count          = [5];            % If 'h' is 'variable': [LRs: Ex. 1.5 2 3], [LNs: Ex. 4 5... + H.S]; If 'fixed': [LNs] as the number of layers
Layer_param.Vs               = 'variable';     % Always: 'variable'
Layer_param.Vs_range         = [100 800];      % numeric bounds:Ex.[150 800] or use: 'Auto'
Layer_param.nu               = 'variable';     % 'variable' or 'fixed'
Layer_param.nu_range         = [0.25 0.35];    % fixed at: 0.XX, or use : [0.25 0.35]; for variable range
Layer_param.rho              = 'fixed';        % 'variable' or 'fixed'
Layer_param.rho_range        = [2000 2100];    % fixed at : XXXX, or use :[1700 2000]; for variable range
Layer_param.Vs_reversal      = 0.5;            % allowable fractional Vs drop (if any)
Layer_param.Vs_rev_lim       = 20;             % Depth below which Vs reversal is not allowed
Layer_param.Vs_halfspace_max = 'yes';          % Half-space velocity always maximum, 'yes'/'no'
Layer_param.DCR_range   = [2];    % DCR:Depth conversion ratio, D_max = max_wavelength / DCR; Fixed at :[X] or provide range[2 3]

% Read description ---------------------------------------
%  Layer Thickness (h) Settings
% Layer_param.h
%   'variable' – Layer thickness values are optimized by DE.
%   'fixed'    – Layer thicknesses remain constant.
%
% Layer_param.h_range
%   If h = 'variable', this determines *how* the thickness varies:
%     'LN_based' – "Layer Number"-based control:
%                  the number of layers is specified in 'h_count',
%                  and their thicknesses are automatically adjusted
%                  to match the total model depth.
%     'LR_based' – "Layer Ratio"-based control: (Cox and Teague - 2016)
%                  thickness ratios between consecutive layers are
%                  defined in h_count (e.g., [1.5 2 3] means the
%                  lower layers are proportionally thicker).
%
% Layer_param.fixed
%   Used only when h = 'fixed':
%     'increasing' – Each layer thickness increases with depth.
%     'equal'      – All layers have equal thickness.
%
% Layer_param.h_count
%   Defines the number or ratios of layers depending on the mode:
%     - If h = 'variable' and h_range = 'LN_based':
%           h_count = [4 5 6]  → Try models with 4, 5, and 6 layers.
%     - If h = 'variable' and h_range = 'LR_based':
%           h_count = [1.5 2 3] → Define layer thickness ratios.
%     - If h = 'fixed':
%           h_count = [5] → Use a fixed 5-layer model.
%
% Example: 1
%   Layer_param.h       = 'variable';
%   Layer_param.h_range = 'LN_based';
%   Layer_param.h_count = [5]; % Supports multiple values: [5 6 7]
% Example: 2
%   Layer_param.h       = 'variable';
%   Layer_param.h_range = 'LR_based';
%   Layer_param.h_count = [1.5]; % Supports multiple values: [1.5 2 3]
% When multiple values are provided, the code runs in a loop and saves all
% corresponding results.
% -------------------------------------------------------------
% Shear-Wave Velocity (Vs) Settings
%
% Layer_param.Vs
%   'variable' – Vs values are optimized during the DE process.
%
% Layer_param.Vs_range
%   Defines the search range or constraint for Vs values:
%     - Numeric range (e.g., [150 800]) → Vs can vary within these bounds.
%     - 'Auto' → Automatically estimated from input data or depth trends.
%
% Layer_param.Vs_halfspace_max
%   Ensures the half-space (bottom layer) has the maximum Vs value.
%     'yes' → Force half-space Vs ≥ all other layers.
%     'no'  → Allow the half-space Vs to vary freely.
%
% The function, "PM_Inv_Objective.m" contain objective function. The user can
% define new objective function by their choice.
% -------------------------------------------------------------
% Note::'Vs_reversal' ensures that the velocity drop between layers is not too abrupt. 
% Example, Vs(1) = 400 m/s, Vs_reversal = 0.1. So, minimum allowable Vs(2) = 400-0.1*400 = 360 m/s
% If there is no indication of velocity reversal, keep this value small.
%% Prepare output folder ================================================================================================================
% Output folder based on fv_file stem
global name
[~, name, ~] = fileparts(Target.name);
outdir       = fullfile(scriptDir, name);   % outputs anchored to the script folder, not pwd
if ~exist(outdir,'dir'); mkdir(outdir); end

% ArchiveRun('PM_Inv_Input.m', outdir);

save(fullfile(outdir, 'input_params.mat'), 'Target','Layer_param','HTLM','DE');

%% Run inversion ========================================================================================================================
for ii = 1 : numel(Layer_param.h_count)

    h_count = Layer_param.h_count(ii);

    for trial = DE.Trial_num

        fprintf('Running %s Profile: %s | h_count=%.3g | trial=%d\n', Target.Measured_component, name, h_count, trial);

        % --- Initial model & bounds ---
        PlotTargetData(Target)
        [initial_param, lbinfo, ubinfo, N_layer] = GetInitialModel(ii, Target, Layer_param);
        save(fullfile(outdir, 'Temp_params.mat'), ...
            'Target','Layer_param','HTLM', 'initial_param', 'lbinfo', 'ubinfo','N_layer');

        % Upper and lower bound
        lb = lbinfo.all;
        ub = ubinfo.all;

        %%%%%%%%%%%% Differential Evolution (DE) %%%%%%%%%%%%

        [initPop] = GetInitialPopulation(lbinfo, ubinfo, Layer_param, DE.NP);

        % Attach the freshly generated initial population to the DE options struct
        DE.InitPop = initPop;

        global param ; param = cell(8, DE.Maxitr);

        % --- Solve ---
        tic
        [xbest, ~, info] = DE_Optimizer(lb, ub, N_layer, DE, HTLM,...
            Layer_param, initial_param, Target);
        toc

        itr = info.gens;

        % --- Save per-run output ---
        save(fullfile(outdir,[name,'_h_count_', num2str(h_count),'_trial_', num2str(trial),'.mat']),...
            'itr','N_layer','param','xbest','Target', 'Layer_param','lbinfo','ubinfo');

        close all

    end
end
