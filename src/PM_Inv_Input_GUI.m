function PM_Inv_Input_GUI
% % =====================================================================================================================================
% ============================= Differential Evolution (DE) for Predominant Mode Inversion (PMinv) ==============================
% ===================================================================================================================================== 
% It has both VCPM (Vertical component predominant mode) and RCPM (Radial
% component predominant mode) Inversion

%  INTERACTIVE INPUT (GUI front-end for PM_Inv_Input.m)
%  Author : [Mrinal Bhaumik]
%  Date   : [Jun_2026]
%  Version: 1.0
% ------------------------------------------------------------
%  A tabbed uifigure that collects the same parameters as PM_Inv_Input.m
%  (Target data, Layering, Model parameters, Constraints & objective,
%  Optimizer) and lets you:
%     - Preview Target : load + preprocess + plot the target curve
%     - Export Config  : push Target/DE/HTLM/Layer_param to the base workspace
%     - Run Inversion  : run the full DE inversion inline (same pipeline as PM_Inv_Input.m)
% ------------------------------------------------------------
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

%% Resolve paths and make the helper functions visible
scriptDir       = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(scriptDir, 'helpers')));

%% Verify the required MATLAB toolboxes are installed
checkRequiredToolboxes();

S               = struct();     % holds every input-control handle
S.scriptDir     = scriptDir;

%% ---------------------------------------------------------------- Main window
fig                  = uifigure('Name','PM_Inv — Interactive Input','Position',[100 100 780 660], 'HandleVisibility','off');   
mainGrid             = uigridlayout(fig, [3 1]);
mainGrid.RowHeight   = {42, '1x', 96};
mainGrid.ColumnWidth = {'1x'};
hdr                  = uilabel(mainGrid, 'Text','Surface-Wave Inversion — Input Setup', ...
    'FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center');
hdr.Layout.Row  = 1;

tg              = uitabgroup(mainGrid);
tg.Layout.Row   = 2;

%% ---------------------------------------------------------------- Tab 1: Target data
t1              = uitab(tg, 'Title','Target Data');
g1              = tabGrid(t1, 6); % creating six buttons
g1.ColumnWidth  = {210, '1x', 90};
addNote(g1, 1, 'Target dispersion (phase velocity vs frequency). Point your cursor on the field for a description.');
S.fv            = addField(g1, 2, 'Dispersion file (fv)', 'text', fullfile(scriptDir,'50m_R0_dispersion_stats.txt'), {}, ...
    ['Path to the target dispersion file (.txt/.csv/.xlsx). Columns: 1 = frequency (Hz); ' ...
     '2 = phase velocity (m/s) or slowness (s/m, auto-detected and inverted); ' ...
     '3 = optional uncertainty (define it in "3rd-column meaning")']);
addBrowse(g1, 2, @() onBrowse('fv'));
S.stat3         = addField(g1, 3, '3rd-column meaning', 'drop', '(none)', {'(none)','COV','STD'}, ...
    ['Define the OPTIONAL 3rd column: "COV" = coef. of variation; "STD" = standard dev. ' ...
     'With "(none)" the 3rd column is taken as the std of velocity as-is, so if your file actually holds a COV, ' ...
     'leaving this at "(none)" makes the weights (1/std^2) blow up (huge misfit)']);
S.component     = addField(g1, 4, 'Measured_component', 'drop','RCPM',{'VCPM','RCPM'}, ...
    ['"VCPM" = vertical-component predominant mode; ' ...
     '"RCPM" = radial-component predominant mode.']);
S.cov           = addField(g1, 5, 'Include_min_COV', 'text', '', {}, ...
    ['Used ONLY when the file has no 3rd (uncertainty) column: builds a std column as ' ...
     '(phase velocity x this COV), e.g. 0.05 or 0.1. Leave empty to skip.']);
S.resmp         = addField(g1, 6, 'Resample (log-scale)', 'drop', 'no', {'no','yes'}, ...
    ['Resample the target onto a logarithmic frequency scale before inversion. "no" keeps the original ' ...
     'sampling; "yes" resamples to "Number of log samples" points.']);
S.number        = addField(g1, 7, 'Number of log samples', 'text', '', {}, ...
    'Number of log-spaced frequency points when Resample = yes (e.g. 40). Ignored when Resample = no; defaults to 40 if left empty.');

%% ---------------------------------------------------------------- Tab 2: Layering + Model parameters
t2              = uitab(tg, 'Title','Layering & Model');
tg2             = uigridlayout(t2, [1 2]);
tg2.ColumnWidth = {'1x','1x'};
tg2.Padding     = [8 8 8 8];

gL              = groupGrid(tg2, 1, 'Layering (h)', ...
    'Point your cursor on the field for a description.', 4);
S.h             = addField(gL, 2, 'Thickness mode (h)',    'drop', 'variable',   {'variable','fixed'}, ...
    'How layer thicknesses are handled ? "variable": thicknesses are inverted. "fixed": thicknesses are held constant — h_count sets the number of layers and "Layering (if fixed)" sets their pattern.');
S.h_range       = addField(gL, 3, 'h_range (if variable)', 'drop', 'LN_based',   {'LN_based','LR_based'}, ...
    'Used only when h = variable. "LN_based" (Layer-Number): h_count = number of layers; thicknesses auto-adjust to the model depth. "LR_based" (Layer-Ratio, Cox & Teague 2016): h_count = thickness ratio between consecutive layers.');
S.h_fixed       = addField(gL, 4, 'Layering (if fixed)',   'drop', 'increasing', {'increasing','equal'}, ...
    'Used only when h = fixed. "increasing": each layer is thicker than the one above. "equal": all layers have equal thickness.');
S.h_count       = addField(gL, 5, 'h_count',               'text', '5', {}, ...
    'Batch input — accepts multiple values and runs each as a separate case in a loop (e.g. "5 6 7", "5:7", or "5,6,8"). Meaning depends on the mode: LN_based = number of layers (+ half-space); LR_based = layer-thickness ratio(s), e.g. 1.5 2 3; fixed = number of layers.');

gR              = groupGrid(tg2, 2, 'Model Parameters', ...
    'Point your cursor on the field for a description. Ranges:: two values "min max" when variable; a single value when fixed.', 6);
S.vs            = addField(gR, 2, 'Vs mode',               'drop', 'variable', {'variable'}, ...
    'Shear-wave velocity is always inverted ("variable").');
S.vs_range      = addField(gR, 3, 'Vs_range',              'text', '100 800', {}, ...
    'Search bounds for Vs (m/s): two numbers "min max" (e.g. 100 800), or "Auto" to estimate the bounds from the target dispersion curve.');
S.nu            = addField(gR, 4, 'nu mode',               'drop', 'variable', {'variable','fixed'}, ...
    'Poisson''s ratio. "variable": nu is inverted — give a two-value range in nu_range. "fixed": nu is held constant — give a single value in nu_range.');
S.nu_range      = addField(gR, 5, 'nu_range',              'text', '0.25 0.35', {}, ...
    'If nu = variable: two values "min max" (e.g. 0.25 0.35). If nu = fixed: a single value (e.g. 0.3).');
S.rho           = addField(gR, 6, 'rho mode',              'drop', 'fixed',    {'variable','fixed'}, ...
    'Mass density. "variable": rho is inverted — give a two-value range. "fixed": rho is held constant — give a single value.');
S.rho_range     = addField(gR, 7, 'rho_range',             'text', '2000', {}, ...
    'If rho = variable: two values "min max" (e.g. 1700 2100). If rho = fixed: a single value (e.g. 2000). Units: kg/m^3.');

%% ---------------------------------------------------------------- Tab 3: Constraints/Objective + Optimizer
t3              = uitab(tg, 'Title','Constraints & Optimizer');
tg3             = uigridlayout(t3, [1 2]);
tg3.ColumnWidth = {'1x','1x'};
tg3.Padding     = [8 8 8 8];

gC              = groupGrid(tg3, 1, 'Constraints', ...
    'Point your cursor on field for a description.', 4);
S.vs_rev        = addField(gC, 2, 'Vs_reversal',           'num',  0.5, {}, ...
    'Maximum allowed fractional Vs drop between consecutive layers. e.g. 0.5 rejects any model where a layer''s Vs falls more than 50% below the layer above it (velocity-inversion control).');
S.vs_rev_lim    = addField(gC, 3, 'Vs_rev_lim (depth, m)', 'num',  20, {}, ...
    'Depth (m) below which Vs reversals are not allowed. Low-velocity zones are tolerated only shallower than this depth; deeper reversals are repaired.');
S.vs_hs         = addField(gC, 4, 'Vs_halfspace_max',      'drop', 'yes', {'yes','no'}, ...
    'If "yes", the half-space (bottom) Vs is forced to be the maximum of the profile, enforcing increasing velocity at the base.');
S.dcr           = addField(gC, 5, 'DCR_range',             'text', '2', {}, ...
    'Depth-conversion ratio (DCR): d_max_limit = max(wavelength) / DCR. Enter one value (e.g. 2) to fix it, or "min max" (e.g. 2 3) to randomise DCR per initialisation. Smaller DCR = deeper model.');

gO              = groupGrid(tg3, 2, 'Optimizer (DE + HTLM)', ...
    'Point your cursor on field for a description. DE controls (top) + forward-model / HTLM settings (bottom).', 11);
S.trial         = addField(gO, 2,  'DE.Trial_num',              'text',  '1:2', {}, ...
    'Independent DE trials to run per layering (e.g. 1:2 = 2 trials, 1:5 = 5). Each trial is a fresh run with a new random population, saved as its own .mat file.');
S.maxitr        = addField(gO, 3,  'DE.Maxitr',                 'num',   100, {}, ...
    'Maximum number of DE generations (iterations) per trial.');
S.np            = addField(gO, 4,  'DE.NP',                     'num',   50, {}, ...
    'Population size (candidate models per generation). Rule of thumb: NP is about 15-20 x the number of inverted variables.');
S.cr            = addField(gO, 5,  'DE.CR',                     'num',   0.9, {}, ...
    'Crossover probability (0-1): fraction of vector components taken from the mutant. Higher = more mixing / exploration.');
S.pbest         = addField(gO, 6,  'DE.pbestFrac',              'num',   0.15, {}, ...
    'JADE current-to-pbest fraction (0-1): the top pbestFrac of the population forms the elite "p-best" pool that mutations steer toward. Smaller = greedier/faster convergence, larger = more exploratory.');
S.tolstall      = addField(gO, 7,  'DE.TolStall',               'num',   1e-7, {}, ...
    'Relative-improvement tolerance for the stall test. A generation counts as "no progress" when the best misfit improves by less than this.');
S.stallgens     = addField(gO, 8,  'DE.StallGens',              'num',   40, {}, ...
    'Early-stop patience: stop the trial if the best misfit has not improved for this many consecutive generations.');
S.verbose       = addField(gO, 9,  'DE.Verbose',                'check', true, {}, ...
    'Print per-generation progress and early-stop messages to the MATLAB command window.');
S.htlm_d        = addField(gO, 10, 'HTLM.d (poly order)',       'num',   5, {}, ...
    'Polynomial order of the thin-layer element shape functions (forward model). Higher = more accurate but slower.');
S.htlm_epw      = addField(gO, 11, 'HTLM.EPW (elem/wavelength)','num',   4, {}, ...
    'Elements per wavelength used to discretize each layer in the forward model. Higher = finer mesh, more accurate, slower.');
S.htlm_pml      = addField(gO, 12, 'HTLM.PML (layers)',         'num',   20, {}, ...
    'Number of PML (perfectly matched layer) elements forming the absorbing bottom half-space boundary.');

%% ---------------------------------------------------------------- Tab 4: Layer limits
tLim            = uitab(tg, 'Title','Layer Limits');
gLim            = uigridlayout(tLim, [3 1]);
gLim.RowHeight   = {36, 26, '1x'};
gLim.ColumnWidth = {'1x'};

barLim          = uigridlayout(gLim, [1 2]);
barLim.Layout.Row   = 1;
barLim.RowHeight    = {'1x'};
barLim.ColumnWidth  = {230, '1x'};
S.computeLimBtn = uibutton(barLim, 'Text','Compute / Refresh limits', 'ButtonPushedFcn', @(~,~) onComputeLimits());
S.computeLimBtn.Layout.Column = 1;
lblLim          = uilabel(barLim, 'FontAngle','italic', 'FontColor',[0.35 0.35 0.35], 'WordWrap','on', ...
    'Text','Computes per-layer lower/upper bounds from the current settings. Edit any cell; Run/Preview then uses the edited limits (while the layer count still matches).');
lblLim.Layout.Column = 2;

lblLim2         = uilabel(gLim, 'FontColor',[0.3 0.3 0.3], ...
    'Text','Rows = layers (L1 = shallowest … HS = half-space). Two columns per parameter: lower and upper bound.');
lblLim2.Layout.Row = 2;

S.limTable      = uitable(gLim);
S.limTable.Layout.Row = 3;

%% ---------------------------------------------------------------- Tab 5: Figures & control
tFig            = uitab(tg, 'Title','Figures & Control');
g6              = uigridlayout(tFig, [2 1]);
g6.RowHeight    = {40, '1x'};
g6.ColumnWidth  = {'1x'};

% --- Control bar ---
ctrlBar         = uigridlayout(g6, [1 3]);
ctrlBar.Layout.Row  = 1;
ctrlBar.RowHeight   = {'1x'};
ctrlBar.ColumnWidth = {110, 110, '1x'};

S.pauseBtn      = uibutton(ctrlBar, 'Text','Pause', 'Enable','off', 'ButtonPushedFcn', @(~,~) onPause());
S.pauseBtn.Layout.Column = 1;
S.stopBtn       = uibutton(ctrlBar, 'Text','Stop', 'Enable','off', 'FontColor',[0.6 0 0], ...
    'ButtonPushedFcn', @(~,~) onStop());
S.stopBtn.Layout.Column = 2;
S.figNote       = uilabel(ctrlBar, 'FontAngle','italic', 'FontColor',[0.35 0.35 0.35], 'WordWrap','on', ...
    'Text','Columns: phase velocity | Vs-depth | Vs-depth zoom | convergence.');
S.figNote.Layout.Column = 3;

% --- Dashboard host: axes are (re)built per run to match the number of targets ---
S.dashHost      = uipanel(g6, 'BorderType','none');
S.dashHost.Layout.Row = 2;
buildDashboard();   % initial dashboard

%% ---------------------------------------------------------------- Tab 6: Post-processing
tPP             = uitab(tg, 'Title','Post-Processing');
gPP             = uigridlayout(tPP, [2 1]);
gPP.RowHeight   = {40, '1x'};
gPP.ColumnWidth = {'1x'};

ppBar           = uigridlayout(gPP, [1 11]);
ppBar.Layout.Row   = 1;
ppBar.RowHeight    = {'1x'};
ppBar.ColumnWidth  = {70, 60, 45, 110, 50, 55, 55, 50, 130, 80, '1x'};

uilabel(ppBar, 'Text','Dmax (m)', 'HorizontalAlignment','right');
S.pp_dmax       = uieditfield(ppBar, 'numeric', 'Value',30, 'Limits',[0 Inf]);
uilabel(ppBar, 'Text','Mode', 'HorizontalAlignment','right');
S.pp_mode       = uidropdown(ppBar, 'Items',{'Best','All','MeanBest'}, 'Value','All', ...
    'Tooltip', ['Best: single lowest-misfit model per h_count.  ' ...
                'All: the Nbest lowest-misfit models from each trial file.  ' ...
                'MeanBest: mean of the best-per-trial models for each h_count.']);
uilabel(ppBar, 'Text','Nbest', 'HorizontalAlignment','right');
S.pp_nbest      = uieditfield(ppBar, 'numeric', 'Value',10, 'Limits',[1 Inf], 'RoundFractionalValues','on', ...
    'Tooltip','Models drawn per trial (used only when Mode = All).');
uilabel(ppBar, 'Text','Modes', 'HorizontalAlignment','right');
S.pp_nmodes     = uieditfield(ppBar, 'numeric', 'Value',3, 'Limits',[0 5], 'RoundFractionalValues','on', ...
    'Tooltip','Number of theoretical modes overlaid in the dispersion panel (0-5).');
S.pp_selBtn     = uibutton(ppBar, 'Text','Select .mat…', 'ButtonPushedFcn', @(~,~) onSelectResults());
S.pp_plotBtn    = uibutton(ppBar, 'Text','Plot', 'FontWeight','bold', ...
    'BackgroundColor',[0.20 0.45 0.70], 'FontColor','w', 'ButtonPushedFcn', @(~,~) onPlotResults());
S.pp_status     = uilabel(ppBar, 'Text','Select saved result .mat files, then Plot.', ...
    'FontAngle','italic', 'FontColor',[0.35 0.35 0.35]);

ppAx            = uigridlayout(gPP, [1 2]);
ppAx.Layout.Row = 2;
ppAx.ColumnWidth = {'1x','1x'};
S.pp_axDisp     = uiaxes(ppAx); S.pp_axDisp.Layout.Column = 1; title(S.pp_axDisp,'Dispersion');
S.pp_axProf     = uiaxes(ppAx); S.pp_axProf.Layout.Column = 2; title(S.pp_axProf,'Vs profiles');
S.pp_files      = {};   S.pp_path = '';

%% ---------------------------------------------------------------- Bottom bar: actions + status
bottom          = uigridlayout(mainGrid, [2 3]);
bottom.Layout.Row   = 3;
bottom.RowHeight    = {36, '1x'};
bottom.ColumnWidth  = {'1x','1x','1x'};

S.previewBtn    = uibutton(bottom, 'Text','Preview Target', 'ButtonPushedFcn', @(~,~) onPreview());
S.previewBtn.Layout.Row = 1; S.previewBtn.Layout.Column = 1;

S.exportBtn     = uibutton(bottom, 'Text','Export Config', 'ButtonPushedFcn', @(~,~) onExport());
S.exportBtn.Layout.Row = 1; S.exportBtn.Layout.Column = 2;

S.runBtn        = uibutton(bottom, 'Text','Run Inversion', 'FontWeight','bold', ...
    'BackgroundColor',[0.20 0.55 0.30], 'FontColor','w', 'ButtonPushedFcn', @(~,~) onRun());
S.runBtn.Layout.Row = 1; S.runBtn.Layout.Column = 3;

S.status        = uilabel(bottom, 'Text','Ready.', 'FontColor',[0.2 0.2 0.2]);
S.status.Layout.Row = 2; S.status.Layout.Column = [1 3];

% ================================================================ Callbacks (nested)

    function onBrowse(~)
        [f, p] = uigetfile({'*.txt;*.csv;*.xlsx','Data files (*.txt,*.csv,*.xlsx)'; '*.*','All files'});
        if isequal(f, 0); return; end
        S.fv.Value = fullfile(p, f);
    end

    function onPreview()
        global WMI_RUN %#ok<GVMIS>
        tg.SelectedTab = tFig;
        lockUI(true);
        try
            setStatus('Loading & preprocessing target…');
            C = gatherInputs();
            buildDashboard();
            WMI_RUN = struct('pause',false, 'stop',false, 'convAx',S.convAx, 'dash',S.dash);
            Target = PreprocessTargetData(C.Target, C.Include_min_COV, C.Resample, C.Number, C.Stat_3rd_Col);
            PlotTargetData(Target);
            setStatus('Target preview plotted in the Figures tab.');
        catch ME
            setStatus(['Error: ' ME.message]);
            uialert(fig, ME.message, 'Preview error');
        end
        lockUI(false);
    end

    function onExport()
        try
            C = gatherInputs();
            assignin('base','Target',          C.Target);
            assignin('base','DE',              C.DE);
            assignin('base','HTLM',            C.HTLM);
            assignin('base','Layer_param',     C.Layer_param);
            assignin('base','Include_min_COV', C.Include_min_COV);
            assignin('base','Resample',        C.Resample);
            assignin('base','Number',          C.Number);
            assignin('base','Stat_3rd_Col',    C.Stat_3rd_Col);
            setStatus('Config exported to base workspace (Target, DE, HTLM, Layer_param, …).');
        catch ME
            setStatus(['Error: ' ME.message]);
            uialert(fig, ME.message, 'Export error');
        end
    end

    function onRun()
        global WMI_RUN %#ok<GVMIS>
        tg.SelectedTab = tFig;
        lockUI(true);
        try
            setStatus('Gathering inputs…');
            C = gatherInputs();
            buildDashboard();
            WMI_RUN = struct('pause',false, 'stop',false, 'convAx',S.convAx, 'dash',S.dash);
            S.pauseBtn.Text = 'Pause';
            runInversion(C);
            if WMI_RUN.stop
                setStatus('Inversion stopped by user.');
            else
                setStatus('Inversion complete.');
            end
        catch ME
            setStatus(['Error: ' ME.message]);
            uialert(fig, ME.message, 'Inversion error');
        end
        lockUI(false);
    end

    function onPause()
        global WMI_RUN %#ok<GVMIS>
        if isempty(WMI_RUN); return; end
        WMI_RUN.pause = ~WMI_RUN.pause;
        if WMI_RUN.pause
            S.pauseBtn.Text = 'Resume';
            setStatus('Paused (waits after the current generation).');
        else
            S.pauseBtn.Text = 'Pause';
            setStatus('Resumed.');
        end
    end

    function onStop()
        global WMI_RUN %#ok<GVMIS>
        if isempty(WMI_RUN); return; end
        WMI_RUN.stop  = true;
        WMI_RUN.pause = false;
        setStatus('Stopping after the current generation…');
    end

    function onComputeLimits()
        global WMI_RUN %#ok<GVMIS>
        savedRun = WMI_RUN;
        lockUI(true);
        try
            setStatus('Computing layer limits…');
            C = gatherInputs();
            Target = PreprocessTargetData(C.Target, C.Include_min_COV, C.Resample, C.Number, C.Stat_3rd_Col);

            WMI_RUN = struct('silent', true);              % suppress GetInitialModel's limit plots
            [~, lbinfo, ubinfo, N_layer] = GetInitialModel(1, Target, C.Layer_param);
            WMI_RUN = savedRun;                            % restore

            canon  = {'h','vs','nu','rho'};
            labels = struct('h','h', 'vs','Vs', 'nu','nu', 'rho','rho');
            nRow   = N_layer + 1;
            params = {}; colNames = {}; Data = [];
            for k = 1:numel(canon)
                p = canon{k};
                if ~isfield(lbinfo, p); continue; end
                lo = lbinfo.(p)(:); up = ubinfo.(p)(:);
                if numel(lo) < nRow                        % pad h (finite layers only) to include the HS row
                    lo(end+1:nRow,1) = NaN; up(end+1:nRow,1) = NaN;
                end
                params{end+1}   = p;                          %#ok<AGROW>
                Data            = [Data, lo, up];             %#ok<AGROW>
                colNames{end+1} = [labels.(p) ' (low)'];      %#ok<AGROW>
                colNames{end+1} = [labels.(p) ' (up)'];       %#ok<AGROW>
            end
            rowNames = arrayfun(@(i) sprintf('L%d', i), 1:N_layer, 'UniformOutput', false);
            rowNames{end+1} = 'HS';

            S.limTable.Data           = Data;
            S.limTable.ColumnName     = colNames;
            S.limTable.RowName        = rowNames;
            S.limTable.ColumnEditable = true(1, size(Data,2));
            S.limits = struct('active',true, 'N_layer',N_layer, 'params',{params});
            setStatus(sprintf('Limits computed for %d layers (+ half-space). Edit cells as needed; Run/Preview will use them.', N_layer));
        catch ME
            WMI_RUN = savedRun;
            setStatus(['Error: ' ME.message]);
            uialert(fig, ME.message, 'Layer-limits error');
        end
        lockUI(false);
    end

    function [lbinfo, ubinfo] = applyLayerLimits(lbinfo, ubinfo, N_layer)
        % Override the computed bounds with the GUI-edited table (order-preserving),
        % but only when limits were computed and the layer count still matches.
        if ~isfield(S,'limits') || isempty(S.limits) || ~S.limits.active || S.limits.N_layer ~= N_layer
            return;
        end
        Data   = S.limTable.Data;
        params = S.limits.params;
        col = 1; allLb = []; allUb = [];
        for k = 1:numel(params)
            p  = params{k};
            lo = Data(:, col); up = Data(:, col+1); col = col + 2;
            if strcmp(p, 'h'); lo = lo(1:N_layer); up = up(1:N_layer); end
            lbinfo.(p) = lo(:).'; ubinfo.(p) = up(:).';
            allLb = [allLb lbinfo.(p)]; allUb = [allUb ubinfo.(p)]; %#ok<AGROW>
        end
        lbinfo.all = allLb; ubinfo.all = allUb;
    end

    % ---- Assemble all GUI values into the config structs ----
    function C = gatherInputs()
        C.Target.fv_file = strtrim(S.fv.Value);
        C.Target.Measured_component = S.component.Value;

        C.Include_min_COV = parseVec(S.cov.Value);
        C.Resample        = S.resmp.Value;
        C.Number          = parseNumOrEmpty(S.number.Value);
        st = S.stat3.Value;
        if strcmp(st,'(none)'); C.Stat_3rd_Col = []; else; C.Stat_3rd_Col = st; end

        C.DE.Trial_num = parseVec(S.trial.Value);
        C.DE.Maxitr    = S.maxitr.Value;
        C.DE.TolStall  = S.tolstall.Value;
        C.DE.StallGens = S.stallgens.Value;
        C.DE.NP        = S.np.Value;
        C.DE.CR        = S.cr.Value;
        C.DE.pbestFrac = S.pbest.Value;
        C.DE.Verbose   = logical(S.verbose.Value);

        C.HTLM.d   = S.htlm_d.Value;
        C.HTLM.EPW = S.htlm_epw.Value;
        C.HTLM.PML = S.htlm_pml.Value;

        C.Layer_param.h                = S.h.Value;
        C.Layer_param.h_range          = S.h_range.Value;
        C.Layer_param.fixed            = S.h_fixed.Value;
        C.Layer_param.h_count          = parseVec(S.h_count.Value);
        C.Layer_param.Vs               = S.vs.Value;
        C.Layer_param.Vs_range         = parseRangeOrKeyword(S.vs_range.Value);
        C.Layer_param.nu               = S.nu.Value;
        C.Layer_param.nu_range         = parseVec(S.nu_range.Value);
        C.Layer_param.rho              = S.rho.Value;
        C.Layer_param.rho_range        = parseVec(S.rho_range.Value);
        C.Layer_param.Vs_reversal      = S.vs_rev.Value;
        C.Layer_param.Vs_rev_lim       = S.vs_rev_lim.Value;
        C.Layer_param.Vs_halfspace_max = S.vs_hs.Value;
        C.Layer_param.DCR_range = parseVec(S.dcr.Value);
    end

    % ---- Run the DE inversion (mirrors PM_Inv_Input.m) ----
    function runInversion(C)
        global name param WMI_RUN %#ok<GVMIS>   % name/param for helpers; WMI_RUN (set by onRun) for GUI pause/stop

        % Ensure a parallel pool exists (the DE solver uses parfor)
        if license('test','Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
            setStatus('Starting parallel pool…'); parpool;
        end

        setStatus('Preprocessing target data…');
        Target      = PreprocessTargetData(C.Target, C.Include_min_COV, C.Resample, C.Number, C.Stat_3rd_Col);
        Layer_param = C.Layer_param;
        HTLM        = C.HTLM;
        DE          = C.DE;

        % Output folder based on fv_file stem
        [~, name, ~] = fileparts(Target.name);
        outdir = fullfile(S.scriptDir, name);
        if ~exist(outdir,'dir'); mkdir(outdir); end
        save(fullfile(outdir, 'input_params.mat'), 'Target','Layer_param','HTLM','DE');

        for ii = 1 : numel(Layer_param.h_count)
            h_count = Layer_param.h_count(ii);
            for trial = DE.Trial_num
                if ~isempty(WMI_RUN) && WMI_RUN.stop; return; end   % user pressed Stop
                setStatus(sprintf('Running %s Profile: %s | h_count=%.3g | trial=%d …', Target.Measured_component, name, h_count, trial));
                fprintf('Running %s Profile: %s | h_count=%.3g | trial=%d\n', Target.Measured_component, name, h_count, trial);

                PlotTargetData(Target);
                [initial_param, lbinfo, ubinfo, N_layer, d_max_limit] = GetInitialModel(ii, Target, Layer_param);
                [lbinfo, ubinfo] = applyLayerLimits(lbinfo, ubinfo, N_layer);   % apply GUI-edited limits if present
                PlotLayerLimits(lbinfo, ubinfo, Layer_param, d_max_limit);      % plot the (edited) limit envelope
                save(fullfile(outdir, 'Temp_params.mat'), ...
                    'Target','Layer_param','HTLM','initial_param','lbinfo','ubinfo','N_layer');

                lb = lbinfo.all;
                ub = ubinfo.all;

                initPop    = GetInitialPopulation(lbinfo, ubinfo, Layer_param, DE.NP);
                DE.InitPop = initPop;

                param = cell(8, DE.Maxitr);

                tic
                [xbest, ~, info] = DE_Optimizer(lb, ub, N_layer, DE, HTLM, Layer_param, initial_param, Target);
                toc

                itr = info.gens;
                save(fullfile(outdir, [name,'_h_count_',num2str(h_count),'_trial_',num2str(trial),'.mat']), ...
                    'itr','N_layer','param','xbest','Target','Layer_param','lbinfo','ubinfo');

                close all   % closes analysis figures (GUI is HandleVisibility 'off', so it survives)
            end
        end
    end

    % ---- Small UI helpers (nested: need S/fig) ----
    function setStatus(msg)
        S.status.Text = msg; drawnow;
    end

    function lockUI(running)
        if running; s = 'off'; else; s = 'on'; end
        S.runBtn.Enable = s; S.exportBtn.Enable = s; S.previewBtn.Enable = s;
        if running; s2 = 'on'; else; s2 = 'off'; end   % Pause/Stop are active only while running
        S.pauseBtn.Enable = s2; S.stopBtn.Enable = s2;
        drawnow;
    end

    function onSelectResults()
        [f, p] = uigetfile('*.mat', 'Select PM_Inv result MAT-files', 'MultiSelect','on');
        if isequal(f, 0); return; end
        if ischar(f); f = {f}; end
        S.pp_files = f; S.pp_path = p;
        S.pp_status.Text = sprintf('%d file(s) selected. Click Plot.', numel(f));
    end

    function onPlotResults()
        if isempty(S.pp_files)
            S.pp_status.Text = 'Select .mat files first.'; return;
        end
        o = struct('Dmax',   S.pp_dmax.Value,   'Mode',  S.pp_mode.Value, ...
                   'Nbest',  S.pp_nbest.Value,  'nModes', S.pp_nmodes.Value);
        S.pp_status.Text = 'Plotting…'; drawnow;
        try
            PM_Inv_PlotResults(S.pp_files, S.pp_path, o, S.pp_axDisp, S.pp_axProf);
            S.pp_status.Text = sprintf('Plotted %d file(s)  [Mode: %s].', numel(S.pp_files), S.pp_mode.Value);
        catch ME
            S.pp_status.Text = ['Error: ' ME.message];
            uialert(fig, ME.message, 'Post-processing error');
        end
    end

    function buildDashboard()
        % (Re)create the dashboard axes:
        %   phase velocity | Vs-depth | Vs-depth zoom | convergence
        delete(S.dashHost.Children);
        gd = uigridlayout(S.dashHost, [2 4]);
        gd.RowHeight   = {'1x','1x'};
        gd.ColumnWidth = repmat({'1x'}, 1, 4);
        gd.Padding     = [6 6 6 6];

        D = struct();
        D.disp   = mkAxes(gd, [1 2], 1);      % phase velocity - frequency
        D.vs     = mkAxes(gd, [1 2], 2);      % Vs - depth
        D.vszoom = mkAxes(gd, [1 2], 3);      % Vs - depth (zoom)
        cax      = mkAxes(gd, [1 2], 4);      % convergence

        title(cax, 'DE Convergence'); xlabel(cax,'Generation'); ylabel(cax,'Best misfit');
        S.dash   = D;
        S.convAx = cax;
    end

end  % ===================================================================== main

% ---------------------------------------------------------------- Local helpers

function checkRequiredToolboxes()
% Verify the MATLAB release and toolboxes this GUI / inversion needs, and point
% the user to the install page for anything missing. Stops launch if a REQUIRED
% toolbox (or too-old MATLAB) is found; only warns for recommended ones.

    % --- MATLAB release (uifigure / uigridlayout need R2018b) ---
    if verLessThan('matlab','9.5')
        m = sprintf(['PM_Inv GUI needs MATLAB R2018b or newer (uigridlayout).\n' ...
                     'You are running R%s.'], version('-release'));
        errordlg(m, 'PM_Inv — MATLAB too old', 'modal');
        error('PM_Inv_Input_GUI:oldMATLAB', '%s', m);
    end

    % --- Toolboxes: {display name, a function it provides, install URL, mandatory?} ---
    reqs = {
        'Statistics and Machine Learning Toolbox', 'unifrnd', ...
            'https://www.mathworks.com/products/statistics.html',         true
        'Parallel Computing Toolbox',              'parpool', ...
            'https://www.mathworks.com/products/parallel-computing.html', false
        };

    lines = {}; hardMissing = false;
    for k = 1:size(reqs,1)
        if exist(reqs{k,2}, 'file') == 0        % toolbox function not found -> not installed
            if reqs{k,4}; tag = 'REQUIRED'; hardMissing = true; else; tag = 'recommended'; end
            lines{end+1} = sprintf('  - %s  (%s)\n        install: %s', ...
                reqs{k,1}, tag, reqs{k,3}); %#ok<AGROW>
        end
    end
    if isempty(lines); return; end

    msg = sprintf(['The following MATLAB toolbox(es) are not installed:\n\n%s\n\n' ...
                   'Use  Home > Add-Ons > Get Add-Ons,  or open the link(s) above.'], ...
                   strjoin(lines, '\n\n'));
    fprintf(2, '\n[PM_Inv GUI] %s\n\n', msg);   % Command Window (links are clickable there)

    if hardMissing
        errordlg(msg, 'PM_Inv — missing required toolbox', 'modal');
        error('PM_Inv_Input_GUI:missingToolbox', ...
              'A required toolbox is not installed (see the Command Window message).');
    else
        warndlg(msg, 'PM_Inv — optional toolbox missing');
    end
end

function g = tabGrid(parentTab, nFields)
% A 2-column (label | control) grid with one note row on top and a spacer at the bottom.
    g = uigridlayout(parentTab, [nFields+2, 2]);
    g.RowHeight   = [{28}, repmat({30}, 1, nFields), {'1x'}];
    g.ColumnWidth = {210, '1x'};
    g.Scrollable  = 'on';
    g.Padding     = [12 12 12 12];
end

function g = groupGrid(parentGrid, col, titleText, noteText, nFields)
% A titled group panel placed in column `col` of a parent grid, holding a
% [label | control] grid (note row on top, spacer at bottom). Its row layout
% matches tabGrid so the same addField row numbers (2..nFields+1) can be reused.
    p = uipanel(parentGrid, 'Title', titleText, 'FontWeight','bold');
    p.Layout.Row = 1; p.Layout.Column = col;
    g = uigridlayout(p, [nFields+2, 2]);
    g.RowHeight   = [{26}, repmat({30}, 1, nFields), {'1x'}];
    g.ColumnWidth = {185, '1x'};
    g.Scrollable  = 'on';
    g.Padding     = [8 8 8 8];
    addNote(g, 1, noteText);
end

function addNote(g, row, txt)
    lbl = uilabel(g, 'Text', txt, 'FontAngle','italic', 'FontColor',[0.35 0.35 0.35], ...
        'WordWrap','on');
    lbl.Layout.Row = row; lbl.Layout.Column = [1 2];
end

function ctrl = addField(g, row, labelTxt, kind, def, items, tip)
    if nargin < 6; items = {}; end
    if nargin < 7; tip = ''; end          % optional hover description (tooltip)
    lbl = uilabel(g, 'Text', labelTxt);
    lbl.Layout.Row = row; lbl.Layout.Column = 1;
    switch kind
        case 'text'
            ctrl = uieditfield(g, 'text', 'Value', def);
        case 'num'
            ctrl = uieditfield(g, 'numeric', 'Value', def);
        case 'drop'
            ctrl = uidropdown(g, 'Items', items, 'Value', def);
        case 'check'
            ctrl = uicheckbox(g, 'Text', 'enabled', 'Value', def);
    end
    ctrl.Layout.Row = row; ctrl.Layout.Column = 2;
    if ~isempty(tip)
        lbl.Tooltip  = tip;               % description shows on hover (label or control)
        ctrl.Tooltip = tip;
    end
end

function addBrowse(g, row, cb)
    b = uibutton(g, 'Text','Browse…', 'ButtonPushedFcn', @(~,~) cb());
    b.Layout.Row = row; b.Layout.Column = 3;
end

function ax = mkAxes(parent, r, c)
    ax = uiaxes(parent);
    ax.Layout.Row = r; ax.Layout.Column = c;
end

% ---- Parsing helpers ----
function v = parseVec(txt)
    txt = strtrim(txt);
    if isempty(txt); v = []; else; v = str2num(txt); end %#ok<ST2NM>
end

function v = parseNumOrEmpty(txt)
    txt = strtrim(txt);
    if isempty(txt); v = []; else; v = str2double(txt); end
end

function v = parseRangeOrKeyword(txt)
    txt = strtrim(txt);
    num = str2num(txt); %#ok<ST2NM>
    if isempty(num); v = txt; else; v = num; end
end
