% =====================================================================================================================================
% ============================= PM_Inv Post-Processing : plot saved inversion results ==============================
% =====================================================================================================================================
%  Author : [Mrinal Bhaumik]
%  Date   : [Jun_2025]
%  Version: 1.0
% ------------------------------------------------------------
%  Description:
%  Standalone post-processor for PM_Inv_Input.m results. Prompts for the saved
%  per-run .mat files, then draws two panels: the dispersion curves (first modes
%  + predominant mode vs the target) and the Vs-depth profiles. Choose how the
%  ensemble is reduced with the options below.
% ------------------------------------------------------------

clc; clear; close all

% -------- options --------
opts.Dmax   = 30;        % maximum plotting depth (m)
opts.Mode   = 'Best';    % 'All' | 'Best' | 'MeanBest'
opts.Nbest  = 1;         % models per trial (used only when Mode = 'All')
opts.nModes = 3;         % number of theoretical modes to overlay in the dispersion panel
%
%   'All'      : the Nbest lowest-misfit models from EACH selected trial file.
%   'Best'     : the single lowest-misfit model per h_count (across its trials).
%   'MeanBest' : the mean of the best-per-trial models for each h_count.

% Make the helper functions visible (relative to this script)
scriptDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(scriptDir, 'helpers')));

% -------- select the saved result files --------
[files, path] = uigetfile('*.mat', 'Select PM_Inv result MAT-files', 'MultiSelect', 'on');
if isequal(files, 0); disp('No files selected.'); return; end

% -------- draw --------
figure('Color','w', 'Name','PM_Inv Post-Processing', 'Position',[100 100 1000 480]);
axDisp = subplot(1,2,1);
axProf = subplot(1,2,2);
PM_Inv_PlotResults(files, path, opts, axDisp, axProf);
