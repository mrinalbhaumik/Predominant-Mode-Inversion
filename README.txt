
=====================================================================================================================================
============================= Differential Evolution (DE) for Predominant Mode Inversion (PMinv) ==============================
===================================================================================================================================== 
It has both VCPM (Vertical component predominant mode) and RCPM (Radial component predominant mode) Inversion

************************ OVERVIEW ************************

- Two ways to run, both giving the SAME inversion:
    * GUI    : PM_Inv_Input_GUI.m   (recommended)
    * Script : PM_Inv_Input.m       (edit values at the top, then Run)

************************ REQUIREMENTS ************************

- MATLAB R2019a or newer (uses uifigure / uigridlayout).
- Parallel Computing Toolbox   (DE uses parfor; runs serially if absent).
- Statistics and Machine Learning Toolbox.
- The GUI checks these at start-up and points to installers if any is missing.

************************ KEY FILES ************************

- PM_Inv_Input_GUI.m ........ interactive GUI (recommended).
- PM_Inv_Input.m ............ script version of the same inversion.
- PM_Inv_PostProcessing.m ... plot saved results (also built into the GUI).
- helpers\ .................. supporting functions (auto-added to the path).
- 50m_R0_dispersion_stats.txt  example target dispersion file.

************************ HOW TO RUN THE GUI ************************

- Open MATLAB and set the folder as the Current Folder.
- In the Command Window type:  PM_Inv_Input_GUI (or open PM_Inv_Input_GUI.m and press Run).
- Fill the tabs left to right (Point your cursor on any field for a description):
    * Target Data ............ pick the dispersion file; set "3rd-column meaning" (COV/STD) and component (VCPM/RCPM).
    * Layering & Model ....... set h_count (batch OK, e.g. 5:7) and Vs/nu/rho.
    * Constraints & Optimizer  DE settings, DCR, and Vs constraints.
    * Layer Limits ........... (optional) Compute, then hand-edit per-layer bounds.
- Click "Preview Target" to check the loaded curve.
- Click "Run Inversion" - watch the live dispersion / Vs / convergence panels
  in "Figures & Control"; use Pause / Stop any time.
- Results are saved automatically to a new subfolder named after the target
  file (one .mat per h_count and trial, plus input_params.mat).

************************ VIEW / PLOT RESULTS ************************

- GUI: open the "Post-Processing" tab -> Select .mat -> Plot.
- Or run PM_Inv_PostProcessing.m and select the saved .mat files.
- Options: Dmax (depth), Mode (Best / All / MeanBest), Nbest, number of modes.

************************ QUICK NOTES ************************
- Dispersion file columns: 1 = frequency (Hz), 2 = phase velocity (or slowness),
  3 = optional uncertainty. If column 3 is a COV, set "3rd-column meaning" = COV
  (otherwise the misfit weighting is wrong).
- h_count runs in BATCH: e.g. 5:7 or 5,6,8 launches several models in one go.
- Accuracy depends strongly on the parameterisation - try a few settings/trials.

************************ CONTACT ************************
- Mrinal Bhaumik  -  mrinal.bhaumik2012@gmail.com  /  mrinal.bhaumik@usu.edu
- The authors provide no warranty; use at your own risk.
  For issues, feedback, collaboration opportunities, or any explanation or suggestions, 
  please feel free to reach out via email. We welcome discussion and user input 
  to further improve this package.

************************ IMPORTANT NOTE ON PERFORMANCE:************************
The performance and accuracy of the inversion heavily depend on the model parameterization. 
It is strongly recommended to experiment with different parameterizations and perform 
multiple trial runs to achieve stable and meaningful results.
===============================================================================
