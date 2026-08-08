# PM_Inv: Predominant Mode Surface-Wave Inversion

PM_Inv is a MATLAB package for surface-wave inversion based on the **Predominant Mode (PM)** framework. This method is particularly useful when the target dispersion data is constructed combining multiple active shot, together with ambient-noise passive data. The software estimates shear-wave velocity (Vs) profiles from measured surface-wave dispersion data using a **Hybrid Thin-Layer Method (HTLM)** forward solver coupled with an adaptive **Differential Evolution (DE)** optimization algorithm.

The package implements both:

- **VCPM** – Vertical Component Predominant Mode inversion
- **RCPM** – Radial Component Predominant Mode inversion

---

# Overview

PM_Inv provides two interfaces that execute the **same inversion engine** and produce identical results.

- **GUI:** `PM_Inv_Input_GUI.m` *(recommended)*
- **Script:** `PM_Inv_Input.m`

---

# Requirements

- MATLAB **R2019a** or later
- Parallel Computing Toolbox *(recommended; the code automatically runs in serial mode if unavailable)*
- Statistics and Machine Learning Toolbox

The GUI automatically checks for the required toolboxes during start-up and notifies the user if any dependency is missing.

---

# Repository Structure

```text
src/
    PM_Inv_Input_GUI.m
    PM_Inv_Input.m
    PM_Inv_PostProcessing.m
    helpers/

Examples/
    Example datasets and inversion results

Input/
    Example dispersion files

docs/
    User manual
```

---

# Quick Start

## Graphical User Interface (Recommended)

1. Open MATLAB.
2. Set the repository as the **Current Folder**.
3. Run

```matlab
PM_Inv_Input_GUI
```

4. Complete the GUI tabs from left to right:

- **Target Data** – Select the dispersion file, define the uncertainty type (COV/STD), and choose the measured component (VCPM or RCPM).
- **Layering & Model** – Define the layer parameterization (`h_count`), initial model, and material properties.
- **Constraints & Optimizer** – Configure Differential Evolution parameters, DCR, and geological constraints.
- **Layer Limits** *(optional)* – Compute and manually edit layer-wise search bounds.

5. Click **Preview Target** to verify the imported dispersion curve.
6. Click **Run Inversion** to start the inversion.

During the inversion, PM_Inv continuously displays:

- Measured and predicted dispersion curves
- Inverted Vs profile
- Differential Evolution convergence history

The inversion can be **paused** or **terminated** at any time. Partial results are automatically saved.

---

# Post-Processing

Previously saved inversion results can be visualized without rerunning the inversion.

### GUI

Open the **Post-Processing** tab and

1. Select one or more `.mat` files.
2. Choose the plotting mode.
3. Click **Plot**.

### Script

Run

```matlab
PM_Inv_PostProcessing
```

Available plotting modes include

- **Best**
- **All**
- **MeanBest**

Additional options include

- Dmax
- Nbest
- Number of theoretical modes

---

# Input Dispersion File

The input dispersion file should contain

| Column | Description |
|---------|-------------|
| 1 | Frequency (Hz) |
| 2 | Phase velocity (m/s) or slowness |
| 3 *(optional)* | Measurement uncertainty (STD or COV) |

If the third column contains the **Coefficient of Variation (COV)**, set **Third-column Interpretation = COV** in the GUI.

---

# Notes

- Multiple layer parameterizations can be evaluated simultaneously using MATLAB indexing (e.g., `5:7`).
- The inversion results depend strongly on the selected parameterization. Users are encouraged to evaluate several `h_count` values and independent inversion trials.
- Results are automatically organized into folders named after the target dataset.
- PM_Inv is an actively developed research software package and will continue to evolve with ongoing research. Future releases will include bug fixes, performance improvements, and additional inversion capabilities.

---

# Documentation

A detailed user guide is available in

```
PM_Inv_UserManual.pdf
```

---

# Citation

If you use **PM_Inv** in your research, please cite the accompanying publication.

[1]. Bhaumik, M., & Cox, B. R. (2026). Radial-Component Predominant-Mode Inversion of Rayleigh Waves: Application to DAS-based Site Characterization. arXiv preprint arXiv:2605.16717. https://doi.org/10.48550/arXiv.2605.16717
> *(Paper citation will be updated after publication.)*
> 
[2]. Bhaumik, M., & Cox, B. R. (2026). Predominant-mode inversion of surface waves: Inherently addressing inconspicuous low frequency mode jumps. Engineering Geology, 108834. https://doi.org/10.1016/j.enggeo.2026.108834

[3]. Bhaumik, M., & Naskar, T. (2024). Computation of Surface Wave's Dominating Mode for Stratified Media. Indian Geotechnical Journal, 55(3), 1699-1714. https://doi.org/10.1007/s40098-024-01045-x


A Zenodo DOI for the software will be provided with the first official release.

---

# Contact

**Mrinal Bhaumik, Ph.D.**

Utah State University

📧 mrinal.bhaumik@usu.edu

📧 mrinal.bhaumik2012@gmail.com

Questions, bug reports, feature requests, collaboration opportunities, and feedback are welcome.

---

# Disclaimer

This software is provided **"as is"**, without warranty of any kind, express or implied. While every effort has been made to verify the correctness of the implementation, the authors assume no responsibility for errors or consequences arising from the use of this software.
