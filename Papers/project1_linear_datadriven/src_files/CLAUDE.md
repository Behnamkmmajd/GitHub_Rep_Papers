# Project Guidelines — Linear Data-Driven Methods for Multiphase Rising Bubble Flows

## Overview

PhD research project (NTNU, Dept. of Chemical Engineering) investigating the efficiency of **POD** (Proper Orthogonal Decomposition) and **SPDMD** (Sparsity-Promoting Dynamic Mode Decomposition) as reduced-order models for state estimation of buoyancy-driven rising bubbles. The study compares three reference frames—**Eulerian**, **Lagrangian**, and **Moving-Eulerian**—across five bubble cases in 2D axisymmetric and 3D Cartesian configurations. Key rise characteristics evaluated: **bubble interface**, **rise velocity**, **circularity/sphericity**, and **trajectory**.

Primary finding: **ME-POD** (Moving-Eulerian POD) offers the best balance of accuracy and cost; Lagrangian yields lowest reconstruction error but suffers interpolation artifacts; Eulerian is least efficient.

## Repository Structure

```
Codes/
  DMD/           # Fortran 95 SPDMD pipeline (dmd → amp → reconstruction)
  POD/           # Fortran 95 POD pipeline (pod → reconstruction)
  2D_LPT_modified_Hallers/  # Python/MATLAB Lagrangian particle tracking (2D)
database/        # CSV snapshot datasets per case/framework/observable
data_preprocess/ # Python diagnostic scripts for DMD preprocessing analysis
Figures_Results/ # MATLAB plotting scripts for paper figures (2D/ and 3D/)
LateX_Paper/     # Elsevier cas-sc manuscript (On Data_Driven Methods.tex)
```

## Build and Compile

No Makefile exists. Build via VS Code tasks or manually:

```bash
# DMD programs (from Codes/DMD/)
mpifort -g -I/usr/include -I/usr/include/hdf5/openmpi \
  -I/usr/local/include/fortran_stdlib/GNU-11.4.0 \
  -J./Module ./Module/Sparity.f95 <source>.f95 \
  -o <binary> -llapack -lblas -lfftw3 -lhdf5_fortran -lhdf5 \
  -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi \
  -L/usr/local/lib -lfortran_stdlib

# POD programs (from Codes/POD/) — same but use FFTW_sub_module.f95
mpifort -g ... -J./Module ./Module/FFTW_sub_module.f95 <source>.f95 ...
```

**Dependencies:** GFortran/mpifort (GNU 11.4), LAPACK/BLAS, FFTW3, HDF5 (OpenMPI), fortran_stdlib, OpenMP (POD only)

## Pipeline Execution Order

### DMD
1. `./dmd` — reads `input_DMD.txt` + CSV snapshots → produces `Results_DMD_*.h5` (SVD, eigenvalues, modes, SPDMD matrices P/q/G)
2. `./amp` — reads same HDF5, sweeps γ values → appends sparse amplitudes B_p to HDF5
3. `./reconstruction_opt_error` — computes ‖G − X_DMD‖² error per γ → appends to HDF5
4. `./reconstruction_opt_memory` — full spatial reconstruction (partitioned), re-adds preprocessing → appends xDMD to HDF5

### POD
1. `./pod` — reads `input_POD.txt` + CSV snapshots → produces `Results_POD_*.h5` (eigenvalues, Ψ, Φ, FFT)
2. `./reconstruction` — truncated reconstruction at specified rank → appends xPOD to HDF5

### Post-processing
- `Codes/*/Results/post/*.py` — Python scripts to reshape partitioned HDF5 data into full snapshots
- `Figures_Results/2D/*.m` and `Figures_Results/3D/*.m` — MATLAB scripts generating paper figures

## Database Naming Convention

Folders: `dataRe<Re>Bo<Bo>_<FrameworkObservable>`

| Code | Framework | Observable |
|------|-----------|------------|
| `EG` | Eulerian | Gas volume fraction |
| `EUV` | Eulerian | U,V velocity |
| `LEG_interp` | Moving-Eulerian | Gas fraction (interpolated from DNS onto analytically-tracked grid) |
| `LEUV_interp` / `LEUVW_interp` | Moving-Eulerian | Velocity (2D/3D, interpolated from DNS) |
| `LGHa_Vlinear_Glinear` | Lagrangian (Hallermeier particles) | Gas fraction |
| `LPHa_Vlinear_Glinear` | Lagrangian | Positions (x,y) |
| `LUVHa_Vlinear_Glinear` | Lagrangian | Velocity |

Each folder contains `input_DMD.txt`, `input_POD.txt`, and a `data/` directory with `Stacked_*.csv` snapshot files.

## Physics Cases

| Case | Re | Bo | Coord. System | Solver | Snapshots | Behaviour | Status |
|------|----|----|---------------|--------|-----------|-----------|--------|
| C₁ | 20 | 20 | 2D axisym | SI-LSM_col | 400 | Small deformation, rectilinear | Active |
| C₂ | 40 | 40 | 2D axisym | SI-LSM_col | 400 | Moderate deformation | **Skipped** (no results generated) |
| C₃ | 100 | 100 | 2D axisym | SI-LSM_col | 400 | High deformation, rectilinear | Active |
| C₄ | 168.5 | 2.2 | 3D Cartesian | CDI | 490 | Low deformation, zigzag trajectory | Active |
| C₅ | 168.5 | 12 | 3D Cartesian | CDI | 600 | High deformation, zigzag trajectory | Active |

## Key Algorithms

- **POD:** Method of snapshots (Sirovich 1987). Covariance R = DᵀD → eigendecomposition → spatial modes Φ = DΨΛ⁻¹/². FFT of temporal coefficients via FFTW3.
- **DMD:** Tu et al. Algorithm-1 → reduced dynamics Ã = Uᵣᵀ D″ Wᵣ Σᵣ⁻¹ → eigendecomposition of Ã.
- **SPDMD amplitudes:** ADMM optimization (Jovanović 2014) with warm-starting, adaptive ρ, polishing step.
- **Preprocessing options** (flags in input files): mean subtraction, linear detrending, FD position→velocity conversion, eigenvalue projection onto unit circle.
- **Memory management:** Large spatial fields partitioned row-wise (`np2` blocks) to avoid RAM overflow.

## Result File Naming

`Results_<Method>_<Observable>_<dim>new_<start>_<end>_DR<dr>_VR<vr>_Re_<Re>_Bo_<Bo>_<variant>.h5`

Variants: `mean0`/`mean1` (mean subtraction off/on), `der1` (FD derivative preprocessing)

## Input File Format

### input_DMD.txt (key parameters)
Line-delimited: filename, nt, dt, Time_phase, ns, nv, np, np2, r_mode, r_answer, nx/ny/nz, gamma_n, L_start/L_stop, rho/rho_min/rho_max, mem_man, subtract_mean, detrend_linear, project_eig, derivative

### input_POD.txt (key parameters)
Line-delimited: filename, nt, Time_phase, ns, nv, np, np2, nms, normal, fs, nx/ny/nz, scaling_p, subtract_mean, detrend_linear, derivative

## Code Conventions

- **Language:** Fortran 95 (free-form `.f95`), implicit none throughout
- **I/O:** All results in HDF5; input snapshots as CSV (`Stacked_*.csv`)
- **LAPACK calls:** `dsyev` (symmetric eigenvalue), `dgesvd` (SVD), `zgeev` (complex eigenvalue), `zpotrf/zpotrs` (Cholesky), `zgesv` (complex linear solve)
- **Partitioned spatial operations:** Both POD and DMD loop over `np2` spatial partitions for mode computation/reconstruction to manage memory
- **Plotting:** MATLAB scripts read HDF5 results + `.mat` grid files; figures follow journal style (Elsevier cas-sc)
- **LPT:** Python (tracking_2d_forward_data_driven.py) with Jupyter notebooks for integration/interpolation subfunctions

## Important Notes

- The paper may not reflect the latest work — code and results are more up-to-date; DMD code has been modified and some results regenerated
- **Moving-Eulerian update:** The grid positions are now found analytically (the grid rises with the bubble centroid), and fields (gas fraction, velocity) are **linearly interpolated** from the DNS Eulerian data onto this moving grid. Previously, the original Eulerian grid was shifted and DNS values used directly without interpolation.
- **Evaluation metrics update:** Instead of relative reconstruction error, the current metrics are **cumulative energy** $E_{\text{cum}}(r)$, **r₉₉ rank** (modes for 99% energy), and **normalised spectral entropy** $H_{\text{norm}}$.
- Case C₂ (Re=40, Bo=40) has not been processed — no database or results exist for it
- Validation dataset (`data_DMD_validation/`) uses synthetic two-mode signal from Kutz et al. (2016) SIAM textbook
- Additional test datasets exist for 2D cylinder flow and impinging jet
- Density and viscosity ratios are fixed at 10 for all cases; Fr=1, We=Bo

## Mean Subtraction Preprocessing for SPDMD

### Why it matters

DMD eigenvalues can have $|\mu| > 1$ (outside the unit circle), especially for Moving-Eulerian and Lagrangian data. The Vandermonde matrix $C$ has entries $\mu_j^{k}$, so $|\mu|>1$ causes exponential growth $|\mu|^{m-1}$ that overflows `float64` and corrupts the SPDMD matrices $P$, $q$ (used by ADMM). Mean subtraction (`mean1`) projects eigenvalues onto the unit circle, eliminating this overflow.

### Key findings (from `data_preprocess/diagnostic_eigenvalues.py`)

| Framework | Observable | mean0 Vandermonde | mean1 Vandermonde | Recommendation |
|-----------|-----------|-------------------|-------------------|----------------|
| Eulerian | Gas fraction (EG) | Stable ($|\mu|<1$) | Stable ($|\mu|=1$) | mean1 required for reconstruction quality |
| Eulerian | Velocity (EUV) | Marginally stable | Stable | mean1 preferred |
| Moving-Eulerian | Gas (LEGi) | **OVERFLOW** (log₁₀V up to 34) | Stable | mean1 **required** |
| Moving-Eulerian | Velocity (LEUVi/LEUVWi) | OVERFLOW (log₁₀V up to 25) | Stable | mean1 **required** |
| Lagrangian | Gas (LG) | Marginally stable | Stable | mean1 preferred |
| Lagrangian | Velocity (LUV) | Marginally stable | Stable | mean1 preferred |
| Lagrangian | Position (LP) | **P CORRUPTED** ($|\mu|$ up to 20) | Marginally stable | FD derivative (`der1`) required |

### Fair error comparison

The standard SPDMD error metric $\epsilon = \|G - X_{\text{DMD}}\|_F^2 / \|\Sigma\|_F^2 \times 100$ uses a **different denominator** for mean0 vs mean1 (mean0 includes mean energy, mean1 only fluctuation energy). This makes mean0 appear to have lower error for non-Eulerian frameworks — a **denominator artifact**.

To compare fairly, scale mean1 errors by $\alpha = \|\Sigma_{m1}\|_F^2 / \|\Sigma_{m0}\|_F^2$ (fluctuation-to-total energy ratio):

- **Eulerian Gas:** mean1 genuinely much better (α ≈ 0.82–0.92; error drops from 54–100% to 10–31%)
- **Moving-Eulerian / Lagrangian:** After fair correction, mean0 and mean1 give **comparable** reconstruction quality (within 2–3×). The apparent 10–300× advantage of mean0 was an artifact.
- **At matched ranks** (3D LEGi): mean1 fair error is actually **lower** than mean0 (Re168Bo2: 0.39% vs 1.30%; Re168Bo12: 0.43% vs 1.16%)

**Alpha is exact (not an approximation):** In all 39 HDF5 files, $r = m-1$ (no SVD truncation), so $\|\Sigma_r\|_F^2 = \|\mathcal{D}\|_F^2$ captures 100% of the energy. Alpha is the exact fluctuation-to-total energy ratio.

**Conclusion:** Mean subtraction costs nothing in reconstruction accuracy but is essential for Vandermonde stability and ADMM convergence. Always use `mean1` for SPDMD.

### Mathematical consistency (manuscript ↔ code)

The manuscript describes SPDMD in full space: $\min_b \|\mathcal{D}' - \varphi B_b V_\gamma\|_F^2$. The code works in reduced space: $\min_b \|G - W \text{diag}(b) C\|_F^2$. These are **exactly equivalent** because $\varphi = U_r W$ and $U_r$ has orthonormal columns, so $\|U_r X\|_F = \|X\|_F$. The error formula $\epsilon = \|G - \hat{G}\|_F^2 / \|\Sigma_r\|_F^2 \times 100\%$ equals the full-space relative error $\|\mathcal{D}' - \hat{\mathcal{D}}'\|_F^2 / \|\mathcal{D}'\|_F^2 \times 100\%$.

Full derivation and variable mapping in `data_preprocess/fair_reconstruction_error.tex` (compiled PDF, 8 pages).

**Notation issue:** Manuscript uses $\gamma_j$ for DMD eigenvalues; code and SPDMD literature use $\mu_j$. This clashes with the SPDMD sparsity penalty $\gamma$. Should change to $\mu_j$ in the manuscript.

### ADMM convergence issues (mean0)

- ME gas/velocity 3D cases (Re168Bo2, Re168Bo12): amp fails to converge with mean0 — ADMM oscillates between nnz values, never satisfying convergence criteria even after 150k+ iterations
- LP (Lagrangian Position): P matrix contains NaN/Inf due to catastrophic Vandermonde overflow ($|\mu|^{m-1} > 10^{269}$), making amp impossible

## Diagnostic Tools

### data_preprocess/diagnostic_eigenvalues.py

Scans all `Codes/DMD/Results/Results_DMD_*.h5` files and produces:
- **diagnostics_summary.csv** — full table of eigenvalue stability, reconstruction errors, fair errors
- **diagnostics_summary.xlsx** — Excel version with formatting (two sheets: mean0-only and all variants)
- **eigenvalue_spectra.png** — grid of complex-plane plots per dataset/variant
- **eigenvalue_comparison.png** — paired mean0 vs mean1 overlay plots

Key columns: `|mu|_max`, `log10(Vand_max)`, `Vand_overflow`, `P_has_NaN/Inf`, `recons_err%`, `recons_err%_fair`, `alpha`

Mean0 sheet extra columns (looked up from paired mean1 results):
- `mean1_err%` — raw mean1 reconstruction error at rank ~10
- `mean1_err%_fair` — alpha-corrected fair error at rank ~10
- `mean1_nzval` — rank used for mean1 error
- `mean1_err%@m0rank` — mean1 error at closest rank to mean0's rank (from HDF5 Recons_vec)
- `mean1_err%_fair@m0rank` — alpha-corrected version of above
- `mean1_nzval@m0rank` — actual rank found in mean1 Recons_vec

### data_preprocess/fair_reconstruction_error.tex

LaTeX document (8 pages) containing:
- SPDMD reconstruction derivation matched to code implementation
- Fair reconstruction error formula derivation (α scaling)
- Proof that α is exact (r = m-1 in all files, no SVD truncation)
- Full-space ↔ reduced-space equivalence proof via U_r orthonormality
- SPDMD P/q quadratic expansion matched to Fortran code
- Variable mapping table (manuscript notation ↔ code variables)
- Manuscript–code consistency verdict table

### HDF5 Result Structure (DMD)

```
/Raw/
  EigenV_Real, EigenV_Imag         # Projected eigenvalues (on unit circle if project_eig=1)
  EigenV_Orig_Real, EigenV_Orig_Imag  # Original eigenvalues (pre-projection)
  subtract_mean_flag, detrend_flag, project_eig_flag, fd_preprocessing_flag
  dt
/Standard_Decomposition/
  SVD_S                            # Singular values (r×r diagonal)
  SVD_VT                           # Right singular vectors
  Vr_Real, Vr_Imag                 # DMD eigenvectors
  C_Real, C_Imag                   # Vandermonde matrix (m-1 × r)
  G                                # Reduced data matrix (r × m-1)
  P_Real, P_Imag                   # P = C^H diag(|b_opt|²) C (SPDMD)
  q_Real, q_Imag                   # q = diag(b_opt*) C g  (SPDMD)
  A                                # Ã reduced dynamics matrix
  B_p_Real, B_p_Imag               # Sparse amplitudes per γ (r × n_gamma)
  nzval                            # Number of nonzero modes per γ (1D, length n_gamma)
  Recons_vec                       # Reconstruction results (n_evaluated × 4):
                                   #   col 0: gamma_index, col 1: nzval (rank),
                                   #   col 2: gamma_value, col 3: error (%)
```
