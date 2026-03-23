# Project Guidelines

## Scope
This workspace combines manuscript sources in Latex/ with research code, datasets, and figure-generation pipelines in src_files/.
Prefer minimal, focused changes. Do not reorganize datasets, generated figures, or archived material unless the task explicitly requires it.

## Architecture
- Latex/ contains the active manuscript, Elsevier class files, and figure assets for the paper.
- src_files/Codes/DMD/ contains the main Fortran 95 SPDMD pipeline.
- src_files/Codes/POD/ contains the main Fortran 95 POD pipeline.
- src_files/Codes/2D_LPT_modified_Hallers/ contains Python and MATLAB code for 2D Lagrangian particle tracking.
- src_files/database/ contains processed case datasets. Each case directory typically includes input_DMD.txt, input_POD.txt, and stacked CSV snapshot data.
- src_files/database_raw/ contains raw extraction utilities and original-case data.
- src_files/data_preprocess_analysis/ contains Python analysis scripts used to inspect preprocessing and reconstruction metrics.
- src_files/Figures_Results/ contains MATLAB plotting scripts for final figures.

## Build And Test
- For manuscript work, compile from Latex/ using pdflatex on On_data_driven.tex unless the task points to a different TeX entrypoint.
- For Fortran work, compile inside the relevant DMD or POD directory with mpifort and the existing HDF5, FFTW3, BLAS, LAPACK, and fortran_stdlib dependencies already implied by the source tree.
- There is no single project-wide build or test runner. Validate changes with the narrowest applicable command or script for the touched area.
- Prefer checking script entrypoints and input files before running expensive data-processing jobs.

## Conventions
- Preserve the current separation between manuscript assets, source code, processed databases, and raw databases.
- Treat src_files/CLAUDE.md as a high-value project reference for scientific context, naming conventions, and workflow details.
- Database folders follow the pattern dataRe<Re>Bo<Bo>_<variant>; result files and input files depend on those names, so avoid renaming folders or files casually.
- Many pipelines assume line-delimited input_DMD.txt or input_POD.txt files with positional parameters. Keep edits format-compatible.
- Scientific outputs are often large HDF5 or CSV artifacts. Avoid regenerating or touching them unless the task requires it.

## Pitfalls
- Case Re40Bo40 is intentionally incomplete in the current study; do not assume every nominal case has generated results.
- DMD preprocessing flags such as mean subtraction, detrending, eigenvalue projection, and derivative handling materially affect stability and reconstruction. Do not change defaults without checking surrounding analysis.
- Reconstruction and mode computations are memory-sensitive and often partition data with np2-style parameters. Prefer small, localized edits over broad refactors in numerical pipelines.
- This repository contains archival and legacy directories. Prefer the active paths under Latex/, src_files/Codes/, src_files/database/, src_files/data_preprocess_analysis/, and src_files/Figures_Results/ unless the task clearly targets archived material.