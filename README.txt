# Source Code for "A Novel Efficient Numerical Method for 3D Bi-Anisotropic Photonic Crystal Simulation"

This repository contains the MATLAB implementation for the paper "A Novel Efficient Numerical Method for 3D Bi-Anisotropic Photonic Crystal Simulation", submitted to Computer Physics Communications (CPC). 

The code is designed to help readers reproduce the algorithms and results presented in the manuscript, including the band structure calculations and light localization simulations.

## Requirements
    Software: MATLAB.
    Dependencies: No external commercial toolboxes are required.

## Directory Structure and Paper Correspondence

The repository is organized to map directly to the sections and figures in the manuscript.

### 1. Main Programs
These scripts are located in the root directory and correspond to the core experiments discussed in the paper:

Section 4.1: Accuracy and efficiency of Algorithm 1:
    band_structure_woodpile_Lanczos_n60_gamma0p2.m: Main script for verifying the accuracy and efficiency of the proposed algorithm on woodpile structures.

Section 4.2: Light localization:
    main_EIG_np_sphere_n60_gamma3p46427.m: Simulates light localization for spherical model.
    main_EIG_np_tri_n60_gamma0p5.m: Simulates light localization for tetrahedron model (gamma = 0.5).
    main_EIG_np_tri_n60_gamma1p009.m: Simulates light localization for tetrahedron model (gamma = 1.009).

### 2. Data Extraction and Visualization ("draw_test/")
The "draw_test" folder contains all the necessary MATLAB scripts to extract the simulated data and plot the figures shown in the manuscript. Each script is named after the corresponding figure.

### 3. Subroutines and Functions
Other folders (such as "fun_Lanczos", "fun_no_inv", "light_locked_Huang", "subroutines", etc.) contain the core functions and solvers called by the main programs. 

## How to Run
1. Open MATLAB and navigate to the root directory of this repository.
2. Ensure that the root directory and all its subfolders are added to the MATLAB path.
3. Run the main scripts to generate the raw data.
4. Navigate to the "draw_test/" folder and run the corresponding plotting scripts to visualize the results.