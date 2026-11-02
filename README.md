# CTRL12-LQR-Fig4-Replication

Team 12 replication project for the EEEN3015J Control Theory module. This repository contains the MATLAB source code for the replication of the LQR tuning result (Fig. 4 and Eqs. 23, 25, 26) from Chirme-Sisa et al. (2024).

## File Structure

- **`CTRL12_replication_figure4.m`**: The main replication script. It reconstructs the linearized state-space model, computes the five LQR candidate gains using the cumulative tuning procedure (Algorithm 1), simulates the closed-loop responses, and generates the penalty values and Fig. 4-style plots.
- **`CTRL12_robustness_settling_band.m`**: The robustness and sensitivity analysis script. It tests how changing the unspecified position settling-time band (1%, 2%, and 5%) affects the selected optimal LQR candidate. It also includes a custom implementation of the continuous LQR solver using the Hamiltonian matrix.

## Requirements & Execution

- **Software**: MATLAB.
- **Toolbox**: 
  - The main script (`CTRL12_replication_figure4.m`) requires the **Control System Toolbox** (for `lqr` and `lsim`).
  - The robustness script (`CTRL12_robustness_settling_band.m`) is written using **base MATLAB** functions and `ode45`, so it runs without additional toolboxes.
- **How to run**: Simply open either `.m` file in MATLAB and click **Run**. The scripts will output the computed gains, penalty tables to the command window, and generate the corresponding response/sensitivity plots.
