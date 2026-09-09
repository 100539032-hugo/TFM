# Phase 1 Optimization: Angle Isolation ($\Phi, \Psi$)

## Overview
Phase 1 isolates the influence of Discrete Draping (DD) orientation angles ($\Phi, \Psi$) on overall structural behavior. A single uniform angle pair is applied across both the upper and lower covers as well as all three spanwise zones, holding the laminate repetition vector constant at $\mathbf{r} = [26, 21, 15]^T$ (Mass = 420.52 kg).

---

## Design Space

| Parameter | Type | Domain / Value | Description |
| :--- | :--- | :--- | :--- |
| **$\Phi$** | Variable | $[0^\circ, 90^\circ]$ | Primary steering angle |
| **$\Psi$** | Variable | $[0^\circ, 90^\circ]$ | Secondary steering angle |
| **$\mathbf{r}$** | Fixed | $[26, 21, 15]^T$ | Spanwise laminate repetitions (Root, Mid, Tip) |
| **Target Responses** | Output | $u_{\max}$ (mm), $U$ (J), $\lambda_b$ | Max displacement, Strain energy, Buckling factor |

---

## File Structure

* **`lhs_part1.m`**: Generates the 20-sample Latin Hypercube Sampling (LHS) DOE matrix for $(\Phi, \Psi) \in [0^\circ, 90^\circ]^2$.
* **`sim1.hm`**: Altair HyperMesh template containing the structural FEA mesh and laminate definitions.
* **`sim1.py`**: Python automation script that updates FEA properties, and creates the .tcl.
* **`model_part1.m`**: MATLAB analysis pipeline that:
  * Trains Gaussian Process Regression (GPR) models with Matérn 5/2 kernels.
  * Computes Leave-One-Out Cross-Validation (LOOCV) metrics ($R^2$, RMSE).
  * Generates response surfaces, parity plots, and prediction uncertainty maps.
  * Conducts constrained minimization of strain energy $U$ subject to buckling stability $\lambda_b \ge 1.0$.

---

## Workflow Execution

1. **Design Sampling**: Run `lhs_part1.m` to generate the design point coordinates.
2. **FEA Automation**: Execute `sim1.py` to batch-process create the .tcl to execute on `sim1.hm` in OptiStruct.
3. **Surrogate Fitting & Optimization**: Run `extract_data.m` to obtain a .xtx with the results and then run `model_part1.m` to build GPR metamodels, perform cross-validation diagnostics, and identify optimal angle combinations.
