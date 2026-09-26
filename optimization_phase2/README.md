

### Phase 2 Optimization: Five Design Variables and Tapering

**Overview**

Phase 2 expands the design space to five dimensions by introducing regional thickness tapering alongside the Double-Double orientation angles (Φ,Ψ). The upper and lower skins remain symmetric (sharing the same angles and thickness profiles), but the total structural mass is no longer fixed. This phase evaluates the multi-objective trade-offs between mass reduction, global stiffness, and buckling stability, while enforcing a strict monotonic spanwise tapering constraint.

**Design Space**

| Parameter | Type | Domain / Value | Description |
| --- | --- | --- | --- |
| **Φ** | Variable | $[0^\circ, 90^\circ]$ | Common primary steering angle

 |
| **Ψ** | Variable | $[0^\circ, 90^\circ]$ | Common secondary steering angle

 |
| **$r_{Root}$** | Variable (Integer) | Defined by LHS bounds | Number of DD block repetitions at the root

 |
| **$r_{Mid}$** | Variable (Integer) | Defined by LHS bounds | Number of DD block repetitions at midspan

 |
| **$r_{Tip}$** | Variable (Integer) | Defined by LHS bounds | Number of DD block repetitions at the tip

 |

*Constraint:* Monotonic physical tapering $r_{Root} \ge r_{Mid} \ge r_{Tip}$ must be satisfied.

**Target Responses**

| Output | Description |
| --- | --- |
| **$m$ (kg)** | Total structural mass

 |
| **$u_{max}$ (mm)** | Maximum displacement

 |
| **$U$ (J)** | Total strain energy (global stiffness)

 |
| **$\lambda_b$** | First linear buckling load factor

 |

**File Structure**

* `lhs_part2.m`: Generates the initial 40-sample 5D Latin Hypercube Sampling (LHS) DOE matrix enforcing the monotonic tapering constraints.


* `sim2.hm`: Altair HyperMesh template containing the structural FEA mesh, properties, and spanwise zones (Root, Mid, Tip).


* `sim2.py`: Python automation script that translates continuous repeat counts to integer DD stacks, updates FEA properties, and creates the `.tcl` execution file.


* `model_part2.m`: MATLAB analysis pipeline that:
* Trains Gaussian Process Regression (GPR) models using Automatic Relevance Determination (ARD) Matérn 5/2 kernels with a linear basis function.


* Performs uncertainty-driven adaptive infill, adding samples in high-variance regions near the stability boundary.


* Computes Leave-One-Out Cross-Validation (LOOCV) metrics ($R^2$, RMSE) to assess generalization.


* Executes a 20,000-point Monte Carlo screening to generate a 3D Pareto Front (Mass vs. Strain Energy vs. Buckling) and identifies the optimal multiobjective "Knee Point".





**Workflow Execution**

1. **Design Sampling:** Run `lhs_part2.m` to generate the geometrically admissible 5D design points.
2. **FEA Automation:** Execute `sim2.py` to batch-process and create the `.tcl` scripts to execute on `sim2.hm` in OptiStruct.
3. **Surrogate Fitting & Infill:** Run `extract_data.m` to obtain a `.txt` with the initial results. Run `model_part2.m` to build initial GPR metamodels, evaluate predictive variance, and define adaptive infill points.
4. **Optimization:** After running the infill FEA, re-run `model_part2.m` to perform cross-validation diagnostics and conduct the Monte Carlo multiobjective screening to locate the Pareto optima.
