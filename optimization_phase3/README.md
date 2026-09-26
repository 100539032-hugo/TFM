
### Phase 3 Optimization: Independent Skin Tailoring (10D)

**Overview**

Phase 3 represents the culmination of the optimization methodology by completely decoupling the upper and lower wingbox covers, expanding the design space to 10 independent continuous variables. This phase evaluates whether assigning different Double-Double (DD) angles and tapering profiles to the tension-dominated lower skin and compression-dominated upper skin yields structural improvements over symmetric designs. To manage the 10D space efficiently, a physics-guided sampling approach is used, allocating the bulk of the material to the critical upper skin to prevent buckling.

**Design Space**

| Parameter | Type | Domain / Value | Description |
| --- | --- | --- | --- |
| **$\phi_{up}$ / $\psi_{up}$** | Variable | $[65^\circ, 90^\circ]$ / $[0^\circ, 20^\circ]$ | Primary/secondary steering angles for the upper skin (compression)|
| **$r_{z, up}$** | Variable (Integer) | Defined by LHS bounds | Upper skin DD block repetitions at Root, Mid, Tip|
| **$\phi_{lo}$ / $\psi_{lo}$** | Variable | $[10^\circ, 40^\circ]$ / $[65^\circ, 90^\circ]$ | Primary/secondary steering angles for the lower skin (tension)|
| **$r_{z, lo}$** | Variable (Integer) | Defined by LHS bounds | Lower skin DD block repetitions at Root, Mid, Tip|

*Constraints:*

1. Independent monotonic physical tapering ($r_{Root, s} \ge r_{Mid, s} \ge r_{Tip, s}$) for both skins $s \in \{up, lo\}$.


2. Minimum combined regional thicknesses ($r_{Root} \ge 40$, $r_{Mid} \ge 30$, $r_{Tip} \ge 22$) and a global aggregate minimum ($r_{\Sigma, 3} \ge 94$).


3. **Physics-guided filter:** The upper skin must contain $\ge 60\%$ of the total aggregate repetitions.



**Target Responses**

| Output | Description |
| --- | --- |
| **$m$ (kg)** | Total structural mass|
| **$u_{max}$ (mm)** | Maximum displacement|
| **$U$ (J)** | Total strain energy (global stiffness)|
| **$\lambda_b$** | First linear buckling load factor|

**File Structure**

* `lhs_part3.m`: Generates the focused 60-sample 10D Latin Hypercube Sampling (LHS) DOE matrix, rigorously enforcing the independent tapering rules and the $\ge 60\%$ upper-skin allocation filter.


* `sim3.hm`: Altair HyperMesh template containing the structural FEA mesh, properties, and spanwise zones decoupled for upper and lower covers.


* `sim3.py`: Python automation script that translates the 10 continuous variables into integer DD stacks, assigns independent PCOMP properties to the top and bottom skins, and creates the `.tcl` execution file.


* `model_part3.m`: MATLAB analysis pipeline that:
* Trains 10D Gaussian Process Regression (GPR) models using Automatic Relevance Determination (ARD) Matérn 5/2 kernels.


* Computes Leave-One-Out Cross-Validation (LOOCV) metrics ($R^2$, RMSE, MAPE) to verify surrogate accuracy.


* Executes a Monte Carlo screening to map the Pareto Front (Mass vs. Strain Energy vs. Buckling).


* Deploys a Genetic Algorithm (GA) to locate the absolute minimum mass optimum in the 10D space subject to $\lambda_b \ge 1.0$.





**Workflow Execution**

1. **Design Sampling:** Run `lhs_part3.m` to generate the 60 geometrically admissible and physics-filtered 10D design points.


2. **FEA Automation:** Execute `sim3.py` to batch-process and create the `.tcl` scripts to execute on `sim3.hm` in OptiStruct.
3. **Surrogate Fitting:** Run `extract_data.m` to parse the results. Run `model_part3.m` to build the highly accurate 10D GPR metamodels and perform LOOCV diagnostics.


4. **Optimization & Verification:** Use `model_part3.m` to run the Genetic Algorithm (GA) and Monte Carlo screening to extract the asymmetric minimum mass and "Knee Point" candidates, followed by final high-fidelity OptiStruct verification.
