# Master's Thesis: Structural Optimization of Aerostructures using Double-Double Composite Laminates

**Author:** Hugo Blanco Quintana  
**Institution:** Universidad Carlos III de Madrid (UC3M) - Master's in Aeronautical Engineering  
**Advisor:** Alberto Racionero Sánchez-Majano

## 📖 Project Overview

This repository contains the source code, automation scripts, and methodology developed for my Master's Thesis. The project focuses on the structural optimization of a representative composite wingbox by replacing conventional laminates with the novel **Double-Double (DD)** laminate architecture. 

The primary objective is to exploit the continuous design space of DD laminates to minimize the structural mass of the wingbox while guaranteeing the stringent aeronautical requirements for global stiffness and linear buckling resistance.

## 🚀 Key Features

The project integrates an advanced computational framework that combines:
*   **FEA Automation:** A progressive pipeline connecting Python, Tcl (HyperMesh), and Altair OptiStruct to generate, run, and evaluate thousands of designs without manual intervention.
*   **Design of Experiments (DoE):** Latin Hypercube Sampling (LHS) enhanced with physics-guided and geometric filters.
*   **Surrogate Modeling (Machine Learning):** Training Gaussian Process Regression (GPR/Kriging) models in MATLAB, featuring adaptive infill updating driven by predictive uncertainty.
*   **Multi-Objective Optimization:** Pareto Front mapping and absolute minimum mass identification (using Genetic Algorithms and Monte Carlo screening) to evaluate the trade-offs between mass, strain energy, and buckling stability.

## 📂 Optimization Phases

The code is structured to replicate the three incremental phases of the research:
1.  **Phase 1 (Angle Isolation):** Optimization of the continuous steering angles ($\phi, \psi$) at a constant thickness to maximize stiffness and directional stability.
2.  **Phase 2 (Symmetric Tapering):** Introduction of a 5D design space to control spanwise thickness reduction (Root, Mid, Tip) while maintaining strict symmetry between the upper and lower skins.
3.  **Phase 3 (Independent 10D Tailoring):** Complete decoupling of the compression-dominated upper skin and tension-dominated lower skin into 10 independent variables, yielding the most highly tailored and lightweight designs.

This work is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND).
