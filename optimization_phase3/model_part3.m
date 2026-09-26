% =========================================================================
% PHASE 3: SURROGATE MODELING, LOOCV DIAGNOSTICS & OPTIMIZATION (10D)
% =========================================================================
% Description:
%   This script builds Gaussian Process Regression (GPR / Kriging) surrogate
%   models with Automatic Relevance Determination (ARD) Matern 5/2 kernel 
%   for a 10D composite wing box design space. It performs Leave-One-Out 
%   Cross-Validation (LOOCV), multi-objective Pareto optimization via Monte 
%   Carlo sampling, multi-algorithm comparison (fmincon & GA), sensitivity 
%   analysis, and adaptive sampling selection for Phase 4.
% =========================================================================

clear; clc; close all;
rng(42, 'twister'); % Set random seed for reproducibility

%% 1. GLOBAL GRAPHICAL CONFIGURATION
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16); % Base font size for figures
fs = 16;

%% 2. DATA LOADING (10D DESIGN SPACE & FEA RESULTS)
% Define lower and upper bounds for the 10 design variables:
% [Phi_up, Psi_up, r_R_u, r_M_u, r_T_u, Phi_lo, Psi_lo, r_R_l, r_M_l, r_T_l]
lb = [65,  0, 24, 18, 13, 10, 65,  8,  6,  4];
ub = [90, 20, 30, 28, 25, 40, 90, 16, 13, 10];

var_names = {'\Phi_{up}', '\Psi_{up}', 'r_{R,u}', 'r_{M,u}', 'r_{T,u}', ...
             '\Phi_{lo}', '\Psi_{lo}', 'r_{R,l}', 'r_{M,l}', 'r_{T,l}'};

% Import 10D Latin Hypercube Sampling (LHS) design points
design_file = 'lhs_phase3_60cases2.txt';
if exist(design_file, 'file')
    opts_txt = detectImportOptions(design_file);
    opts_txt.VariableNamesLine = 4;
    opts_txt.DataLines = [6, Inf];
    data_X = readtable(design_file, opts_txt);
    X_fea = table2array(data_X(:, 3:12));
else
    error('10D design file (%s) not found.', design_file);
end

% Import finite element analysis (FEA) KPI results
fea_file = 'OptiStruct_KPI_Resultados_Phase3.txt'; 
if exist(fea_file, 'file')
    resTable = readtable(fea_file);
    y_mass   = resTable.Masa_Total;
    y_disp   = resTable.Desplazamiento_Max;
    y_strain = resTable.Energia_Deformacion;
    y_blf    = resTable.Factor_Pandeo_BLF;
else
    warning('FEA result file not found. Generating synthetic test fallback data...');
    n_dummy  = size(X_fea, 1);
    y_mass   = 12.5 + sum(X_fea(:, [3:5, 8:10]), 2)*0.18 + randn(n_dummy,1)*0.05;
    y_disp   = (120 ./ (sum(X_fea(:, [3:5, 8:10]), 2) + 1e-3))*1e-3 + randn(n_dummy,1)*1e-5;
    y_strain = 4500 ./ (sum(X_fea(:, [3:5, 8:10]), 2) + 1e-3) + randn(n_dummy,1)*5;
    y_blf    = 0.5 + 0.02 * (X_fea(:,3)*1.2 + X_fea(:,8)*0.8) .* cosd(X_fea(:,1)-75) + rand(n_dummy,1)*0.1;
end

n_samples = length(y_mass);

%% 3. KRIGING (GPR) SURROGATE MODEL TRAINING WITH ARD
% Hyperparameter configuration: Matern 5/2 kernel with Automatic Relevance Determination
opts_ard = {'KernelFunction', 'ARDMatern52', 'BasisFunction', 'linear', 'Standardize', true};

fprintf('Training GPR Metamodels (10D)...\n');
gpr_mass   = fitrgp(X_fea, y_mass, opts_ard{:});
gpr_disp   = fitrgp(X_fea, y_disp, opts_ard{:});
gpr_strain = fitrgp(X_fea, y_strain, opts_ard{:});
gpr_blf    = fitrgp(X_fea, y_blf, opts_ard{:});

%% 4. LOOCV EVALUATION (LEAVE-ONE-OUT CROSS-VALIDATION)
% Performance metrics definitions
calc_r2   = @(y, p) 1 - sum((y - p).^2) / sum((y - mean(y)).^2);
calc_rmse = @(y, p) sqrt(mean((y - p).^2));
calc_mape = @(y, p) mean(abs((y - p) ./ y)) * 100;

% Predict via Leave-One-Out Cross-Validation
p_mass   = kfoldPredict(crossval(gpr_mass, 'KFold', n_samples));
p_disp   = kfoldPredict(crossval(gpr_disp, 'KFold', n_samples));
p_strain = kfoldPredict(crossval(gpr_strain, 'KFold', n_samples));
p_blf    = kfoldPredict(crossval(gpr_blf, 'KFold', n_samples));

% Compute validation metrics
r2_m = calc_r2(y_mass, p_mass);     rmse_m = calc_rmse(y_mass, p_mass);     mape_m = calc_mape(y_mass, p_mass);
r2_d = calc_r2(y_disp, p_disp);     rmse_d = calc_rmse(y_disp*1e3, p_disp*1e3); mape_d = calc_mape(y_disp, p_disp);
r2_s = calc_r2(y_strain, p_strain); rmse_s = calc_rmse(y_strain, p_strain); mape_s = calc_mape(y_strain, p_strain);
r2_b = calc_r2(y_blf, p_blf);       rmse_b = calc_rmse(y_blf, p_blf);       mape_b = calc_mape(y_blf, p_blf);

fprintf('\n=================================================================================\n');
fprintf('                     SURROGATE MODEL METRICS (LOOCV DIAGNOSTICS)                 \n');
fprintf('=================================================================================\n');
fprintf('Total Mass    -> R2: %.4f | RMSE: %.4f kg | MAPE: %.2f%%\n', r2_m, rmse_m, mape_m);
fprintf('Max Displ.    -> R2: %.4f | RMSE: %.4f mm | MAPE: %.2f%%\n', r2_d, rmse_d, mape_d);
fprintf('Strain Energy -> R2: %.4f | RMSE: %.4f J  | MAPE: %.2f%%\n', r2_s, rmse_s, mape_s);
fprintf('Buckling BLF  -> R2: %.4f | RMSE: %.4f    | MAPE: %.2f%%\n', r2_b, rmse_b, mape_b);
fprintf('=================================================================================\n\n');

%% 5. VISUAL DIAGNOSTICS: LOOCV PARITY & RESIDUAL PLOTS
% Parity Plot: Total Mass
figure('Name', 'LOOCV Parity - Total Mass', 'Color', 'w', 'Position', [100 100 700 550]);
plot(y_mass, p_mass, 'bo', 'MarkerFaceColor', [0.2 0.6 1], 'MarkerSize', 8); hold on;
plot([min(y_mass) max(y_mass)], [min(y_mass) max(y_mass)], 'k--', 'LineWidth', 1.5);
title(sprintf('Total Mass ($R^2 = %.4f$, MAPE = $%.2f\\%%$)', r2_m, mape_m), 'Interpreter', 'latex', 'FontSize', fs + 3); 
xlabel('FEA Mass (kg)', 'Interpreter', 'latex', 'FontSize', fs + 2); 
ylabel('GPR Mass (kg)', 'Interpreter', 'latex', 'FontSize', fs + 2); grid on;

% Parity Plot: Max Displacement
figure('Name', 'LOOCV Parity - Max Displacement', 'Color', 'w', 'Position', [120 120 700 550]);
plot(y_disp*1e3, p_disp*1e3, 'sq', 'MarkerFaceColor', [0.9 0.4 0.2], 'MarkerSize', 8); hold on;
plot([min(y_disp) max(y_disp)]*1e3, [min(y_disp) max(y_disp)]*1e3, 'k--', 'LineWidth', 1.5);
title(sprintf('Max Displacement ($R^2 = %.4f$, MAPE = $%.2f\\%%$)', r2_d, mape_d), 'Interpreter', 'latex', 'FontSize', fs + 3); 
xlabel('FEA Displacement (mm)', 'Interpreter', 'latex', 'FontSize', fs + 2); 
ylabel('GPR Predicted Displacement (mm)', 'Interpreter', 'latex', 'FontSize', fs + 2); grid on;

% Parity Plot: Strain Energy
figure('Name', 'LOOCV Parity - Strain Energy', 'Color', 'w', 'Position', [140 140 700 550]);
plot(y_strain, p_strain, 'd', 'MarkerFaceColor', [0.4 0.8 0.3], 'MarkerSize', 8); hold on;
plot([min(y_strain) max(y_strain)], [min(y_strain) max(y_strain)], 'k--', 'LineWidth', 1.5);
title(sprintf('Strain Energy ($R^2 = %.4f$, MAPE = $%.2f\\%%$)', r2_s, mape_s), 'Interpreter', 'latex', 'FontSize', fs + 3); 
xlabel('FEA Strain Energy $U$ (J)', 'Interpreter', 'latex', 'FontSize', fs + 2); 
ylabel('GPR Strain Energy $U$ (J)', 'Interpreter', 'latex', 'FontSize', fs + 2); grid on;

% Parity Plot: Buckling Load Factor
figure('Name', 'LOOCV Parity - Buckling Factor', 'Color', 'w', 'Position', [160 160 700 550]);
plot(y_blf, p_blf, '^', 'MarkerFaceColor', [0.8 0.2 0.8], 'MarkerSize', 8); hold on;
plot([min(y_blf) max(y_blf)], [min(y_blf) max(y_blf)], 'k--', 'LineWidth', 1.5);
title(sprintf('Buckling Factor $\\lambda_b$ ($R^2 = %.4f$, MAPE = $%.2f\\%%$)', r2_b, mape_b), 'Interpreter', 'latex', 'FontSize', fs + 3); 
xlabel('FEA Buckling Load Factor $\lambda_b$', 'Interpreter', 'latex', 'FontSize', fs + 2); 
ylabel('GPR Buckling Load Factor $\lambda_b$', 'Interpreter', 'latex', 'FontSize', fs + 2); grid on;

% Residual Plot: Buckling Load Factor
figure('Name', 'Residual Analysis - Buckling Factor', 'Color', 'w', 'Position', [180 180 700 550]);
res_b = y_blf - p_blf;
stem(p_blf, res_b, 'filled', 'Color', [0.5 0 0.5], 'LineWidth', 1.2); hold on;
yline(0, 'k--', 'LineWidth', 1.5);
title('Buckling: LOOCV Residual Distribution', 'Interpreter', 'latex', 'FontSize', fs + 3);
xlabel('GPR Predicted $\lambda_b$', 'Interpreter', 'latex', 'FontSize', fs + 2); 
ylabel('Residual ($y - \hat{y}$)', 'Interpreter', 'latex', 'FontSize', fs + 2); grid on;

%% 6. MONTE CARLO EXPLORATION & 10D CONSTRAINTS FILTERING
N_mc = 60000;
X_mc_raw = lb + rand(N_mc, 10) .* (ub - lb);

% Discretize ply repetition variables (columns 3-5 and 8-10)
X_mc_disc = X_mc_raw;
X_mc_disc(:, [3:5, 8:10]) = round(X_mc_raw(:, [3:5, 8:10]));

% Physical structural constraints:
% 1. Monotonicity across span (Root >= Mid >= Tip)
mon_up = (X_mc_disc(:,3) >= X_mc_disc(:,4)) & (X_mc_disc(:,4) >= X_mc_disc(:,5));
mon_lo = (X_mc_disc(:,8) >= X_mc_disc(:,9)) & (X_mc_disc(:,9) >= X_mc_disc(:,10));

% 2. Minimum total section repetitions
sum_root = X_mc_disc(:, 3) + X_mc_disc(:, 8); 
sum_mid  = X_mc_disc(:, 4) + X_mc_disc(:, 9); 
sum_tip  = X_mc_disc(:, 5) + X_mc_disc(:, 10);
local_mask = (sum_root >= 40) & (sum_mid >= 30) & (sum_tip >= 22);

% 3. Aggregate structural stiffness constraint
r_sig_3 = sum_root + sum_mid + sum_tip;
aggregate_mask = (r_sig_3 >= 94);

% 4. Load allocation ratio constraint (Upper skin proportion >= 60%)
allocation_ratio = sum(X_mc_disc(:, 3:5), 2) ./ r_sig_3;
allocation_mask = (allocation_ratio >= 0.60);

% Combine all structural feasibility masks
valid_mc_idx = mon_up & mon_lo & local_mask & aggregate_mask & allocation_mask;
X_mc_valid = X_mc_disc(valid_mc_idx, :);

% Surrogate evaluation on valid MC candidates
mc_mass   = predict(gpr_mass, X_mc_valid);
mc_strain = predict(gpr_strain, X_mc_valid);
mc_blf    = predict(gpr_blf, X_mc_valid);

% Filter stable candidates (Buckling Load Factor >= 1.0)
stable_idx  = mc_blf >= 1.0;
stable_X    = X_mc_valid(stable_idx, :);
stable_mass = mc_mass(stable_idx);
stable_blf  = mc_blf(stable_idx);
stable_U    = mc_strain(stable_idx);

% Identify single-objective optima within feasible set
[min_m_val, idx_min_m] = min(stable_mass);
[min_u_val, idx_min_u] = min(stable_U);
[max_b_val, idx_max_b] = max(stable_blf);

%% 7. 3D OBJECTIVE SPACE & PARETO FRONT IDENTIFICATION
P_mass   = stable_mass(:);
P_strain = stable_U(:);
P_blf    = stable_blf(:);

% Formulate minimization problem: [Mass, Strain Energy, -Buckling]
objs = [P_mass, P_strain, -P_blf];
N_pts = size(objs, 1);

% Non-dominated sorting for 3D Pareto Front
is_pareto = true(N_pts, 1);
for i = 1:N_pts
    dominates = all(objs <= objs(i,:), 2) & any(objs < objs(i,:), 2);
    if any(dominates)
        is_pareto(i) = false;
    end
end

pareto_mass   = P_mass(is_pareto);
pareto_strain = P_strain(is_pareto);
pareto_blf    = P_blf(is_pareto);

% Min-Max Normalization to identify Best Compromise (Knee Point)
m_n = (pareto_mass - min(pareto_mass)) / (max(pareto_mass) - min(pareto_mass) + 1e-9);
u_n = (pareto_strain - min(pareto_strain)) / (max(pareto_strain) - min(pareto_strain) + 1e-9);
l_n = (max(pareto_blf) - pareto_blf) / (max(pareto_blf) - min(pareto_blf) + 1e-9);

dist_ideal = sqrt(m_n.^2 + u_n.^2 + l_n.^2);
[~, idx_knee_pareto] = min(dist_ideal);

indices_pareto  = find(is_pareto);
idx_knee_stable = indices_pareto(idx_knee_pareto);

% Focused axis limits for clean rendering
x_min_focus = 340;  x_max_focus = 420;
y_min_focus = 2200; y_max_focus = 4100;
z_min_focus = 0.3;  z_max_focus = 2.0;

figure('Name', 'Fig 4C: 3D Pareto Objective Space', 'Color', 'w', 'Position', [150, 150, 850, 650]);
hold on; box on; grid on;

% 1. Stability Limit Plane (lambda_b = 1.0)
[X_grid, Y_grid] = meshgrid(linspace(x_min_focus, x_max_focus, 15), ...
                            linspace(y_min_focus, y_max_focus, 15));
h_lim = surf(X_grid, Y_grid, ones(size(X_grid)), 'FaceColor', 'r', ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none');

% 2. Monte Carlo Candidate Point Scatter
h_inf = scatter3(mc_mass(~stable_idx), mc_strain(~stable_idx), mc_blf(~stable_idx), ...
                 10, [0.85 0.85 0.85], 'filled', 'MarkerEdgeAlpha', 0.15);
h_stab = scatter3(P_mass(~is_pareto), P_strain(~is_pareto), P_blf(~is_pareto), ...
                  16, [0.4 0.6 0.9], 'filled', 'MarkerFaceAlpha', 0.25);
h_par = scatter3(pareto_mass, pareto_strain, pareto_blf, ...
                 50, [0.9 0.1 0.1], 'filled', 'MarkerEdgeColor', 'k');
h_fea = scatter3(y_mass, y_strain, y_blf, ...
                 65, 'y', '^', 'filled', 'MarkerEdgeColor', 'k');

% 3. Key Optimal Designs Highlighted
h_opt_m = plot3(stable_mass(idx_min_m), stable_U(idx_min_m), stable_blf(idx_min_m), ...
                'p', 'MarkerSize', 18, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k'); 
h_opt_u = plot3(stable_mass(idx_min_u), stable_U(idx_min_u), stable_blf(idx_min_u), ...
                's', 'MarkerSize', 15, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k'); 
h_opt_b = plot3(stable_mass(idx_max_b), stable_U(idx_max_b), stable_blf(idx_max_b), ...
                'd', 'MarkerSize', 15, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k'); 
h_knee  = plot3(stable_mass(idx_knee_stable), stable_U(idx_knee_stable), stable_blf(idx_knee_stable), ...
                'h', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'k'); 

% Axes Formatting
xlabel('Total Mass (kg)', 'Interpreter', 'latex', 'FontSize', fs + 1);
ylabel('Strain Energy $U$ (J)', 'Interpreter', 'latex', 'FontSize', fs + 1);
zlabel('Buckling Load Factor $\lambda_b$', 'Interpreter', 'latex', 'FontSize', fs + 1);
title('3D Objective Trade-off \& Pareto Front Identification', 'Interpreter', 'latex', 'FontSize', fs + 3);

xlim([x_min_focus, x_max_focus]);
ylim([y_min_focus, y_max_focus]);
zlim([z_min_focus, z_max_focus]);
view(-37.5, 22);

legend([h_inf, h_stab, h_par, h_knee, h_opt_m, h_opt_u, h_opt_b, h_fea, h_lim], ...
       {'Infeasible ($\lambda_b < 1.0$)', ...
        'Feasible Domain', ...
        'Pareto Front', ...
        'Best Compromise (Knee)', ...
        'Min Mass Design', ...
        'Min Strain Energy', ...
        'Max $\lambda_b$ Design', ...
        'FEA Samples', ...
        'Stability Limit $\lambda_b = 1.0$'}, ...
       'Location', 'northeast', ...
       'Interpreter', 'latex', ...
       'FontSize', fs - 4);

set(gca, 'FontSize', fs - 1, 'TickLabelInterpreter', 'latex');
hold off;

%% EXPORT KRIGING KPI PREDICTIONS FOR OPTIMAL CASES & KNEE POINT
opt_indices = [idx_min_m, idx_min_u, idx_max_b, idx_knee_stable];
opt_types   = {'MIN_MASS', 'MIN_STRAIN', 'MAX_BUCKL', 'KNEE_POINT'};

opt_disp_mm = predict(gpr_disp, stable_X(opt_indices, :)) * 1e3;

fprintf('\n========================================================================================================\n');
fprintf('                         KRIGING MODEL PREDICTIONS (OPTIMAL KPIS)                                       \n');
fprintf('========================================================================================================\n');
fprintf('%-6s %-12s %-15s %-18s %-18s %-12s\n', ...
    'Pt', 'Type', 'Mass (kg)', 'Strain Energy (J)', 'Max Disp (mm)', 'Buckling (BLF)');
fprintf('--------------------------------------------------------------------------------------------------------\n');

for k = 1:length(opt_indices)
    idx = opt_indices(k);
    fprintf('%-6d %-12s %-15.3f %-18.2f %-18.3f %-12.4f\n', ...
        k, opt_types{k}, stable_mass(idx), stable_U(idx), opt_disp_mm(k), stable_blf(idx));
end
fprintf('========================================================================================================\n\n');

%% 8. MULTI-ALGORITHM OPTIMIZATION & UNCERTAINTY BENCHMARK
obj_fun = @(x) predict(gpr_mass, x);

% Linear inequality constraints for monotonicity: A*x <= b_lin
A = zeros(4, 10);
A(1, 3) = -1; A(1, 4) = 1;  % r_M_u - r_R_u <= 0
A(2, 4) = -1; A(2, 5) = 1;  % r_T_u - r_M_u <= 0
A(3, 8) = -1; A(3, 9) = 1;  % r_M_l - r_R_l <= 0
A(4, 9) = -1; A(4, 10) = 1; % r_T_l - r_M_l <= 0
b_lin = zeros(4, 1);

% Non-linear constraints function definition
nonlcon_10d = @(x) deal([1.0 - predict(gpr_blf, x); ...
                         40 - (x(3)+x(8)); ...            
                         30 - (x(4)+x(9)); ...            
                         22 - (x(5)+x(10)); ...           
                         94 - sum(x([3:5, 8:10])); ...    
                         0.60 - (sum(x(3:5))/sum(x([3:5, 8:10])))], []);

% 1. Monte Carlo Minimum
tic; x_opt_mc = stable_X(idx_min_m, :); mass_opt_mc = min_m_val; t_mc = toc;
blf_opt_mc = stable_blf(idx_min_m);

% 2. Gradient-based SQP Optimization (fmincon)
opts_fmc = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'TolFun', 1e-5);
tic; [x_opt_fmc, mass_opt_fmc] = fmincon(obj_fun, x_opt_mc, A, b_lin, [], [], lb, ub, nonlcon_10d, opts_fmc); t_fmc = toc;
blf_opt_fmc = predict(gpr_blf, x_opt_fmc);

% 3. Genetic Algorithm Optimization (GA)
opts_ga = optimoptions('ga', 'PopulationSize', 100, 'MaxGenerations', 60, 'Display', 'off');
tic; [x_opt_ga, mass_opt_ga] = ga(obj_fun, 10, A, b_lin, [], [], lb, ub, nonlcon_10d, opts_ga); t_ga = toc;
blf_opt_ga = predict(gpr_blf, x_opt_ga);

% Uncertainty Quantification (GPR Standard Deviation)
X_opts = [x_opt_mc; x_opt_fmc; x_opt_ga];
[~, sig_m] = predict(gpr_mass, X_opts);
[b_pred, sig_b] = predict(gpr_blf, X_opts);
ic95_b_low  = b_pred - 1.96*sig_b;
ic95_b_high = b_pred + 1.96*sig_b;

T_alg = table({'Monte Carlo'; 'fmincon (SQP)'; 'GA'}, ...
    [mass_opt_mc; mass_opt_fmc; mass_opt_ga], sig_m, ...
    [blf_opt_mc; blf_opt_fmc; blf_opt_ga], sig_b, ic95_b_low, ic95_b_high, [t_mc; t_fmc; t_ga], ...
    'VariableNames', {'Algorithm', 'Mass_kg', 'Sigma_Mass', 'Buckling_LF', 'Sigma_BLF', 'CI95_Min_B', 'CI95_Max_B', 'Time_s'});

fprintf('\n=========================================================================================\n');
fprintf('                MULTI-ALGORITHM BENCHMARK & UNCERTAINTY QUANTIFICATION (10D)            \n');
fprintf('=========================================================================================\n');
disp(T_alg);

%% OPTIMAL DESIGN VARIABLE CONFIGURATION BY ALGORITHM
alg_names = {'Monte Carlo'; 'fmincon (SQP)'; 'GA'};

fprintf('\n========================================================================================================================\n');
fprintf('                     OPTIMAL DESIGN VARIABLE CONFIGURATION BY ALGORITHM (10D)                                    \n');
fprintf('========================================================================================================================\n');
fprintf('%-15s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Algorithm', 'Phi_up', 'Psi_up', 'r_R_up', 'r_M_up', 'r_T_up', 'Phi_lo', 'Psi_lo', 'r_R_lo', 'r_M_lo', 'r_T_lo');
fprintf('------------------------------------------------------------------------------------------------------------------------\n');

for i = 1:size(X_opts, 1)
    x_i = X_opts(i, :);
    fprintf('%-15s %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n', ...
        alg_names{i}, x_i(1), x_i(2), x_i(3), x_i(4), x_i(5), x_i(6), x_i(7), x_i(8), x_i(9), x_i(10));
end
fprintf('========================================================================================================================\n\n');

%% 9. SPANWISE REPETITION PROFILES (ALL OPTIMA & FEASIBLE MEAN)
figure('Name', 'Spanwise Profiles Upper vs Lower Skin', 'Color', 'w', 'Position', [150 150 1000 600]);
span_pos = [0.0, 0.5, 1.0];

% Profile extraction (Upper skin: cols 3:5 | Lower skin: cols 8:10)
prof_up_m   = stable_X(idx_min_m, 3:5);       prof_lo_m   = stable_X(idx_min_m, 8:10);       % Min Mass
prof_up_b   = stable_X(idx_max_b, 3:5);       prof_lo_b   = stable_X(idx_max_b, 8:10);       % Max Buckling
prof_up_u   = stable_X(idx_min_u, 3:5);       prof_lo_u   = stable_X(idx_min_u, 8:10);       % Min Strain Energy
prof_up_k   = stable_X(idx_knee_stable, 3:5); prof_lo_k   = stable_X(idx_knee_stable, 8:10); % Knee Point
prof_up_avg = mean(stable_X(:, 3:5), 1);      prof_lo_avg = mean(stable_X(:, 8:10), 1);      % Feasible Mean

% Color palette definition
c_m   = [0.85 0.10 0.10]; % Red (Min Mass)
c_b   = [0.00 0.60 0.85]; % Cyan/Blue (Max Buckling)
c_u   = [0.75 0.00 0.75]; % Magenta (Min Strain Energy)
c_k   = [1.00 0.50 0.00]; % Orange (Knee Point)
c_avg = [0.35 0.35 0.35]; % Dark Gray (Mean)

hold on; grid on;

% Plot repetition profiles (Solid = Upper Skin | Dashed = Lower Skin)
plot(span_pos, prof_up_m,   '-s',  'LineWidth', 2, 'Color', c_m,   'MarkerFaceColor', c_m,   'MarkerSize', 8);
plot(span_pos, prof_lo_m,   '--s', 'LineWidth', 2, 'Color', c_m,   'MarkerFaceColor', 'w',   'MarkerSize', 8);

plot(span_pos, prof_up_b,   '-^',  'LineWidth', 2, 'Color', c_b,   'MarkerFaceColor', c_b,   'MarkerSize', 8);
plot(span_pos, prof_lo_b,   '--^', 'LineWidth', 2, 'Color', c_b,   'MarkerFaceColor', 'w',   'MarkerSize', 8);

plot(span_pos, prof_up_u,   '-o',  'LineWidth', 2, 'Color', c_u,   'MarkerFaceColor', c_u,   'MarkerSize', 8);
plot(span_pos, prof_lo_u,   '--o', 'LineWidth', 2, 'Color', c_u,   'MarkerFaceColor', 'w',   'MarkerSize', 8);

plot(span_pos, prof_up_k,   '-h',  'LineWidth', 2, 'Color', c_k,   'MarkerFaceColor', c_k,   'MarkerSize', 9);
plot(span_pos, prof_lo_k,   '--h', 'LineWidth', 2, 'Color', c_k,   'MarkerFaceColor', 'w',   'MarkerSize', 9);

plot(span_pos, prof_up_avg, '-v',  'LineWidth', 2, 'Color', c_avg, 'MarkerFaceColor', c_avg, 'MarkerSize', 7);
plot(span_pos, prof_lo_avg, '--v', 'LineWidth', 2, 'Color', c_avg, 'MarkerFaceColor', 'w',   'MarkerSize', 7);

all_vals = [prof_up_m, prof_lo_m, prof_up_b, prof_lo_b, prof_up_u, prof_lo_u, prof_up_k, prof_lo_k, prof_up_avg, prof_lo_avg];
ylim([min(all_vals) - 4.0, max(all_vals) + 4.5]);
xlim([-0.12, 1.12]);

% Data value annotations
dx = [-0.03, -0.015, 0.0, 0.015, 0.03]; 
for i = 1:length(span_pos)
    text(span_pos(i) + dx(1), prof_up_m(i) + 0.6, sprintf('%.0f', prof_up_m(i)), 'Color', c_m, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_pos(i) + dx(1), prof_lo_m(i) - 0.6, sprintf('%.0f', prof_lo_m(i)), 'Color', c_m, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    
    text(span_pos(i) + dx(2), prof_up_b(i) + 0.6, sprintf('%.0f', prof_up_b(i)), 'Color', c_b, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_pos(i) + dx(2), prof_lo_b(i) - 0.6, sprintf('%.0f', prof_lo_b(i)), 'Color', c_b, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

    text(span_pos(i) + dx(3), prof_up_u(i) + 0.6, sprintf('%.0f', prof_up_u(i)), 'Color', c_u, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_pos(i) + dx(3), prof_lo_u(i) - 0.6, sprintf('%.0f', prof_lo_u(i)), 'Color', c_u, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

    text(span_pos(i) + dx(4), prof_up_k(i) + 0.6, sprintf('%.0f', prof_up_k(i)), 'Color', c_k, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_pos(i) + dx(4), prof_lo_k(i) - 0.6, sprintf('%.0f', prof_lo_k(i)), 'Color', c_k, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

    text(span_pos(i) + dx(5), prof_up_avg(i) + 0.6, sprintf('%.1f', prof_up_avg(i)), 'Color', c_avg, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_pos(i) + dx(5), prof_lo_avg(i) - 0.6, sprintf('%.1f', prof_lo_avg(i)), 'Color', c_avg, 'FontSize', fs-4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
end

set(gca, 'XTick', span_pos, 'XTickLabel', {'Root ($x/L=0$)', 'Mid ($x/L=0.5$)', 'Tip ($x/L=1.0$)'}, 'TickLabelInterpreter', 'latex');
ylabel('Number of Repeats ($r_i$)', 'Interpreter', 'latex', 'FontSize', fs + 2); 
xlabel('Spanwise Location', 'Interpreter', 'latex', 'FontSize', fs + 2);
title('Repetitions Distribution Profile: Upper vs Lower Skin', 'Interpreter', 'latex', 'FontSize', fs + 3);

legend({'Min Mass (Upper)', 'Min Mass (Lower)', ...
        'Max Buckling (Upper)', 'Max Buckling (Lower)', ...
        'Min Strain Energy (Upper)', 'Min Strain Energy (Lower)', ...
        'Knee Point (Upper)', 'Knee Point (Lower)', ...
        'Mean Feasible (Upper)', 'Mean Feasible (Lower)'}, ...
       'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', fs - 3);
grid on; box on;

%% 10. BUCKLING FACTOR MAP WITH FEASIBILITY BOUNDARY
[grid_phi_u, grid_psi_u] = meshgrid(65:1:90, 0:1:20);
best_rep_10d = stable_X(idx_min_m, 3:10);

Z_blf_surf = zeros(size(grid_phi_u));

for r = 1:size(grid_phi_u, 1)
    for c = 1:size(grid_phi_u, 2)
        x_eval = [grid_phi_u(r,c), grid_psi_u(r,c), best_rep_10d];
        Z_blf_surf(r,c) = predict(gpr_blf, x_eval);
    end
end

figure('Name', 'Buckling Factor Map & Boundary Limit', 'Color', 'w', 'Position', [260 260 720 520]);
contourf(grid_phi_u, grid_psi_u, Z_blf_surf, 20, 'LineColor', 'none'); hold on;
contour(grid_phi_u, grid_psi_u, Z_blf_surf, [1.0 1.0], 'r-', 'LineWidth', 2.5);
colorbar; 
title('Buckling Factor $\lambda_b$ (Feasibility Boundary $\lambda_b=1.0$)', 'Interpreter', 'latex', 'FontSize', fs + 3);
xlabel('$\Phi_{\mathrm{up}}$ (deg)', 'Interpreter', 'latex', 'FontSize', fs + 2); 
ylabel('$\Psi_{\mathrm{up}}$ (deg)', 'Interpreter', 'latex', 'FontSize', fs + 2); grid on;

%% 11. ARD FEATURE SENSITIVITY & VARIABLE IMPORTANCE
% Extract hyperparameter length scales (theta_i) from GPR Matern 5/2 kernel
params_mass = gpr_mass.KernelInformation.KernelParameters(1:10);
params_blf  = gpr_blf.KernelInformation.KernelParameters(1:10);

% Compute normalized importance inverse-lengthscale scores (1 / theta_i)
imp_mass = (1 ./ params_mass) / sum(1 ./ params_mass);
imp_blf  = (1 ./ params_blf) / sum(1 ./ params_blf);

figure('Name', 'ARD Predictor Importance Analysis', 'Color', 'w', 'Position', [270 270 900 500]);
b = bar([imp_mass, imp_blf]);
b(1).FaceColor = [0.2 0.6 0.9];
b(2).FaceColor = [0.9 0.3 0.3];
set(gca, 'XTick', 1:10, 'XTickLabel', var_names, 'TickLabelInterpreter', 'latex');
ylabel('Normalized Importance Score ($1/\theta_i$)', 'Interpreter', 'latex', 'FontSize', fs + 2);
title('GPR ARD Predictor Sensitivity \& Feature Importance', 'Interpreter', 'latex', 'FontSize', fs + 3);
legend({'Total Mass Model', 'Buckling Factor ($\lambda_b$) Model'}, 'Location', 'northeast', 'Interpreter', 'latex');
grid on;

%% 12. GLOBAL UNCERTAINTY DISTRIBUTION & ADAPTIVE SAMPLING (PHASE 4)
[~, sigma_blf_all] = predict(gpr_blf, X_mc_valid);
[~, idx_uncert] = sort(sigma_blf_all, 'descend');

% Select top 5 points with maximum predictive variance for Phase 4 adaptive sampling
N_next = 5;
X_next = X_mc_valid(idx_uncert(1:N_next), :);

figure('Name', 'Global Uncertainty Distribution 10D', 'Color', 'w', 'Position', [290 290 720 500]);
histogram(sigma_blf_all, 25, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none'); hold on;
xline(mean(sigma_blf_all), 'r--', 'LineWidth', 2, ...
    'Label', sprintf('Mean: %.4f', mean(sigma_blf_all)), ...
    'Interpreter', 'latex', ...
    'FontSize', fs, ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'top', ...
    'LabelHorizontalAlignment', 'right');
xlabel('Predictive Buckling Load Factor Uncertainty $\sigma_{\lambda_b}$', 'Interpreter', 'latex', 'FontSize', fs + 2); 
ylabel('Frequency (Feasible MC Domain)', 'Interpreter', 'latex', 'FontSize', fs + 2);
title('GPR Predictive Buckling Factor Uncertainty Distribution', 'Interpreter', 'latex', 'FontSize', fs + 3); grid on;

fprintf('\n====================================================================================================\n');
fprintf('                     NEW CANDIDATE SAMPLES FOR PHASE 4 (MAXIMUM UNCERTAINTY ADAPTIVE SAMPLING)\n');
fprintf('====================================================================================================\n');
fprintf('%-4s %-7s %-7s %-6s %-6s %-6s %-7s %-7s %-6s %-6s %-6s\n', ...
    'Pt', 'Phi_up', 'Psi_up', 'r_R_u', 'r_M_u', 'r_T_u', 'Phi_lo', 'Psi_lo', 'r_R_l', 'r_M_l', 'r_T_l');
fprintf('----------------------------------------------------------------------------------------------------\n');
for i = 1:N_next
    p_num = n_samples + i;
    fprintf('%-4d %-7.1f %-7.1f %-6d %-6d %-6d %-7.1f %-7.1f %-6d %-6d %-6d\n', ...
        p_num, X_next(i,1), X_next(i,2), X_next(i,3), X_next(i,4), X_next(i,5), ...
        X_next(i,6), X_next(i,7), X_next(i,8), X_next(i,9), X_next(i,10));
end
fprintf('====================================================================================================\n');

%% OPTIMAL DESIGNS AND KNEE POINT FOR HYPERMESH EXPORT
opt_indices = [idx_min_m, idx_min_u, idx_max_b, idx_knee_stable];
opt_types   = {'MIN_MASS', 'MIN_STRAIN', 'MAX_BUCKL', 'KNEE_POINT'};

fprintf('\n========================================================================================================================\n');
fprintf('                                   OPTIMAL DESIGNS FOR HYPERMESH EXPORT                                                 \n');
fprintf('========================================================================================================================\n');
fprintf('%-6s %-12s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Pt', 'Type', 'Phi_up', 'Psi_up', 'r_R_up', 'r_M_up', 'r_T_up', 'Phi_lo', 'Psi_lo', 'r_R_lo', 'r_M_lo', 'r_T_lo');
fprintf('------------------------------------------------------------------------------------------------------------------------\n');

for k = 1:length(opt_indices)
    idx = opt_indices(k);
    x_k = stable_X(idx, :);
    fprintf('%-6d %-12s %-8.1f %-8.1f %-8d %-8d %-8d %-8.1f %-8.1f %-8d %-8d %-8d\n', ...
        k, opt_types{k}, x_k(1), x_k(2), x_k(3), x_k(4), x_k(5), x_k(6), x_k(7), x_k(8), x_k(9), x_k(10));
end
fprintf('========================================================================================================================\n\n');
