% =========================================================================
% PHASE 1: MULTI-RESPONSE SURROGATE MODELING & DIAGNOSTICS (2D: PHI, PSI)
% Laminate Repetitions: [r_Root, r_Mid, r_Tip] = [26, 21, 15] | Mass = 420.52 kg
% =========================================================================
clear; clc; close all;

% Set global LaTeX formatting and enlarged font sizes for high-res figures
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16); % Increased default axis font size
set(0, 'defaultAxesTitleFontSizeMultiplier', 1.25);

%% 1. INPUT DATA LOADING (20 LHS SAMPLES & OPTISTRUCT FEA RESULTS)
X_p1 = [
    35.0, 21.1; 29.5, 32.8; 15.6,  4.7; 63.0, 67.1;  5.8, 88.7;
    72.1, 69.3; 22.3, 72.8; 46.6, 23.5;  9.2, 38.0;  3.3, 53.8;
    76.7, 42.0; 22.9, 81.7; 44.1, 80.7; 37.2, 29.9; 57.7, 48.3;
    85.0, 61.1; 52.6, 57.1; 70.4, 13.7; 88.2,  3.9; 61.7, 12.5
];
Y_p1 = [
    0.0413158, 3838.44154, 0.9094562;
    0.0441925, 4120.07722, 0.9421997;
    0.0396258, 3644.81819, 0.6038449;
    0.0370649, 3308.41500, 0.8623080;
    0.0295141, 2452.58439, 1.1324130;
    0.0360718, 3127.24435, 0.7658765;
    0.0294471, 2613.83246, 1.3310430;
    0.0388523, 3600.53777, 1.1464970;
    0.0375091, 3463.22367, 0.8972938;
    0.0315553, 2866.68528, 1.1438510;
    0.0334588, 2996.34284, 1.1489320;
    0.0288854, 2521.23342, 1.2949500;
    0.0333188, 2965.85828, 1.0972300;
    0.0443937, 4141.11429, 1.0035030;
    0.0422049, 3901.51453, 1.0396280;
    0.0353584, 3057.03525, 0.8044214;
    0.0422226, 3895.62705, 1.0064910;
    0.0287040, 2531.11677, 1.2555540;
    0.0296652, 2458.04058, 1.1257390;
    0.0301032, 2708.12229, 1.2455460
];

% Raw OptiStruct responses in SI units
y_disp_m   = Y_p1(:, 1); % Max displacement (m)
y_strain_J = Y_p1(:, 2); % Strain Energy (J)
y_blf      = Y_p1(:, 3); % Buckling Load Factor (dimensionless)

%% 2. GAUSSIAN PROCESS REGRESSION (GPR) FITTING
gpr_disp   = fitrgp(X_p1, y_disp_m,   'KernelFunction', 'Matern52', 'Standardize', true);
gpr_strain = fitrgp(X_p1, y_strain_J, 'KernelFunction', 'Matern52', 'Standardize', true);
gpr_blf    = fitrgp(X_p1, y_blf,    'KernelFunction', 'Matern52', 'Standardize', true);

%% 3. LEAVE-ONE-OUT CROSS-VALIDATION (LOOCV)
cv_disp   = crossval(gpr_disp,   'KFold', length(y_disp_m));
cv_strain = crossval(gpr_strain, 'KFold', length(y_strain_J));
cv_blf    = crossval(gpr_blf,    'KFold', length(y_blf));

p_disp_m   = kfoldPredict(cv_disp);
p_strain_J = kfoldPredict(cv_strain);
p_blf      = kfoldPredict(cv_blf);

% Performance metric calculation handles
calc_r2   = @(y,p) 1 - sum((y - p).^2) / sum((y - mean(y)).^2);
calc_rmse = @(y,p) sqrt(mean((y - p).^2));

r2_d = calc_r2(y_disp_m, p_disp_m);
r2_s = calc_r2(y_strain_J, p_strain_J);
r2_b = calc_r2(y_blf, p_blf);

rmse_d_m = calc_rmse(y_disp_m, p_disp_m);
rmse_s_J = calc_rmse(y_strain_J, p_strain_J);
rmse_b   = calc_rmse(y_blf, p_blf);

% Convert displacement quantities to millimeters for display
y_disp_mm = 1e3 * y_disp_m;
p_disp_mm = 1e3 * p_disp_m;
rmse_d_mm = 1e3 * rmse_d_m;

fprintf('=========================================================\n');
fprintf('         PHASE 1 METAMODEL LOOCV VALIDATION METRICS       \n');
fprintf('=========================================================\n');
fprintf('Max Displacement    -> R^2: %.4f | RMSE: %.3f mm\n', r2_d, rmse_d_mm);
fprintf('Strain Energy       -> R^2: %.4f | RMSE: %.2f J\n', r2_s, rmse_s_J);
fprintf('Buckling Factor BLF -> R^2: %.4f | RMSE: %.4f\n', r2_b, rmse_b);
fprintf('=========================================================\n\n');

%% 4. CONTINUOUS GRID PREDICTION (91x91 MESH)
[grid_phi, grid_psi] = meshgrid(0:1:90, 0:1:90);
X_grid = [grid_phi(:), grid_psi(:)];

Z_disp_m = reshape(predict(gpr_disp, X_grid), size(grid_phi));
Z_disp_mm = 1e3 * Z_disp_m;
Z_strain_J = reshape(predict(gpr_strain, X_grid), size(grid_phi));
[Z_blf_pred, Z_blf_sd] = predict(gpr_blf, X_grid);
Z_blf = reshape(Z_blf_pred, size(grid_phi));
Z_sd  = reshape(Z_blf_sd, size(grid_phi));

%% 5. SURROGATE RESPONSE SURFACES (INDEPENDENT FIGURES)
% Figure 1: Max Displacement Surface
figure('Name', 'Surrogate Surface - Max Displacement', 'Position', [100, 100, 700, 550]);
contourf(grid_phi, grid_psi, Z_disp_mm, 20, 'LineColor', 'none');
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$u_{\max}$ (mm)';
cb.Label.FontSize = 16;
title('Max Displacement $u_{\max}$ (mm)', 'FontSize', 20);
xlabel('$\Phi$ (deg)', 'FontSize', 18); ylabel('$\Psi$ (deg)', 'FontSize', 18);
xlim([0 90]); ylim([0 90]); xticks(0:10:90); yticks(0:10:90); grid on;

% Figure 2: Strain Energy Surface
figure('Name', 'Surrogate Surface - Strain Energy', 'Position', [150, 100, 700, 550]);
contourf(grid_phi, grid_psi, Z_strain_J, 20, 'LineColor', 'none');
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$U$ (J)';
cb.Label.FontSize = 16;
title('Strain Energy $U$ (J)', 'FontSize', 20);
xlabel('$\Phi$ (deg)', 'FontSize', 18); ylabel('$\Psi$ (deg)', 'FontSize', 18);
xlim([0 90]); ylim([0 90]); xticks(0:10:90); yticks(0:10:90); grid on;

% Figure 3: Buckling Factor Surface
figure('Name', 'Surrogate Surface - Buckling Factor', 'Position', [150, 100, 600, 500]);
contourf(grid_phi, grid_psi, Z_blf, 20, 'LineColor', 'none'); 
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\lambda_b$';
cb.Label.FontSize = 16; hold on;
[~, h_blf] = contour(grid_phi, grid_psi, Z_blf, [1.0 1.0], 'r-', 'LineWidth', 2.5);
title('Buckling Factor $\lambda_b$', 'FontSize', 20);
xlabel('$\Phi$ (deg)', 'FontSize', 18); ylabel('$\Psi$ (deg)', 'FontSize', 18);
xlim([0 90]); ylim([0 90]); xticks(0:10:90); yticks(0:10:90);
legend(h_blf, 'Boundary $\lambda_b = 1.0$', 'Location', 'bestoutside', 'FontSize', 14);
grid on;

%% 6. PARITY PLOTS (INDEPENDENT FIGURES)
% Figure 4: Parity Displacement
figure('Name', 'Parity Plot - Displacement', 'Position', [250, 100, 650, 500]);
plot(y_disp_mm, p_disp_mm, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 9);
hold on;
lims = [min(y_disp_mm), max(y_disp_mm)];
plot(lims, lims, 'k--', 'LineWidth', 2);
xlabel('FEA True $u_{\max}$ (mm)', 'FontSize', 18);
ylabel('GPR Predicted $u_{\max}$ (mm)', 'FontSize', 18);
title(sprintf('Displacement ($R^2 = %.3f$)', r2_d), 'FontSize', 20);
grid on;

% Figure 5: Parity Strain Energy
figure('Name', 'Parity Plot - Strain Energy', 'Position', [300, 100, 650, 500]);
plot(y_strain_J, p_strain_J, 'sq', 'MarkerFaceColor', 'm', 'MarkerSize', 9);
hold on;
lims = [min(y_strain_J), max(y_strain_J)];
plot(lims, lims, 'k--', 'LineWidth', 2);
xlabel('FEA True Strain Energy (J)', 'FontSize', 18);
ylabel('GPR Predicted Strain Energy (J)', 'FontSize', 18);
title(sprintf('Strain Energy ($R^2 = %.3f$)', r2_s), 'FontSize', 20);
grid on;

% Figure 6: Parity Buckling Factor
figure('Name', 'Parity Plot - Buckling Factor', 'Position', [350, 100, 650, 500]);
plot(y_blf, p_blf, '^', 'MarkerFaceColor', 'r', 'MarkerSize', 9); hold on;
lims = [min(y_blf), max(y_blf)];
plot(lims, lims, 'k--', 'LineWidth', 2);
title(sprintf('Buckling Factor ($R^2 = %.3f$)', r2_b), 'FontSize', 20);
xlabel('FEA True $\lambda_b$', 'FontSize', 18); 
ylabel('GPR Predicted $\lambda_b$', 'FontSize', 18);
grid on;

%% 7. PREDICTION UNCERTAINTY MAP
figure('Name', 'Uncertainty Map', 'Position', [150, 100, 600, 500]);
contourf(grid_phi, grid_psi, Z_sd, 20, 'LineColor', 'none'); 
hold on;
cb = colorbar;
cb.Label.Interpreter = 'latex';
ylabel(cb, '$\sigma$', 'FontSize', 16, 'Interpreter', 'latex');
h_pts = plot(X_p1(:,1), X_p1(:,2), 'r+', 'MarkerSize', 14, 'LineWidth', 2.5);
title('Buckling Uncertainty $\sigma(\Phi, \Psi)$', 'FontSize', 20, 'Interpreter', 'latex');
xlabel('$\Phi$ (deg)', 'FontSize', 18, 'Interpreter', 'latex'); 
ylabel('$\Psi$ (deg)', 'FontSize', 18, 'Interpreter', 'latex');
% Legend positioned horizontally at the bottom outside
legend(h_pts, 'LHS FEA Samples', 'Location', 'bestoutside', 'Orientation', 'horizontal', 'FontSize', 14);
grid on;

%% 8. CONSTRAINED OPTIMAL SEARCH
figure('Name', 'Constrained Optimal Search', 'Position', [450, 100, 850, 580]);
strain_levels = linspace(min(Z_strain_J(:)), max(Z_strain_J(:)), 20);
contourf(grid_phi, grid_psi, Z_strain_J, strain_levels, 'LineColor', 'none');
hold on; 

[~, h_bound] = contour(grid_phi, grid_psi, Z_blf, [1.0 1.0], 'r-', 'LineWidth', 2.5);

% Mask regions violating the buckling constraint (lambda_b < 1.0)
Z_strain_constrained = Z_strain_J;
Z_strain_constrained(Z_blf < 1.0) = NaN;
[min_strain_val, min_idx] = min(Z_strain_constrained(:), [], 'omitnan');
opt_phi = grid_phi(min_idx);
opt_psi = grid_psi(min_idx);
h_opt = plot(opt_phi, opt_psi, 'p', 'MarkerSize', 18, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');

cb = colorbar;
cb.Label.String = 'Strain Energy $U$ (J)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 16;
clim([min(Z_strain_J(:)), max(Z_strain_J(:))]);

title('Strain Energy $U$ (J) with Stability Boundary', 'FontSize', 20, 'Interpreter', 'latex');
xlabel('$\Phi$ (deg)', 'FontSize', 18, 'Interpreter', 'latex'); 
ylabel('$\Psi$ (deg)', 'FontSize', 18, 'Interpreter', 'latex');
xlim([0 90]); ylim([0 90]); 
xticks(0:10:90); yticks(0:10:90);

% Legend positioned horizontally at the bottom outside
legend_text = {'Boundary $\lambda_b = 1.0$', sprintf('Optimum: (%.1f$^\\circ$, %.1f$^\\circ$)', opt_phi, opt_psi)};
legend([h_bound, h_opt], legend_text, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 14, 'Interpreter', 'latex');
colormap('parula');
grid on;

fprintf('--- PHASE 1 LOCAL OPTIMIZATION FINDINGS ---\n');
fprintf('Optimal Angles (Phi*, Psi*) : (%.1f deg, %.1f deg)\n', opt_phi, opt_psi);
fprintf('Minimum Strain Energy       : %.2f J\n', min_strain_val);
fprintf('Predicted Buckling Factor    : %.3f\n\n', Z_blf(min_idx));

%% 9. 1D PARAMETRIC SLICE WITH 95% CONFIDENCE INTERVAL
psi_fixed = 30;
phi_sweep = (0:0.2:90)';
X_sweep = [phi_sweep, repmat(psi_fixed, length(phi_sweep), 1)];
[y_sweep_mean, y_sweep_sd] = predict(gpr_blf, X_sweep);
c_upper = y_sweep_mean + 1.96 * y_sweep_sd;
c_lower = y_sweep_mean - 1.96 * y_sweep_sd;

figure('Name', '1D Parametric Slice', 'Position', [500, 100, 850, 500]);
h_fill = fill([phi_sweep; flipud(phi_sweep)], [c_upper; flipud(c_lower)], [0.8 0.85 1], 'LineStyle', 'none'); hold on;
h_line = plot(phi_sweep, y_sweep_mean, 'b-', 'LineWidth', 2.5);
h_crit = yline(1.0, 'r--', 'LineWidth', 2);
title(sprintf('1D Parametric Slice at $\\Psi = %d^\\circ$ (95\\%% CI)', psi_fixed), 'FontSize', 20);
xlabel('$\Phi$ (deg)', 'FontSize', 18); ylabel('Buckling Factor $\lambda_b$', 'FontSize', 18);
xlim([0 90]); xticks(0:10:90);
legend([h_line, h_fill, h_crit], {'GPR Mean Prediction', '95\% Confidence Interval', 'Critical Limit $\lambda_b = 1.0$'}, ...
    'Location', 'northeastoutside', 'FontSize', 14);
grid on;