% =========================================================================
% PHASE 2: INTEGRAL SURROGATE MODELING, DIAGNOSTICS & OPTIMIZATION
% =========================================================================
clear; clc; close all;

%% Global LaTeX Interpreter Configuration for Graphics
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 10);
set(0, 'defaultAxesTitleFontSizeMultiplier', 1.15);

%% 1. DATA LOADING (50 LHS SAMPLES & OPTISTRUCT FEA)
X_p2 = [
    56.4, 66.4, 23.53, 21.92, 11.72;  23.7, 16.7, 29.75, 28.27, 26.13;
    72.7,  8.9, 29.30, 18.69, 13.02;  66.4,  0.3, 16.19, 15.27, 13.44;
    28.4, 23.7, 28.03, 25.46, 15.01;  25.6, 19.8, 27.57, 17.99, 11.60;
    63.8, 20.3, 25.81, 20.46, 12.25;  33.3, 36.2, 24.36, 23.49, 16.79;
    88.7, 32.3, 22.28, 20.58, 11.38;  49.6, 73.4, 27.07, 22.25, 12.78;
    9.3,  42.9, 22.88, 22.76, 20.06;  69.8, 79.0, 21.54, 20.58, 17.93;
    54.8, 65.1, 28.81, 25.45, 15.69;  12.6, 58.3, 19.77, 17.49, 13.35;
    21.5, 78.2, 28.59, 12.01, 10.13;  85.5,  6.0, 21.46, 19.32, 10.92;
    19.5,  2.6, 26.99, 22.17, 21.87;   5.0, 40.5, 29.71, 14.83, 13.82;
    41.3, 30.9, 24.60, 18.44, 12.68;   1.2, 54.6, 26.27, 16.16, 11.45;
    42.9, 26.8, 29.48, 24.35, 18.99;  33.8, 13.8, 21.23, 20.34, 13.73;
    77.6, 49.1, 27.91, 27.52, 10.02;  46.3, 50.7, 26.23, 25.41, 15.55;
    7.2,  89.4, 16.64, 14.89, 13.98;  53.7, 84.3, 19.16, 10.67, 10.28;
    13.9, 60.8, 29.31, 25.80, 23.19;  86.1, 76.2, 26.89, 20.51, 15.39;
    3.2,  34.0, 18.15, 16.17, 15.71;  62.3, 12.4, 24.91, 24.63, 11.43;
    68.4, 70.7, 29.55, 24.45, 14.36;  59.5, 54.0, 27.25, 19.47, 12.96;
    31.0, 86.7, 28.84, 23.81, 18.54;  80.9, 68.1, 19.57, 17.66, 11.56;
    40.0, 82.0, 21.19, 17.29, 14.94;  47.3,  9.0, 25.82, 14.12, 10.84;
    16.2, 59.6, 23.76, 16.50, 12.18;  76.4, 29.1, 28.35, 26.62, 23.19;
    37.0, 47.0, 27.20, 19.83, 14.34;  82.7, 41.7, 22.81, 18.42, 17.39;
    89.0, 86.2, 29.03, 26.87, 25.60;  80.7, 81.2, 29.34, 27.17, 25.68;
    68.2, 81.2, 29.40, 27.68, 25.85;  14.0, 7.3, 15.26, 12.61, 10.01;
    11.8,   1.0, 15.08, 12.89, 10.73;      86.6,  88.2, 28.73, 11.57, 10.40;
      3.6,  67.7, 28.68, 27.00, 10.27;
     86.9,  68.8, 27.61, 10.89, 10.25;
      3.6,  89.4, 29.22, 23.92, 10.86;
      6.6,  89.9, 28.50, 25.75, 11.39;
];

resTable = readtable('C:\Users\Hugo\OneDrive\Escritorio\surrogate_2\OptiStruct_KPI_Resultados_Phase2_v3.txt');
y_mass   = resTable.Masa_Total;
y_disp   = resTable.Desplazamiento_Max;
y_strain = resTable.Energia_Deformacion;
y_blf    = resTable.Factor_Pandeo_BLF;

lb = [0,  0,  15, 10, 10];
ub = [90, 90, 30, 28, 26];

%% 2. GPR SURROGATE MODELS WITH ARD AND LINEAR BASIS FUNCTIONS
opts_ard = {'KernelFunction', 'ARDMatern52', 'BasisFunction', 'linear', 'Standardize', true};

% Automatic hyperparameter optimization options
opts_auto = { ...
    'OptimizeHyperparameters', {'KernelFunction', 'KernelScale', 'Sigma'}, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'MaxObjectiveEvaluations', 30, ...
        'Verbose', 0) ...
};

gpr_blf = fitrgp(X_p2, y_blf, 'BasisFunction', 'linear', 'Standardize', true, opts_auto{:});

gpr_mass   = fitrgp(X_p2, y_mass, opts_ard{:});
gpr_disp   = fitrgp(X_p2, y_disp, opts_ard{:});
gpr_strain = fitrgp(X_p2, y_strain, opts_ard{:});

%% 3. LOOCV EVALUATION (LEAVE-ONE-OUT CROSS-VALIDATION)
n_samples = height(resTable);
calc_r2   = @(y, p) 1 - sum((y - p).^2) / sum((y - mean(y)).^2);
calc_rmse = @(y, p) sqrt(mean((y - p).^2));

p_mass   = kfoldPredict(crossval(gpr_mass, 'KFold', n_samples));
p_disp   = kfoldPredict(crossval(gpr_disp, 'KFold', n_samples));
p_strain = kfoldPredict(crossval(gpr_strain, 'KFold', n_samples));
p_blf    = kfoldPredict(crossval(gpr_blf, 'KFold', n_samples));

r2_m = calc_r2(y_mass, p_mass);     rmse_m = calc_rmse(y_mass, p_mass);
r2_d = calc_r2(y_disp, p_disp);     rmse_d = calc_rmse(y_disp*1e3, p_disp*1e3);
r2_s = calc_r2(y_strain, p_strain); rmse_s = calc_rmse(y_strain, p_strain);
r2_b = calc_r2(y_blf, p_blf);       rmse_b = calc_rmse(y_blf, p_blf);

%% 4. MONTE CARLO SEARCH & FEASIBLE DESIGN SPACE
N_mc = 20000;
rng(42, 'twister'); % Fix pseudorandom seed for strict reproducibility
X_mc = [rand(N_mc,1)*90, rand(N_mc,1)*90, ...
        lb(3)+rand(N_mc,1)*(ub(3)-lb(3)), ...
        lb(4)+rand(N_mc,1)*(ub(4)-lb(4)), ...
        lb(5)+rand(N_mc,1)*(ub(5)-lb(5))];

% Apply physical geometric constraints (r_Root >= r_Mid >= r_Tip)
valid_idx = (X_mc(:,3) >= X_mc(:,4)) & (X_mc(:,4) >= X_mc(:,5));
X_mc = X_mc(valid_idx, :);

mc_mass   = predict(gpr_mass, X_mc);
mc_strain = predict(gpr_strain, X_mc);
mc_blf    = predict(gpr_blf, X_mc);

stable_idx  = mc_blf >= 1.0;
stable_X    = X_mc(stable_idx, :);
stable_mass = mc_mass(stable_idx);
stable_blf  = mc_blf(stable_idx);
stable_U    = mc_strain(stable_idx);

[min_mass_val, idx_min_m]      = min(stable_mass);
[min_strain_val, idx_min_u]    = min(stable_U);
[max_stability_val, idx_max_b] = max(stable_blf);

p0_r_root = min(stable_X(:,3));  p5_r_root = prctile(stable_X(:,3), 5);  p50_r_root = prctile(stable_X(:,3), 50);
p0_r_mid  = min(stable_X(:,4));  p5_r_mid  = prctile(stable_X(:,4), 5);  p50_r_mid  = prctile(stable_X(:,4), 50);
p0_r_tip  = min(stable_X(:,5));  p5_r_tip  = prctile(stable_X(:,5), 5);  p50_r_tip  = prctile(stable_X(:,5), 50);

%% 5. CONSOLE OUTPUT
fprintf('\n=================================================================\n');
fprintf('                SURROGATE MODEL METRICS (LOOCV)                 \n');
fprintf('=================================================================\n');
fprintf('Mass       -> R2: %.4f | RMSE: %.4f kg\n', r2_m, rmse_m);
fprintf('Disp Max   -> R2: %.4f | RMSE: %.4f mm\n', r2_d, rmse_d);
fprintf('Strain E.  -> R2: %.4f | RMSE: %.4f J\n', r2_s, rmse_s);
fprintf('Buckling L.-> R2: %.4f | RMSE: %.4f\n', r2_b, rmse_b);
fprintf('-----------------------------------------------------------------\n');
fprintf('Monte Carlo Valid Samples  : %d\n', length(X_mc));
fprintf('Monte Carlo Feasible Points: %d / %d (%.2f%%)\n', sum(stable_idx), length(X_mc), (sum(stable_idx)/length(X_mc))*100);
fprintf('-----------------------------------------------------------------\n');
fprintf('OPTIMAL DESIGN CANDIDATES:\n');
fprintf('1. Min Mass Design       : Mass = %.2f kg | Buckling = %.3f | Strain E. = %.2f J\n', ...
        min_mass_val, stable_blf(idx_min_m), stable_U(idx_min_m));
fprintf('   Params: Phi = %.1f deg, Psi = %.1f deg, r_Root = %.2f, r_Mid = %.2f, r_Tip = %.2f\n', ...
        stable_X(idx_min_m,1), stable_X(idx_min_m,2), stable_X(idx_min_m,3), stable_X(idx_min_m,4), stable_X(idx_min_m,5));

fprintf('2. Min Strain E. Design  : Mass = %.2f kg | Buckling = %.3f | Strain E. = %.2f J\n', ...
        stable_mass(idx_min_u), stable_blf(idx_min_u), min_strain_val);
fprintf('   Params: Phi = %.1f deg, Psi = %.1f deg, r_Root = %.2f, r_Mid = %.2f, r_Tip = %.2f\n', ...
        stable_X(idx_min_u,1), stable_X(idx_min_u,2), stable_X(idx_min_u,3), stable_X(idx_min_u,4), stable_X(idx_min_u,5));

fprintf('3. Max Buckling Design   : Mass = %.2f kg | Buckling = %.3f | Strain E. = %.2f J\n', ...
        stable_mass(idx_max_b), max_stability_val, stable_U(idx_max_b));
fprintf('   Params: Phi = %.1f deg, Psi = %.1f deg, r_Root = %.2f, r_Mid = %.2f, r_Tip = %.2f\n', ...
        stable_X(idx_max_b,1), stable_X(idx_max_b,2), stable_X(idx_max_b,3), stable_X(idx_max_b,4), stable_X(idx_max_b,5));

fprintf('-----------------------------------------------------------------\n');
fprintf('STABILITY THRESHOLDS (lambda_b >= 1.0):\n');
fprintf('r_Root -> Min: %.2f | P5 Safe: %.2f | P50 Median: %.2f\n', p0_r_root, p5_r_root, p50_r_root);
fprintf('r_Mid  -> Min: %.2f | P5 Safe: %.2f | P50 Median: %.2f\n', p0_r_mid, p5_r_mid, p50_r_mid);
fprintf('r_Tip  -> Min: %.2f | P5 Safe: %.2f | P50 Median: %.2f\n', p0_r_tip, p5_r_tip, p50_r_tip);
fprintf('=================================================================\n\n');

%% =========================================================================
%% FIGURES OF THE STUDY
%% =========================================================================

%% FIGURE 1: LOOCV PARITY PLOTS
figure('Name', 'Fig 1: LOOCV Parity Plots', 'Position', [80, 80, 1000, 680]);
subplot(2,2,1); plot(y_mass, p_mass, 'bo', 'MarkerFaceColor', [0.2 0.6 1], 'MarkerSize', 6); hold on;
plot([min(y_mass) max(y_mass)], [min(y_mass) max(y_mass)], 'k--', 'LineWidth', 1.5);
title(sprintf('Total Mass ($R^2 = %.4f$)', r2_m)); xlabel('FEA Mass (kg)'); ylabel('GPR Mass (kg)'); grid on;

subplot(2,2,2); plot(y_disp*1e3, p_disp*1e3, 'sq', 'MarkerFaceColor', [0.9 0.4 0.2], 'MarkerSize', 6); hold on;
plot([min(y_disp) max(y_disp)]*1e3, [min(y_disp) max(y_disp)]*1e3, 'k--', 'LineWidth', 1.5);
title(sprintf('Max Displacement ($R^2 = %.4f$)', r2_d)); xlabel('FEA $u_{\\max}$ (mm)'); ylabel('GPR $u_{\\max}$ (mm)'); grid on;

subplot(2,2,3); plot(y_strain, p_strain, 'd', 'MarkerFaceColor', [0.4 0.8 0.3], 'MarkerSize', 6); hold on;
plot([min(y_strain) max(y_strain)], [min(y_strain) max(y_strain)], 'k--', 'LineWidth', 1.5);
title(sprintf('Strain Energy ($R^2 = %.4f$)', r2_s)); xlabel('FEA Strain Energy $U$ (J)'); ylabel('GPR Strain Energy $U$ (J)'); grid on;

subplot(2,2,4); plot(y_blf, p_blf, '^', 'MarkerFaceColor', [0.8 0.2 0.8], 'MarkerSize', 6); hold on;
plot([min(y_blf) max(y_blf)], [min(y_blf) max(y_blf)], 'k--', 'LineWidth', 1.5);
title(sprintf('Buckling Load Factor ($R^2 = %.4f$)', r2_b)); xlabel('FEA $\lambda_b$'); ylabel('GPR $\lambda_b$'); grid on;

%% FIGURE 2A: RESPONSE SURFACES - STRAIN ENERGY U
figure('Name', 'Fig 2A: Response Surfaces Strain Energy', 'Position', [100, 100, 1100, 360]);
[grid_phi, grid_psi] = meshgrid(0:2:90, 0:2:90);
taper_configs_3 = [16, 12, 10; 24, 18, 13; 29, 26, 23];
taper_titles_3  = {'Thin $[16, 12, 10]$', 'Nominal $[24, 18, 13]$', 'Thick $[29, 26, 23]$'};

for t = 1:3
    subplot(1,3,t);
    X_g2D = [grid_phi(:), grid_psi(:), repmat(taper_configs_3(t,:), numel(grid_phi), 1)];
    Z_U = reshape(predict(gpr_strain, X_g2D), size(grid_phi));
    contourf(grid_phi, grid_psi, Z_U, 20, 'LineColor', 'none'); cb = colorbar;
    cb.Label.Interpreter = 'latex'; cb.Label.String = '$U$ (J)';
    title(['Strain Energy $U$ - ' taper_titles_3{t}]);
    xlabel('$\Phi$ (deg)'); ylabel('$\Psi$ (deg)'); grid on;
end

%% FIGURE 2B: RESPONSE SURFACES - BUCKLING LOAD FACTOR (CENTERED TEXT AT CENTROID)
figure('Name', 'Fig 2B: Response Surfaces Buckling Factor (4 Configs)', 'Position', [50, 150, 1400, 350]);
taper_configs_4 = [15, 12, 10; 20, 16, 12; 25, 20, 15; 29, 26, 23];
taper_titles_4  = {'Thin $[15, 12, 10]$', 'Low-Interm $[20, 16, 12]$', ...
                   'High-Interm $[25, 20, 15]$', 'Thick $[29, 26, 23]$'};

for t = 1:4
    subplot(1,4,t);
    X_g2D = [grid_phi(:), grid_psi(:), repmat(taper_configs_4(t,:), numel(grid_phi), 1)];
    Z_BLF = reshape(predict(gpr_blf, X_g2D), size(grid_phi));
    contourf(grid_phi, grid_psi, Z_BLF, 20, 'LineColor', 'none'); hold on; cb = colorbar;
    cb.Label.Interpreter = 'latex'; cb.Label.String = '$\lambda_b$';
    
    [C_b, h_b2] = contour(grid_phi, grid_psi, Z_BLF, [1.0 1.0], 'r-', 'LineWidth', 2.5);
    title(['$\lambda_b$ - ' taper_titles_4{t}]);
    xlabel('$\Phi$ (deg)'); ylabel('$\Psi$ (deg)'); grid on;
    
    feasible_mask = Z_BLF >= 1.0;
    if any(feasible_mask(:)) && ~all(feasible_mask(:))
        legend(h_b2, 'Limit $\lambda_b = 1.0$', 'Location', 'northeast');
        % Geometric centroid calculation of the feasible region
        x_center = mean(grid_phi(feasible_mask));
        y_center = mean(grid_psi(feasible_mask));
        text(x_center, y_center, '\textbf{Feasible Region}', ...
             'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center', ...
             'BackgroundColor', [1 1 1 0.75], 'Interpreter', 'latex');
    elseif all(feasible_mask(:))
        text(45, 45, {'\textbf{Fully Stable Domain:}', 'All $\lambda_b \ge 1.0$'}, 'HorizontalAlignment', 'center', ...
             'Color', 'k', 'FontSize', 8.5, 'BackgroundColor', [0.2 0.8 0.3 0.6], 'Interpreter', 'latex');
    else
        text(45, 45, {'\textbf{Unstable Domain:}', 'All $\lambda_b < 1.0$'}, 'HorizontalAlignment', 'center', ...
             'Color', 'w', 'FontSize', 8.5, 'BackgroundColor', [0.7 0.1 0.1 0.6], 'Interpreter', 'latex');
    end
end

%% FIGURE 3A: MAIN EFFECTS OF ANGLES ON BUCKLING (FIXED ANGLES AND REPETITIONS IN TITLES)
figure('Name', 'Fig 3A: Main Effects Angles', 'Position', [150, 150, 800, 420]);
base_design = [35.0, 60.0, 25.0, 20.0, 15.0];

for v = 1:2
    subplot(1, 2, v);
    sweep = linspace(lb(v), ub(v), 60)';
    X_sw = repmat(base_design, 60, 1); 
    X_sw(:, v) = sweep;
    [y_sw, sd_sw] = predict(gpr_blf, X_sw);
    
    h_fill = fill([sweep; flipud(sweep)], [y_sw + 1.96*sd_sw; flipud(y_sw - 1.96*sd_sw)], ...
                  [0.85 0.9 1], 'LineStyle', 'none'); hold on;
    h_line = plot(sweep, y_sw, 'b-', 'LineWidth', 2.0);
    h_crit = yline(1.0, 'r--', 'LineWidth', 1.5);
    
    % Definition of labels and fixed angle per subplot
    if v == 1
        ang_name  = '$\Phi$';
        ang_unit  = '$\Phi$ (deg)';
        fixed_str = sprintf('$\\Psi = %.0f^\\circ$', base_design(2));
    else
        ang_name  = '$\Psi$';
        ang_unit  = '$\Psi$ (deg)';
        fixed_str = sprintf('$\\Phi = %.0f^\\circ$', base_design(1));
    end
    
    xlim([lb(v), ub(v)]);
    title(sprintf('Effect of %s on $\\lambda_b$ (fixed %s)', ang_name, fixed_str), ...
          'Interpreter', 'latex', 'FontSize', 10);
    xlabel(ang_unit, 'Interpreter', 'latex'); 
    ylabel('Buckling Load Factor $\lambda_b$', 'Interpreter', 'latex');
    legend([h_line, h_fill, h_crit], {'GPR Mean', '95\% CI', 'Limit $\lambda_b = 1.0$'}, ...
           'Location', 'best', 'Interpreter', 'latex');
    grid on;
end

% Overall title specifying fixed base ply repetitions
sgtitle(sprintf('Main Effects of Angles on $\\lambda_b$ ($r_{\\textrm{Root}}=%d, r_{\\textrm{Mid}}=%d, r_{\\textrm{Tip}}=%d$)', ...
        base_design(3), base_design(4), base_design(5)), ...
        'Interpreter', 'latex', 'FontSize', 12);

%% FIGURE 3B: MAIN EFFECTS OF REPETITIONS (GRADUATED X-AXIS)
figure('Name', 'Fig 3B: Main Effects Repetitions (Graduated X-Axis)', 'Position', [150, 180, 950, 420]);
rep_vars = [3, 4, 5];
rep_names = {'$r_{\textrm{Root}}$', '$r_{\textrm{Mid}}$', '$r_{\textrm{Tip}}$'};

for v = 1:3
    subplot(1, 3, v);
    var_idx = rep_vars(v);
    sweep = linspace(lb(var_idx), ub(var_idx), 60)';
    X_sw = repmat(base_design, 60, 1); 
    X_sw(:, var_idx) = sweep;
    [y_sw, sd_sw] = predict(gpr_blf, X_sw);
    
    h_fill = fill([sweep; flipud(sweep)], [y_sw + 1.96*sd_sw; flipud(y_sw - 1.96*sd_sw)], ...
                  [0.85 0.9 1], 'LineStyle', 'none'); hold on;
    h_line = plot(sweep, y_sw, 'b-', 'LineWidth', 2.0);
    h_crit = yline(1.0, 'r--', 'LineWidth', 1.5);
    
    xlim([lb(var_idx), ub(var_idx)]);
    set(gca, 'XTick', lb(var_idx):2:ub(var_idx));
    
    y_min_val = min(y_sw - 1.96*sd_sw);
    y_max_val = max(y_sw + 1.96*sd_sw);
    yticks_vec = floor(y_min_val*10)/10 : 0.1 : ceil(y_max_val*10)/10;
    set(gca, 'YTick', yticks_vec);
    
    title(sprintf('Effect of %s on $\\lambda_b$', rep_names{v}), ...
          'Interpreter', 'latex', 'FontSize', 10);
    xlabel(rep_names{v}, 'Interpreter', 'latex'); 
    ylabel('Buckling Load Factor $\lambda_b$', 'Interpreter', 'latex');
    legend([h_line, h_fill, h_crit], {'GPR Mean', '95\% CI', 'Limit $\lambda_b = 1.0$'}, ...
           'Location', 'best', 'Interpreter', 'latex');
    grid on;
end

% Overall title specifying fixed reference angles
sgtitle(sprintf('Main Effects of Repetitions on $\\lambda_b$ (Fixed Angles: $\\Phi = %.1f^\\circ, \\Psi = %.1f^\\circ$)', ...
        base_design(1), base_design(2)), ...
        'Interpreter', 'latex', 'FontSize', 12);

%% FIGURE 4A: TRADEOFF BUCKLING VS MASS
figure('Name', 'Fig 4A: Buckling vs Mass colored by Strain Energy', 'Position', [200, 200, 850, 480]);
h_inf  = scatter(mc_mass(~stable_idx), mc_blf(~stable_idx), 15, [0.85 0.85 0.85], 'filled'); hold on;
h_stab = scatter(stable_mass, stable_blf, 20, stable_U, 'filled');
cb = colorbar; cb.Label.Interpreter = 'latex'; cb.Label.String = 'Strain Energy $U$ (J)';
h_fea  = plot(y_mass, y_blf, 'k^', 'MarkerSize', 6, 'MarkerFaceColor', 'y', 'LineWidth', 1.0);
h_opt  = plot(stable_mass(idx_min_m), stable_blf(idx_min_m), 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
h_lim  = yline(1.0, 'r--', 'LineWidth', 2.0);
title('Buckling Load Factor $\lambda_b$ vs. Total Mass');
xlabel('Total Mass (kg)'); ylabel('Buckling Load Factor $\lambda_b$');
legend([h_inf, h_stab, h_fea, h_opt, h_lim], ...
       {'Infeasible Domain ($\lambda_b < 1.0$)', 'Feasible MC ($\lambda_b \ge 1.0$)', ...
        'FEA Training Samples (40 LHS)', 'Optimal Min Mass Design', 'Stability Limit $\lambda_b = 1.0$'}, ...
       'Location', 'northwest'); grid on;

%% FIGURE 4B: TRADEOFF BUCKLING VS STRAIN ENERGY U
figure('Name', 'Fig 4B: Buckling vs Strain Energy colored by Mass', 'Position', [220, 220, 850, 480]);
h_inf2  = scatter(mc_strain(~stable_idx), mc_blf(~stable_idx), 15, [0.85 0.85 0.85], 'filled'); hold on;
h_stab2 = scatter(stable_U, stable_blf, 20, stable_mass, 'filled');
cb = colorbar; cb.Label.Interpreter = 'latex'; cb.Label.String = 'Total Mass (kg)';
h_fea2  = plot(y_strain, y_blf, 'k^', 'MarkerSize', 6, 'MarkerFaceColor', 'y', 'LineWidth', 1.0);
h_opt2  = plot(min_strain_val, stable_blf(idx_min_u), 's', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k');
h_lim2  = yline(1.0, 'r--', 'LineWidth', 2.0);
title('Buckling Load Factor $\lambda_b$ vs. Strain Energy $U$');
xlabel('Strain Energy $U$ (J)'); ylabel('Buckling Load Factor $\lambda_b$');
legend([h_inf2, h_stab2, h_fea2, h_opt2, h_lim2], ...
       {'Infeasible Domain ($\lambda_b < 1.0$)', 'Feasible MC ($\lambda_b \ge 1.0$)', ...
        'FEA Training Samples (40 LHS)', 'Optimal Min Strain Energy Design', 'Stability Limit $\lambda_b = 1.0$'}, ...
       'Location', 'northwest'); grid on;

%% FIGURE 5: KRIGING UNCERTAINTY MAP
figure('Name', 'Fig 5: Uncertainty Map', 'Position', [240, 240, 650, 480]);
X_g2D_nom = [grid_phi(:), grid_psi(:), repmat(base_design(3:5), numel(grid_phi), 1)];
[~, Z_SD] = predict(gpr_blf, X_g2D_nom);
Z_SD = reshape(Z_SD, size(grid_phi));
contourf(grid_phi, grid_psi, Z_SD, 20, 'LineColor', 'none'); hold on;
cb = colorbar; cb.Label.Interpreter = 'latex'; cb.Label.String = 'GPR Std Dev $\sigma_{\lambda_b}$';
h_pts = plot(X_p2(:,1), X_p2(:,2), 'r+', 'MarkerSize', 10, 'LineWidth', 1.8);
title('Buckling Uncertainty $\sigma(\Phi, \Psi)$');
xlabel('$\Phi$ (deg)'); ylabel('$\Psi$ (deg)');
legend(h_pts, 'LHS FEA Samples', 'Location', 'northeast'); grid on;

%% FIGURE 6: ZONE MEAN & TAPERING GRADIENT EFFECT
figure('Name', 'Fig 6: Zone Tapering & Tapering Gradient Effect', 'Position', [260, 260, 980, 400]);
subplot(1,2,1);
feas_reps   = X_mc(stable_idx, 3:5);
infeas_reps = X_mc(~stable_idx, 3:5);
mean_feas   = mean(feas_reps, 1);
mean_infeas = mean(infeas_reps, 1);
b_bar = bar([mean_feas; mean_infeas]', 'grouped');
b_bar(1).FaceColor = [0.2 0.7 0.3];
b_bar(2).FaceColor = [0.8 0.3 0.3];
set(gca, 'XTickLabel', {'$r_{\textrm{Root}}$', '$r_{\textrm{Mid}}$', '$r_{\textrm{Tip}}$'});
title('Mean Ply Repetitions: Feasible vs. Infeasible');
ylabel('Number of Repetitions ($r_i$)');
legend({'Feasible ($\lambda_b \ge 1.0$)', 'Infeasible ($\lambda_b < 1.0$)'}, 'Location', 'northeast'); grid on;

subplot(1,2,2);
delta_r = X_p2(:, 3) - X_p2(:, 5);
scatter(delta_r, y_blf, 55, y_mass, 'filled', 'MarkerEdgeColor', 'k');
cb = colorbar; cb.Label.Interpreter = 'latex'; cb.Label.String = 'Total Mass (kg)';
yline(1.0, 'r--', 'LineWidth', 1.5);
title('Buckling Load Factor vs. Tapering');
xlabel('Tapering $\Delta r = (r_{\textrm{Root}} - r_{\textrm{Tip}})$');
ylabel('Buckling Load Factor $\lambda_b$'); grid on;

%% FIGURE 6B: SPANWISE THICKNESS PROFILES
figure('Name', 'Fig 6B: Spanwise Thickness Profiles', 'Position', [280, 280, 850, 450]);
span_positions = [0.0, 0.5, 1.0];
prof_opt_m = stable_X(idx_min_m, 3:5);
prof_opt_u = stable_X(idx_min_u, 3:5);
prof_max_b = stable_X(idx_max_b, 3:5);
prof_mean  = mean(stable_X(:, 3:5), 1);

plot(span_positions, prof_opt_m, 'g-s', 'LineWidth', 2.2, 'MarkerFaceColor', 'g', 'MarkerSize', 8); hold on;
plot(span_positions, prof_opt_u, 'm-d', 'LineWidth', 2.2, 'MarkerFaceColor', 'm', 'MarkerSize', 8);
plot(span_positions, prof_max_b, 'c-^', 'LineWidth', 2.2, 'MarkerFaceColor', 'c', 'MarkerSize', 8);
plot(span_positions, prof_mean, 'k--o', 'LineWidth', 1.8, 'MarkerFaceColor', 'k', 'MarkerSize', 6);

for k = 1:3
    text(span_positions(k), prof_opt_m(k) + 0.45, sprintf('%.1f', prof_opt_m(k)), ...
         'Color', [0 0.5 0], 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_positions(k), prof_opt_u(k) - 0.45, sprintf('%.1f', prof_opt_u(k)), ...
         'Color', [0.6 0 0.6], 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_positions(k), prof_max_b(k) + 0.45, sprintf('%.1f', prof_max_b(k)), ...
         'Color', [0 0.5 0.7], 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    text(span_positions(k), prof_mean(k) - 0.45, sprintf('%.1f', prof_mean(k)), ...
         'Color', 'k', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
end
set(gca, 'XTick', span_positions, 'XTickLabel', {'Root ($x/L=0$)', 'Mid ($x/L=0.5$)', 'Tip ($x/L=1.0$)'});
xlim([-0.12, 1.12]); ylim([8, 32]);
title('Spanwise Laminate Thickness Profiles ($r_i$) for Optimal Designs');
ylabel('Ply Repetition Count ($r_i$)');
xlabel('Normalized Spanwise Location ($x/L$)');
legend({'Optimal Min Mass', 'Optimal Min Strain Energy', 'Max Stability Design', 'Feasible Mean Profile'}, ...
       'Location', 'northeast'); grid on;

%% FIGURE 6C: MULTI-PANEL CONTOUR MAPS (r_Root VS r_Tip)
figure('Name', 'Fig 6C: Multi-panel Contours r_Root vs r_Tip', 'Position', [150, 200, 1250, 420]);
r_mid_vals = [16, 20, 24];

for m_idx = 1:3
    r_mid_curr = r_mid_vals(m_idx);
    subplot(1, 3, m_idx);
    
    % Mesh adapted to the physical domain of each subplot
    r_root_vec = r_mid_curr:0.2:30;
    r_tip_vec  = 10:0.2:r_mid_curr;
    [grid_r_root, grid_r_tip] = meshgrid(r_root_vec, r_tip_vec);
    
    X_g_6c = [repmat(35, numel(grid_r_root), 1), repmat(60, numel(grid_r_root), 1), ...
              grid_r_root(:), repmat(r_mid_curr, numel(grid_r_root), 1), grid_r_tip(:)];
          
    Z_blf_6c = reshape(predict(gpr_blf, X_g_6c), size(grid_r_root));
    
    contourf(grid_r_root, grid_r_tip, Z_blf_6c, 20, 'LineColor', 'none'); hold on;
    cb = colorbar; cb.Label.Interpreter = 'latex'; cb.Label.String = '$\lambda_b$';
    
    [C_6c, h_6c] = contour(grid_r_root, grid_r_tip, Z_blf_6c, [1.0 1.0], 'r-', 'LineWidth', 2.5);
    
    title(sprintf('$r_{\\textrm{Mid}} = %d$', r_mid_curr));
    xlabel('Root Repetitions $r_{\textrm{Root}}$'); 
    ylabel('Tip Repetitions $r_{\textrm{Tip}}$');
    
    % Adjust limits to exact boundaries of physical domain
    xlim([r_mid_curr, 30]); 
    ylim([10, r_mid_curr]); 
    grid on;
    
    if ~isempty(C_6c) && size(C_6c,2) > 0
        legend(h_6c, 'Limit $\lambda_b = 1.0$', 'Location', 'northwest');
    end
end
sgtitle('Buckling Factor Contours ($r_{\textrm{Root}}$ vs $r_{\textrm{Tip}}$)', 'Interpreter', 'latex', 'FontSize', 12);

%% FIGURE 6D: MULTI-PANEL SPECIFIC STRUCTURAL EFFICIENCY MAPS
figure('Name', 'Fig 6D: Multi-panel Specific Efficiency Maps', 'Position', [100, 150, 1350, 420]);
r_tip_vals = [10, 13, 16];

% 1. Evaluation grid and global precomputation of values
[grid_r_root_d, grid_r_mid_d] = meshgrid(15:0.2:30, 10:0.2:28);
Z_eff_all  = cell(1,3);
Z_blf_all  = cell(1,3);
Z_mass_all = cell(1,3);

cmin = Inf; cmax = -Inf;

for t_idx = 1:3
    r_tip_curr = r_tip_vals(t_idx);
    
    % Input matrix for GPR predictions
    X_g_6d = [repmat(35, numel(grid_r_root_d), 1), repmat(60, numel(grid_r_root_d), 1), ...
              grid_r_root_d(:), grid_r_mid_d(:), repmat(r_tip_curr, numel(grid_r_root_d), 1)];

    Z_blf_mat  = reshape(predict(gpr_blf, X_g_6d), size(grid_r_root_d));
    Z_mass_mat = reshape(predict(gpr_mass, X_g_6d), size(grid_r_root_d));
    Z_eff_mat  = Z_blf_mat ./ Z_mass_mat;

    % Mask unphysical geometries (r_Root < r_Mid or r_Mid < r_Tip)
    invalid_domain = (grid_r_root_d < grid_r_mid_d) | (grid_r_mid_d < r_tip_curr);
    Z_eff_mat(invalid_domain) = NaN;

    Z_eff_all{t_idx}  = Z_eff_mat;
    Z_blf_all{t_idx}  = Z_blf_mat;
    Z_mass_all{t_idx} = Z_mass_mat;

    % Global color range (excluding NaNs)
    valid_vals = Z_eff_mat(~isnan(Z_eff_mat));
    if ~isempty(valid_vals)
        cmin = min(cmin, min(valid_vals));
        cmax = max(cmax, max(valid_vals));
    end
end

% 2. Generation of 3 subplots
for t_idx = 1:3
    r_tip_curr  = r_tip_vals(t_idx);
    x_min_panel = max(15, r_tip_curr); % Dynamic lower bound to prevent blank space
    
    subplot(1, 3, t_idx);
    
    Z_eff_mat  = Z_eff_all{t_idx};
    Z_blf_mat  = Z_blf_all{t_idx};
    Z_mass_mat = Z_mass_all{t_idx};
    
    leg_handles = [];
    leg_labels  = {};

    % Neutral background for non-physical regions
    set(gca, 'Color', [0.95 0.95 0.95]);

    % Efficiency contours with unified scale
    levels = linspace(cmin, cmax, 25);
    contourf(grid_r_root_d, grid_r_mid_d, Z_eff_mat, levels, 'LineColor', 'none'); hold on;
    clim([cmin, cmax]); % Identical color scale across all 3 panels
    
    cb = colorbar; 
    cb.Label.Interpreter = 'latex'; 
    cb.Label.String = 'Specific Efficiency ($\lambda_b / \mathrm{Mass}$) [kg$^{-1}$]';

    % A. Critical Buckling boundary (\lambda_b = 1.0)
    [C_b6d, h_b6d] = contour(grid_r_root_d, grid_r_mid_d, Z_blf_mat, [1.0 1.0], 'r-', 'LineWidth', 2.2);
    if ~isempty(C_b6d) && size(C_b6d, 2) > 0
        leg_handles = [leg_handles; h_b6d];
        leg_labels{end+1} = 'Limit $\lambda_b = 1.0$';
    end

    % B. Optimum stable point (\lambda_b >= 1.0) and informative text box
    Z_eff_stable = Z_eff_mat;
    Z_eff_stable(Z_blf_mat < 1.0) = NaN; 
    
    if any(~isnan(Z_eff_stable(:)))
        [~, max_idx] = max(Z_eff_stable(:));
        x_opt = grid_r_root_d(max_idx);
        y_opt = grid_r_mid_d(max_idx);
        m_opt = Z_mass_mat(max_idx);
        l_opt = Z_blf_mat(max_idx);
        
        h_peak = plot(x_opt, y_opt, 'kp', ...
                      'MarkerFaceColor', 'y', 'MarkerSize', 12, 'LineWidth', 1.2, 'Clipping', 'off');
        leg_handles = [leg_handles; h_peak];
        leg_labels{end+1} = 'Optimum Stable Point';
        
        % Floating label for Mass and Lambda
        txt_opt = sprintf('\\textbf{Optimum:}\n$m = %.2f$ kg\n$\\lambda_b = %.2f$', m_opt, l_opt);
        text(x_opt - 0.5, y_opt, txt_opt, ...
             'Interpreter', 'latex', 'FontSize', 8, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
             'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.3 0.3 0.3], 'Margin', 3);
    end

    % C. Physical geometric limit (r_Root = r_Mid) adapted to active range
    h_bound = plot([x_min_panel, 28], [x_min_panel, 28], 'k--', 'LineWidth', 1.5);
    leg_handles = [leg_handles; h_bound];
    leg_labels{end+1} = 'Limit ($r_{\textrm{Root}}=r_{\textrm{Mid}}$)';

    % Formatting, limits and grid
    title(sprintf('$r_{\\textrm{Tip}} = %d$', r_tip_curr));
    xlabel('Root Repetitions $r_{\textrm{Root}}$'); 
    ylabel('Mid Repetitions $r_{\textrm{Mid}}$');
    xlim([x_min_panel, 30]); 
    ylim([r_tip_curr, 28]); 
    grid on;

    % Dynamic legend
    if ~isempty(leg_handles)
        legend(leg_handles, leg_labels, 'Location', 'northwest', 'FontSize', 8);
    end
end

sgtitle('Specific Efficiency Maps ($\lambda_b / \mathrm{Mass}$) Across Tip Zones ($r_{\textrm{Tip}}$)', ...
        'Interpreter', 'latex', 'FontSize', 12);

%% FIGURE 7: INDIVIDUAL SENSITIVITY BY ZONE
figure('Name', 'Fig 7: Zone Repetitions Sensitivity', 'Position', [340, 340, 950, 400]);

subplot(1,3,1);
scatter(X_p2(:,3), y_blf, 45, y_mass, 'filled', 'MarkerEdgeColor', 'k'); hold on;
yline(1.0, 'r--', 'LineWidth', 1.5);
title('Root Zone ($r_{\textrm{Root}}$)'); xlabel('$r_{\textrm{Root}}$'); ylabel('Buckling Factor $\lambda_b$'); grid on;

subplot(1,3,2);
scatter(X_p2(:,4), y_blf, 45, y_mass, 'filled', 'MarkerEdgeColor', 'k'); hold on;
yline(1.0, 'r--', 'LineWidth', 1.5);
title('Mid Zone ($r_{\textrm{Mid}}$)'); xlabel('$r_{\textrm{Mid}}$'); ylabel('Buckling Factor $\lambda_b$'); grid on;

subplot(1,3,3);
scatter(X_p2(:,5), y_blf, 45, y_mass, 'filled', 'MarkerEdgeColor', 'k'); hold on;
cb = colorbar; cb.Label.Interpreter = 'latex'; cb.Label.String = 'Mass (kg)';
yline(1.0, 'r--', 'LineWidth', 1.5);
title('Tip Zone ($r_{\textrm{Tip}}$)'); xlabel('$r_{\textrm{Tip}}$'); ylabel('Buckling Factor $\lambda_b$'); grid on;

sgtitle('Zone Repetitions Sensitivity ($\lambda_b$ vs. $r_i$)', 'Interpreter', 'latex', 'FontSize', 13);

%% FIGURE 8: PERCENTILE THRESHOLDS FOR STABILITY
figure('Name', 'Fig 8: Minimum Ply Thresholds for Stability', 'Position', [360, 360, 850, 460]);
min_total_plies = min(sum(stable_X(:,3:5), 2));

b = bar([p0_r_root, p5_r_root, p50_r_root; ...
         p0_r_mid,  p5_r_mid,  p50_r_mid; ...
         p0_r_tip,  p5_r_tip,  p50_r_tip], 'grouped');
b(1).FaceColor = [0.8 0.3 0.3]; 
b(2).FaceColor = [0.2 0.6 0.4]; 
b(3).FaceColor = [0.2 0.4 0.8];

set(gca, 'XTickLabel', {'$r_{\textrm{Root}}$', '$r_{\textrm{Mid}}$', '$r_{\textrm{Tip}}$'});
title('Required Repetition Thresholds for Structural Stability ($\lambda_b \ge 1.0$)');
ylabel('Repetitions Count');
legend({'Absolute Minimum', '5th Percentile (Conservative Safe Boundary)', '50th Percentile (Median Feasible)'}, ...
       'Location', 'northwest'); grid on;

txt_stats = {
    '\textbf{Predictive Stability Thresholds:}';
    sprintf('$\\bullet$ Min $r_{\\textrm{Root}}$: %.1f repeats (P5 Safe: $\\ge$ %.1f | P50: %.1f)', p0_r_root, p5_r_root, p50_r_root);
    sprintf('$\\bullet$ Min $r_{\\textrm{Mid}}$ : %.1f repeats (P5 Safe: $\\ge$ %.1f | P50: %.1f)', p0_r_mid, p5_r_mid, p50_r_mid);
    sprintf('$\\bullet$ Min $r_{\\textrm{Tip}}$ : %.1f repeats (P5 Safe: $\\ge$ %.1f | P50: %.1f)', p0_r_tip, p5_r_tip, p50_r_tip);
    sprintf('$\\bullet$ Min Cumulative Repetition Sum: %d repeats', round(min_total_plies))
};
annotation('textbox', [0.44, 0.48, 0.45, 0.35], 'String', txt_stats, ...
           'Interpreter', 'latex', 'FitBoxToText', 'on', 'BackgroundColor', [0.95 0.95 0.95 0.85]);

%% 4.1. MULTI-ALGORITHM OPTIMIZATION (MONTE CARLO vs. FMINCON vs. GA)

% 1. Objective function definition (Minimize Mass)
obj_fun = @(x) predict(gpr_mass, x);

% 2. Geometric linear constraints: r_Root >= r_Mid >= r_Tip (A*x <= b)
% Matrix A: [-r_Root + r_Mid <= 0; -r_Mid + r_Tip <= 0]
A = [ 0, 0, -1,  1,  0;
      0, 0,  0, -1,  1 ];
b = [0; 0];

% 3. Nonlinear stability constraint: lambda_b >= 1.0 --> (1.0 - lambda_b <= 0)
nonlcon = @(x) deal(1.0 - predict(gpr_blf, x), []); 

% -------------------------------------------------------------------------
% A. ALGORITHM 1: MONTE CARLO (GLOBAL EXPLORATORY SEARCH)
% -------------------------------------------------------------------------
tic;
x_opt_mc    = stable_X(idx_min_m, :);
mass_opt_mc = min_mass_val;
blf_opt_mc  = stable_blf(idx_min_m);
t_mc = toc;

% -------------------------------------------------------------------------
% B. ALGORITHM 2: FMINCON (GRADIENT-BASED REFINEMENT / SQP)
% -------------------------------------------------------------------------
% Use Monte Carlo optimum as starting point (x0)
x0 = x_opt_mc;

opts_fmc = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'off', ...
    'TolFun', 1e-6, ...
    'TolX', 1e-6);

tic;
[x_opt_fmc, mass_opt_fmc, exitflag_fmc] = fmincon(obj_fun, x0, A, b, [], [], lb, ub, nonlcon, opts_fmc);
t_fmc = toc;
blf_opt_fmc = predict(gpr_blf, x_opt_fmc);

% -------------------------------------------------------------------------
% C. ALGORITHM 3: GENETIC ALGORITHM (GA - STOCHASTIC HEURISTIC SEARCH)
% -------------------------------------------------------------------------
opts_ga = optimoptions('ga', ...
    'PopulationSize', 120, ...
    'MaxGenerations', 80, ...
    'Display', 'off', ...
    'UseParallel', false);

tic;
[x_opt_ga, mass_opt_ga, exitflag_ga] = ga(obj_fun, 5, A, b, [], [], lb, ub, nonlcon, opts_ga);
t_ga = toc;
blf_opt_ga = predict(gpr_blf, x_opt_ga);

% -------------------------------------------------------------------------
% D. COMPARATIVE RESULTS TABLE
% -------------------------------------------------------------------------
Method      = {'Monte Carlo (MC)'; 'fmincon (SQP)'; 'Genetic Alg. (GA)'};
Mass_kg     = [mass_opt_mc; mass_opt_fmc; mass_opt_ga];
Buckling_LF = [blf_opt_mc; blf_opt_fmc; blf_opt_ga];
Phi_deg     = [x_opt_mc(1); x_opt_fmc(1); x_opt_ga(1)];
Psi_deg     = [x_opt_mc(2); x_opt_fmc(2); x_opt_ga(2)];
r_Root      = [x_opt_mc(3); x_opt_fmc(3); x_opt_ga(3)];
r_Mid       = [x_opt_mc(4); x_opt_fmc(4); x_opt_ga(4)];
r_Tip       = [x_opt_mc(5); x_opt_fmc(5); x_opt_ga(5)];
Time_s      = [t_mc; t_fmc; t_ga];

T_comp = table(Method, Mass_kg, Buckling_LF, Phi_deg, Psi_deg, r_Root, r_Mid, r_Tip, Time_s);

fprintf('\n========================================================================================\n');
fprintf('                           OPTIMIZATION ALGORITHMS COMPARISON                           \n');
fprintf('========================================================================================\n');
disp(T_comp);
fprintf('========================================================================================\n');

%% 6. ACTIVE LEARNING - UNCERTAINTY SAMPLING
% =========================================================================
% 1. Use the Monte Carlo candidate pool (already filtered for feasibility)
X_candidates = X_mc; 

% 2. Evaluate population standard deviation uncertainty for buckling
[~, sigma_blf] = predict(gpr_blf, X_candidates);

% 3. Sort candidates from highest to lowest uncertainty
[sigma_sorted, idx_sort] = sort(sigma_blf, 'descend');

% 4. Select N new samples to evaluate
N_new = 5; 
X_new = X_candidates(idx_sort(1:N_new), :);

% 5. Format 1: Detailed table with FEA (rounded) and continuous values
num_current_samples = size(X_p2, 1);
fprintf('\n================================================================================================================\n');
fprintf('                          SUGGESTED NEW SAMPLES (MAXIMUM UNCERTAINTY - BUCKLING)\n');
fprintf('================================================================================================================\n');
fprintf('Point  Phi(deg)   Psi(deg)   r_Root(FEA)  r_Mid(FEA)   r_Tip(FEA)   r_Root(cont)   r_Mid(cont)    r_Tip(cont)   \n');
fprintf('----------------------------------------------------------------------------------------------------------------\n');
for i = 1:N_new
    p_num = num_current_samples + i;
    phi = X_new(i,1);
    psi = X_new(i,2);
    r_r_cont = X_new(i,3);
    r_m_cont = X_new(i,4);
    r_t_cont = X_new(i,5);
    
    % Rounding for OptiStruct input
    r_r_fea = round(r_r_cont);
    r_m_fea = round(r_m_cont);
    r_t_fea = round(r_t_cont);
    
    % Print formatted row aligned with headers
    fprintf('%-6d %-10.1f %-10.1f %-12d %-12d %-12d %-14.2f %-14.2f %-14.2f\n', ...
        p_num, phi, psi, r_r_fea, r_m_fea, r_t_fea, r_r_cont, r_m_cont, r_t_cont);
end
fprintf('----------------------------------------------------------------------------------------------------------------\n');

% 6. Format 2: Direct output ready to copy and paste into matrix X_p2
fprintf('\n---> Copy and paste these lines at the end of your X_p2 matrix:\n');
for i = 1:N_new
    fprintf('    %5.1f, %5.1f, %5.2f, %5.2f, %5.2f;\n', ...
        X_new(i,1), X_new(i,2), X_new(i,3), X_new(i,4), X_new(i,5));
end
fprintf('\n');

%% GLOBAL 5D UNCERTAINTY METRICS (EVALUATION ON MONTE CARLO POOL)
% =========================================================================
% 1. Extract standard deviation across the valid design space
[~, sigma_global_blf] = predict(gpr_blf, X_mc);

% 2. Calculate global statistics
sigma_mean = mean(sigma_global_blf);
sigma_max  = max(sigma_global_blf);

% 3. Output results to console
fprintf('\n=================================================================\n');
fprintf('       GLOBAL 5D UNCERTAINTY METRICS (BUCKLING)\n');
fprintf('=================================================================\n');
fprintf(' Current training samples           : %d\n', size(X_p2, 1));
fprintf(' Mean Uncertainty (Mean Sigma)      : %8.5f\n', sigma_mean);
fprintf(' Max Uncertainty (Max Sigma)        : %8.5f\n', sigma_max);
fprintf('=================================================================\n');

%% 5D UNCERTAINTY VISUALIZATION: COMPARATIVE HISTOGRAM
% =========================================================================
history_file = 'sigma_history.mat';

figure('Name', '5D Uncertainty Evolution', 'Color', 'w');
hold on; box on; grid on;

% 1. Load and plot previous state if history file exists (40 samples)
if isfile(history_file)
    prev_data = load(history_file);
    histogram(prev_data.sigma_blf_40, 'BinWidth', 0.005, 'FaceColor', [0.8 0.2 0.2], ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Initial DoE (40 Samples)');
end

% 2. Plot current state
current_state = sprintf('Current State (%d Samples)', size(X_p2, 1));
histogram(sigma_global_blf, 'BinWidth', 0.005, 'FaceColor', [0.2 0.6 0.8], ...
    'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', current_state);

% 3. Save initial state if evaluating at baseline stage (40 samples)
if size(X_p2, 1) == 40
    sigma_blf_40 = sigma_global_blf;
    save(history_file, 'sigma_blf_40');
end

% 4. Figure formatting
xlabel('Predicted Uncertainty \sigma (Buckling Factor)', 'Interpreter', 'tex');
ylabel('Frequency (Number of 5D Configurations)', 'Interpreter', 'tex');
title('Global Uncertainty Reduction', 'Interpreter', 'tex');
legend('Location', 'northeast');
set(gca, 'FontSize', 11);
hold off;
