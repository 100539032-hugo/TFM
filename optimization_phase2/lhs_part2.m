% =========================================================================
% Data Generation - Part 2: 3-Zone Structural Tapering (5D LHS)
% =========================================================================
close all
clear all
clc

% -------------------------------------------------------------------------
% 1. 5D DESIGN SPACE CONFIGURATION
% -------------------------------------------------------------------------
n_samples = 40; 

% Bounds: [phi (deg), psi (deg), r_min, r_min, r_min]
lb_angles = [0, 0];
ub_angles = [90, 90];
r_min = 10; % Minimum repetitions (Tip)
r_max = 30; % Maximum repetitions (Root)

lb = [lb_angles, r_min, r_min, r_min];
ub = [ub_angles, r_max, r_max, r_max];

% Generate 5D LHS: [phi, psi, r_raw1, r_raw2, r_raw3]
[X_scaled, X_normalized] = lhsdesign_modified(n_samples, lb, ub);

% Enforce Tapering Rule: r_Root >= r_Mid >= r_Tip
r_raw = X_scaled(:, 3:5);
r_sorted = sort(r_raw, 2, 'descend'); % Sort each row descending
X_scaled(:, 3:5) = r_sorted;

% Discretize 'r' for FEA simulation (rounded to integer plies)
X_fea = X_scaled;
X_fea(:, 3:5) = round(X_scaled(:, 3:5));

% -------------------------------------------------------------------------
% 2. MATHEMATICAL & QUALITY EVALUATION (5D & 2D)
% -------------------------------------------------------------------------
disp('===========================================================');
disp('      5D LHS SAMPLING QUALITY EVALUATION                   ');
disp('===========================================================');

% A. Correlation Matrix (Global 5D)
R = corrcoef(X_scaled);
disp('Correlation Matrix R (Off-diagonal max correlation):');
max_off_diag_corr = max(abs(R(triu(true(size(R)), 1))));
fprintf('  - Max off-diagonal correlation            : %.4f\n', max_off_diag_corr);
fprintf('  - Correlation phi - psi                   : %.4f\n', R(1,2));
fprintf('  - Correlation r_Root - r_Mid              : %.4f\n', R(3,4));
fprintf('  - Correlation r_Mid - r_Tip               : %.4f\n\n', R(4,5));

% B. Minimum Distance Criterion in 5D Normalized Space
dist_5d = pdist(X_normalized);
min_dist_5d = min(dist_5d);
mean_dist_5d = mean(dist_5d);
uniformity_ratio_5d = min_dist_5d / mean_dist_5d;

fprintf('Minimum Distance Criterion (5D Normalized Space):\n');
fprintf('  - Minimum 5D distance                    : %.4f\n', min_dist_5d);
fprintf('  - Mean inter-point 5D distance           : %.4f\n', mean_dist_5d);
fprintf('  - Space-Filling Ratio (d_min / d_mean)   : %.4f\n\n', uniformity_ratio_5d);

% C. Minimum Distance Criterion in 2D Angles Only (Phi vs Psi)
distances_angles = pdist(X_scaled(:, 1:2));
min_distance_angles = min(distances_angles);

% Identify the two closest points in angle space
dist_matrix_angles = squareform(distances_angles);
dist_matrix_angles(dist_matrix_angles == 0) = inf; % Ignore self-distance
[~, linear_idx] = min(dist_matrix_angles(:));
[point1, point2] = ind2sub(size(dist_matrix_angles), linear_idx);

fprintf('Minimum Distance Criterion in Angles (phi - psi):\n');
fprintf('  - Minimum distance between two points    : %.4f deg\n', min_distance_angles);
fprintf('  - Closest points identified              : Point %d and Point %d\n\n', point1, point2);

% D. Discretization Error (Continuous vs FEA Integer Plies)
r_diff = abs(X_scaled(:, 3:5) - X_fea(:, 3:5));
mae_r = mean(r_diff(:));
max_err_r = max(r_diff(:));

fprintf('Impact of FEA Discretization (Continuous vs Discrete r):\n');
fprintf('  - Mean Absolute Error (MAE)              : %.4f plies\n', mae_r);
fprintf('  - Maximum Rounding Error                 : %.4f plies\n\n', max_err_r);

% E. Tapering Step Statistics (Ply Drop-offs across zones)
drop_root_mid = X_fea(:, 3) - X_fea(:, 4);
drop_mid_tip  = X_fea(:, 4) - X_fea(:, 5);
flat_zones    = sum(drop_root_mid == 0 | drop_mid_tip == 0);

fprintf('Structural Tapering Statistics (Ply Drop-offs):\n');
fprintf('  - Avg. drop Root -> Mid                  : %.2f plies\n', mean(drop_root_mid));
fprintf('  - Avg. drop Mid -> Tip                   : %.2f plies\n', mean(drop_mid_tip));
fprintf('  - Max total drop (Root -> Tip)           : %d plies\n', max(X_fea(:,3) - X_fea(:,5)));
fprintf('  - Untapered sub-zone samples (r_i = r_j) : %d / %d (%.1f%%)\n', ...
    flat_zones, n_samples, (flat_zones/n_samples)*100);
disp('===========================================================');

% -------------------------------------------------------------------------
% 3. FIGURE 1: DESIGN SPACE & MINIMUM DISTANCE EVALUATION
% -------------------------------------------------------------------------
figure('Color', 'w', 'Name', 'Angle Space Evaluation')
plot(X_scaled(:,1), X_scaled(:,2), 'b*', 'MarkerSize', 8, 'LineWidth', 1.2)
hold on

% Highlight the two closest points detected in angle space
plot(X_scaled([point1, point2], 1), X_scaled([point1, point2], 2), ...
    'ro', 'MarkerSize', 10, 'LineWidth', 1.5)

title('Latin Hypercube Sampling - Experiment Part 2', ...
    'Interpreter', 'latex', 'FontSize', 16)
xlabel('$\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\psi$ (deg)', 'Interpreter', 'latex', 'FontSize', 14)

xlim([0 90]); ylim([0 90]);
xticks(0:10:90); yticks(0:10:90);

bin_size = (ub_angles(1) - lb_angles(1)) / n_samples;
xline(0:bin_size:90, ':', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
yline(0:bin_size:90, ':', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');

grid on; axis square; box on;
legend('LHS Points', 'Closest Points', 'Location', 'northeast')

% -------------------------------------------------------------------------
% 4. FIGURE 2: MARGINAL UNIFORMITY OF ANGLES (SCATTERHIST)
% -------------------------------------------------------------------------
figure('Color', 'w', 'Name', 'Angle Scatterhist')
[hAxes] = scatterhist(X_scaled(:,1), X_scaled(:,2), ...
    'NBins', [20 20], ...
    'Location', 'SouthEast', ...
    'Direction', 'out', ...
    'Color', [0 0.4470 0.7410], ...
    'Marker', '*', ...
    'MarkerSize', 8);

sgtitle('Angle Marginal Uniformity ($\phi - \psi$)', 'Interpreter', 'latex', 'FontSize', 16);
xlim(hAxes(1), [0 90]); ylim(hAxes(1), [0 90]);
xticks(hAxes(1), 0:10:90); yticks(hAxes(1), 0:10:90);
xlabel(hAxes(1), '$\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 13)
ylabel(hAxes(1), '$\psi$ (deg)', 'Interpreter', 'latex', 'FontSize', 13)
grid(hAxes(1), 'on')

% -------------------------------------------------------------------------
% 5. FIGURE 3: TAPERING PROFILES AND ZONE DISTRIBUTIONS
% -------------------------------------------------------------------------
figure('Color', 'w', 'Name', 'Tapering Profiles')

% Subplot A: Individual profiles and mean profile with markers
subplot(1,2,1)
zones = [1, 2, 3];
h_samples = plot(zones, X_fea(:, 3:5)', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hold on
mean_vals = mean(X_fea(:, 3:5), 1);
h_mean = plot(zones, mean_vals, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

for z = 1:3
    text(zones(z), mean_vals(z) + 0.9, sprintf('%.2f', mean_vals(z)), ...
        'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 10);
end

title('Repetition Profiles ($r$)', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Structural Zone', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Discrete Repetitions ($r$)', 'Interpreter', 'latex', 'FontSize', 12)
xticks([1 2 3])
xticklabels({'Root (Z1)', 'Mid (Z2)', 'Tip (Z3)'})
ylim([r_min-1, r_max+2.5])
grid on; axis square;
legend([h_samples(1), h_mean], {'LHS Samples', 'Average Profile'}, 'Location', 'northeast')

% Subplot B: Thickness distribution per zone with complete legend
subplot(1,2,2)
boxplot(X_fea(:, 3:5), 'Labels', {'Root (Z1)', 'Mid (Z2)', 'Tip (Z3)'})
title('Repetition Distribution by Zone', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Structural Zone', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Discrete Repetitions ($r$)', 'Interpreter', 'latex', 'FontSize', 12)
grid on; axis square;

hold on;
h_med = plot(NaN, NaN, 'r-', 'LineWidth', 2);
h_box = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);
h_out = plot(NaN, NaN, 'r+', 'MarkerSize', 7, 'LineWidth', 1.5);
legend([h_med, h_box, h_out], ...
    {'Red Line: Median', 'Box: IQR (25th-75th %)', 'Red +: Outliers (>1.5 IQR)'}, ...
    'Location', 'northeast');

sgtitle('Structural Tapering Behavior', 'Interpreter', 'latex', 'FontSize', 16)

% -------------------------------------------------------------------------
% 6. FIGURE 4: 3D TAPERING VOLUME FOR BOTH PHI AND PSI
% -------------------------------------------------------------------------
figure('Color', 'w', 'Name', '3D Tapering Space Dual', 'Position', [100 100 1000 450])

% Subplot A: Colored by Phi
subplot(1,2,1)
scatter3(X_fea(:,3), X_fea(:,4), X_fea(:,5), 50, X_scaled(:,1), 'filled')
colormap(gca, parula);
cb1 = colorbar; cb1.Label.String = '\phi (deg)';
xlabel('r_{Root} (Zone 1)'); ylabel('r_{Mid} (Zone 2)'); zlabel('r_{Tip} (Zone 3)');
title('Colored by Orientation Angle $\phi$', 'Interpreter', 'latex', 'FontSize', 13)
grid on; axis square; view(135, 25);

% Subplot B: Colored by Psi
subplot(1,2,2)
scatter3(X_fea(:,3), X_fea(:,4), X_fea(:,5), 50, X_scaled(:,2), 'filled')
colormap(gca, parula);
cb2 = colorbar; cb2.Label.String = '\psi (deg)';
xlabel('r_{Root} (Zone 1)'); ylabel('r_{Mid} (Zone 2)'); zlabel('r_{Tip} (Zone 3)');
title('Colored by Orientation Angle $\psi$', 'Interpreter', 'latex', 'FontSize', 13)
grid on; axis square; view(135, 25);

sgtitle('3D Tapering Volume ($r_{Root} \ge r_{Mid} \ge r_{Tip}$)', 'Interpreter', 'latex', 'FontSize', 16)

% -------------------------------------------------------------------------
% 7. FIGURE 5: 5D CORRELATION MATRIX HEATMAP (LATEX INTERPRETER)
% -------------------------------------------------------------------------
figure('Color', 'w', 'Name', 'Correlation Heatmap')
var_names = {'$\phi$', '$\psi$', '$r_{\mathrm{Root}}$', '$r_{\mathrm{Mid}}$', '$r_{\mathrm{Tip}}$'};

% Define key colors for symmetric V-shaped colormap (R = 0: White, |R| -> 1: Red)
red_danger  = [0.85, 0.15, 0.15]; % Red for high correlation (|R| -> 1)
white_ideal = [1.00, 1.00, 1.00]; % Pure white for zero correlation (R = 0)

n_half = 128; % Resolution per half-range
map_neg = [linspace(red_danger(1), white_ideal(1), n_half)', ...
           linspace(red_danger(2), white_ideal(2), n_half)', ...
           linspace(red_danger(3), white_ideal(3), n_half)'];
map_pos = [linspace(white_ideal(1), red_danger(1), n_half)', ...
           linspace(white_ideal(2), red_danger(2), n_half)', ...
           linspace(white_ideal(3), red_danger(3), n_half)'];
symmetric_red_map = [map_neg; map_pos];

% Render correlation matrix
imagesc(R);
colormap(symmetric_red_map);
clim([-1 1]);

% Colorbar with LaTeX formatting
cb = colorbar;
cb.Label.String = 'Correlation Coefficient ($R$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

% Configure axes ticks with LaTeX interpreter
xticks(1:5); yticks(1:5);
xticklabels(var_names); yticklabels(var_names);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

title('5D Correlation Matrix Heatmap', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('Variables', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Variables', 'Interpreter', 'latex', 'FontSize', 12);
axis square; box on;

% Contrast-adaptive numerical values overlay with LaTeX support
for i = 1:5
    for j = 1:5
        val = R(i,j);
        if abs(val) > 0.55
            txt_color = 'w'; % White text on dark red background
        else
            txt_color = 'k'; % Black text on light/white background
        end
        text(j, i, sprintf('%.2f', val), 'HorizontalAlignment', 'center', ...
            'Color', txt_color, 'FontWeight', 'bold', 'FontSize', 11, ...
            'Interpreter', 'latex');
    end
end

% -------------------------------------------------------------------------
% 8. EXPORT DATA TO .TXT FILE (OPTISTRUCT / FEA)
% -------------------------------------------------------------------------
X_rounded_angles = round(X_scaled(:, 1:2), 1);
filename = 'lhs_tapered_3zones.txt';
fileID = fopen(filename, 'w');

fprintf(fileID, '================================================================================================================\n');
fprintf(fileID, '  3-ZONE TAPERED LHS POINTS (LAMINATE MANUFACTURING & FEA SIMULATION)\n');
fprintf(fileID, '  Angles shared across zones | Tapering rule: r_Root >= r_Mid >= r_Tip\n');
fprintf(fileID, '================================================================================================================\n');
fprintf(fileID, '%-6s %-10s %-10s %-12s %-12s %-12s %-14s %-14s %-14s\n', ...
    'Point', 'Phi(deg)', 'Psi(deg)', 'r_Root(FEA)', 'r_Mid(FEA)', 'r_Tip(FEA)', ...
    'r_Root(cont)', 'r_Mid(cont)', 'r_Tip(cont)');
fprintf(fileID, '----------------------------------------------------------------------------------------------------------------\n');

for i = 1:n_samples
    fprintf(fileID, '%-6d %-10.1f %-10.1f %-12d %-12d %-12d %-14.2f %-14.2f %-14.2f\n', ...
        i, X_rounded_angles(i,1), X_rounded_angles(i,2), ...
        X_fea(i,3), X_fea(i,4), X_fea(i,5), ...
        X_scaled(i,3), X_scaled(i,4), X_scaled(i,5));
end

fclose(fileID);
disp(['Design points successfully exported to: ', filename]);
