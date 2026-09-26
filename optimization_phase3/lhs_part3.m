% =========================================================================
% Data Generation - Phase 3: Focused 10D Asymmetric LHS Design (56 + 4)
% Multi-Start Optimization, Advanced Diagnostics & Quality Assessment
% =========================================================================
close all;
clearvars;
clc;
rng(42, 'twister');

%% 1. SAMPLING BUDGET, DESIGN BOUNDS & ANCHORS
n_asym_samples = 56;    % New asymmetric FEA cases for Phase 3
n_anchors      = 4;     % Verified Phase 2 symmetric anchors (Cases 51-54)
n_oversample   = 5000;  % Candidate pool size for filter evaluation
n_candidates   = 100;   % Multi-start LHS candidate evaluations

% Variable Order (10D):
% [1:phi_up, 2:psi_up, 3:r_Root_up, 4:r_Mid_up, 5:r_Tip_up, ...
%  6:phi_lo, 7:psi_lo, 8:r_Root_lo, 9:r_Mid_lo, 10:r_Tip_lo]

lb_angles_up = [65,  0];   ub_angles_up = [90, 20]; 
lb_angles_lo = [10, 65];   ub_angles_lo = [40, 90]; 

lb_r_up = [24, 18, 13];    ub_r_up = [30, 28, 25];  
lb_r_lo = [ 8,  6,  4];    ub_r_lo = [16, 13, 10];  

lb = [lb_angles_up, lb_r_up, lb_angles_lo, lb_r_lo];
ub = [ub_angles_up, ub_r_up, ub_angles_lo, ub_r_lo];

% Verified Phase 2 Symmetric Anchors (4 cases)
anchors_P2 = [
    28.1, 73.1, 25, 19, 14,  28.1, 73.1, 25, 19, 14; % Min-mass anchor
    75.5,  4.5, 30, 28, 22,  75.5,  4.5, 30, 28, 28; % Min-energy anchor
    88.7,  2.1, 30, 26, 25,  88.7,  2.1, 30, 26, 25; % Max-buckling anchor
    18.4, 76.4, 30, 22, 19,  18.4, 76.4, 30, 22, 19  % Knee point
];

%% 2. MULTI-START LHS SEARCH (OPTIMIZE R_MAX & SPACE-FILLING)
best_R_max = Inf;
best_X_fea_56 = [];
best_X_scaled_56 = [];
best_yield = 0;

fprintf('Executing multi-start LHS optimization over %d iterations...\n', n_candidates);

for k = 1:n_candidates
    % Generate raw Latin Hypercube Sample in [0, 1] normalized space
    X_norm_raw = lhsdesign(n_oversample, 10);
    X_scaled_raw = lb + X_norm_raw .* (ub - lb);

    % Discretize 'r' (repetitions) to integer values BEFORE applying structural filters
    X_fea_raw = X_scaled_raw;
    X_fea_raw(:, [3:5, 8:10]) = round(X_scaled_raw(:, [3:5, 8:10]));

    % Structural Monotonicity Filters (Eqs. 6.23 - 6.26)
    mon_up = (X_fea_raw(:,3) >= X_fea_raw(:,4)) & (X_fea_raw(:,4) >= X_fea_raw(:,5));
    mon_lo = (X_fea_raw(:,8) >= X_fea_raw(:,9)) & (X_fea_raw(:,9) >= X_fea_raw(:,10));

    % Local Minimum Repetition Thresholds
    sum_root = X_fea_raw(:, 3) + X_fea_raw(:, 8); 
    sum_mid  = X_fea_raw(:, 4) + X_fea_raw(:, 9); 
    sum_tip  = X_fea_raw(:, 5) + X_fea_raw(:, 10);
    local_mask = (sum_root >= 40) & (sum_mid >= 30) & (sum_tip >= 22);

    % Aggregate Repetition Filter
    r_sigma_3_raw = sum_root + sum_mid + sum_tip;
    aggregate_mask = (r_sigma_3_raw >= 94);

    % Upper Skin Allocation Ratio Filter
    r_sum_up_raw = X_fea_raw(:, 3) + X_fea_raw(:, 4) + X_fea_raw(:, 5);
    allocation_ratio_raw = r_sum_up_raw ./ r_sigma_3_raw;
    allocation_mask = (allocation_ratio_raw >= 0.60);

    % Combined Validation Mask
    valid_mask = mon_up & mon_lo & local_mask & aggregate_mask & allocation_mask;
    valid_indices = find(valid_mask);

    if length(valid_indices) < n_asym_samples
        continue;
    end

    % Extract candidate sample subset
    sel_k = valid_indices(1:n_asym_samples);
    X_scaled_cand = X_scaled_raw(sel_k, :);
    X_fea_cand = X_fea_raw(sel_k, :);

    % Correlation and distance metrics calculation
    R_k = corrcoef(X_scaled_cand);
    R_off_k = abs(R_k(triu(true(10), 1)));
    R_max_k = max(R_off_k);

    X_norm_cand = (X_scaled_cand - lb) ./ (ub - lb);
    dist_k = pdist(X_norm_cand);
    d_min_k = min(dist_k);

    % Selection criterion: Minimum R_max with space-filling constraint d_min >= 0.40
    if d_min_k >= 0.40 && R_max_k < best_R_max
        best_R_max = R_max_k;
        best_X_fea_56 = X_fea_cand;
        best_X_scaled_56 = X_scaled_cand;
        best_yield = (length(valid_indices) / n_oversample) * 100;
    end
end

if isempty(best_X_fea_56)
    error('Could not find a valid LHS candidate satisfying constraints. Increase n_candidates or lower d_min threshold.');
end

X_scaled_56 = best_X_scaled_56;
X_fea_56    = best_X_fea_56;
X_fea_60    = [anchors_P2; X_fea_56];

%% 3. STATISTICAL EVALUATION & SAMPLING METRICS
X_normalized_56 = (X_scaled_56 - lb) ./ (ub - lb);

% A. Space-Filling Metrics (10D Normalized Domain)
dist_10d  = pdist(X_normalized_56);
d_min     = min(dist_10d);
d_mean    = mean(dist_10d);
d_max     = max(dist_10d);
phi_ratio = d_min / d_mean;

% B. Correlation Metrics
R_56       = corrcoef(X_scaled_56);
R_off_diag = abs(R_56(triu(true(size(R_56)), 1)));
max_corr   = max(R_off_diag);
rms_corr   = sqrt(mean(R_off_diag.^2));

% Variance Inflation Factor (VIF)
vif_vec = diag(inv(R_56));
max_vif = max(vif_vec);

% Centered L2 Discrepancy (CL2)
[N, s] = size(X_normalized_56);
term1 = (13/12)^s;
term2 = (2/N) * sum(prod(1 + 0.5*abs(X_normalized_56 - 0.5) - 0.5*abs(X_normalized_56 - 0.5).^2, 2));
dist_matrix_cl2 = zeros(N, N);
for i = 1:N
    for j = 1:N
        dist_matrix_cl2(i,j) = prod(1 + 0.5*abs(X_normalized_56(i,:) - 0.5) + ...
            0.5*abs(X_normalized_56(j,:) - 0.5) - 0.5*abs(X_normalized_56(i,:) - X_normalized_56(j,:)));
    end
end
term3 = (1/(N^2)) * sum(dist_matrix_cl2(:));
CL2   = sqrt(term1 - term2 + term3);

% C. Discretization Error (Continuous vs Discrete Repetitions)
r_diff    = abs(X_scaled_56(:, [3:5, 8:10]) - X_fea_56(:, [3:5, 8:10]));
mae_r     = mean(r_diff(:));
max_err_r = max(r_diff(:));

% D. Identification of the Closest Point Pair (10D d_min)
D_mat = squareform(dist_10d);
D_mat(1:size(D_mat,1)+1:end) = Inf; % Ignore main diagonal (self-distance)
[~, min_linear_idx] = min(D_mat(:));
[pt1_idx, pt2_idx]  = ind2sub(size(D_mat), min_linear_idx);

% Output Metrics to Console
disp('===================================================================');
disp('          OPTIMIZED LHS SAMPLING QUALITY EVALUATION (PHASE 3)      ');
disp('===================================================================');
fprintf('Filter Yield Rate                            : %.2f%%\n', best_yield);
fprintf('Space-Filling Metrics (10D Normalized Domain):\n');
fprintf('  - Minimum Distance (d_min)                 : %.4f\n', d_min);
fprintf('  - Mean Distance (d_mean)                   : %.4f\n', d_mean);
fprintf('  - Uniformity Ratio (d_min / d_mean)        : %.4f\n', phi_ratio);
fprintf('  - Centered L2 Discrepancy (CL2)            : %.4f\n', CL2);
fprintf('Correlation Metrics (10D Matrix):\n');
fprintf('  - Maximum Off-Diagonal Correlation (|R|_max) : %.4f\n', max_corr);
fprintf('  - RMS Correlation                          : %.4f\n', rms_corr);
fprintf('  - Maximum Variance Inflation Factor (VIF)  : %.4f\n', max_vif);
fprintf('Discretization Error (Repetition Rounding):\n');
fprintf('  - Mean Absolute Error (MAE)                : %.4f plies\n', mae_r);
fprintf('  - Maximum Rounding Error                   : %.4f plies\n', max_err_r);
disp('===================================================================');

%% 4. DIAGNOSTIC VISUALIZATIONS
font_size  = 18;
title_size = 24;
label_size = 18;

var_names = {'$\phi_{\mathrm{up}}$', '$\psi_{\mathrm{up}}$', '$r_{\mathrm{R,up}}$', '$r_{\mathrm{M,up}}$', '$r_{\mathrm{T,up}}$', ...
             '$\phi_{\mathrm{lo}}$', '$\psi_{\mathrm{lo}}$', '$r_{\mathrm{R,lo}}$', '$r_{\mathrm{M,lo}}$', '$r_{\mathrm{T,lo}}$'};

% --- Figure 1: Correlation Matrix Heatmap ---
figure('Color', 'w', 'Name', 'Correlation Matrix', 'Position', [100 100 700 600]);
red_danger  = [0.85, 0.15, 0.15];
white_ideal = [1.00, 1.00, 1.00];

n_half = 128;
map_neg = [linspace(red_danger(1), white_ideal(1), n_half)', ...
           linspace(red_danger(2), white_ideal(2), n_half)', ...
           linspace(red_danger(3), white_ideal(3), n_half)'];
map_pos = [linspace(white_ideal(1), red_danger(1), n_half)', ...
           linspace(white_ideal(2), red_danger(2), n_half)', ...
           linspace(white_ideal(3), red_danger(3), n_half)'];
symmetric_red_map = [map_neg; map_pos];

imagesc(R_56); colormap(symmetric_red_map); clim([-1 1]); 
cb = colorbar;
ylabel(cb, 'Correlation Coefficient ($R$)', 'Interpreter', 'latex', 'FontSize', font_size);
cb.FontSize = font_size;

xticks(1:10); yticks(1:10); 
xticklabels(var_names); yticklabels(var_names);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
title('10D Correlation Matrix', 'Interpreter', 'latex', 'FontSize', title_size);
axis square; grid off;

% --- Figure 2: 10D Euclidean Distance Distribution ---
figure('Color', 'w', 'Name', 'Distance Density', 'Position', [150 150 700 550]);
histogram(dist_10d, 20, 'FaceColor', [0.2 0.5 0.7], 'EdgeColor', 'k', 'Normalization', 'pdf'); hold on;
xline(d_min, 'r--', sprintf('d_{min} = %.3f', d_min), 'LineWidth', 2, 'Interpreter', 'tex', 'FontSize', font_size);
xline(d_mean, 'k--', sprintf('d_{mean} = %.3f', d_mean), 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'Interpreter', 'tex', 'FontSize', font_size);
title('10D Normalized Euclidean Distance Distribution', 'Interpreter', 'latex', 'FontSize', title_size);
xlabel('10D Normalized Distance', 'Interpreter', 'latex', 'FontSize', label_size); 
ylabel('Probability Density', 'Interpreter', 'latex', 'FontSize', label_size);
set(gca, 'FontSize', font_size);
grid on; box on; axis square;

% --- Figure 3: Upper Skin Angle Subspace ---
figure('Color', 'w', 'Name', 'Upper Angle Subspace', 'Position', [200 200 750 550]);

% 1. LHS Points (56 cases)
h_lhs = plot(X_scaled_56(:,1), X_scaled_56(:,2), 'b*', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;

% 2. Highlight closest point pair in 10D space
h_pair = plot(X_scaled_56([pt1_idx, pt2_idx], 1), X_scaled_56([pt1_idx, pt2_idx], 2), 'ro', ...
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(X_scaled_56([pt1_idx, pt2_idx], 1), X_scaled_56([pt1_idx, pt2_idx], 2), 'r--', 'LineWidth', 1.5);

title('Upper Skin Angle Subspace', 'Interpreter', 'latex', 'FontSize', title_size);
xlabel('$\phi_{\mathrm{up}}$ (deg)', 'Interpreter', 'latex', 'FontSize', label_size); 
ylabel('$\psi_{\mathrm{up}}$ (deg)', 'Interpreter', 'latex', 'FontSize', label_size);
xlim([65 90]); ylim([0 20]); grid on; grid minor; box on;
set(gca, 'FontSize', font_size);
legend([h_lhs, h_pair], {'LHS Points', sprintf('Closest Pair ($d_{\\mathrm{min}} = %.4f$)', d_min)}, ...
    'Interpreter', 'latex', 'Location', 'bestoutside', 'FontSize', font_size-4);

% --- Figure 4: Lower Skin Angle Subspace ---
figure('Color', 'w', 'Name', 'Lower Angle Subspace', 'Position', [250 250 750 550]);

% 1. LHS Points (56 cases) - Columns 6 (phi_lo) & 7 (psi_lo)
h_lhs = plot(X_scaled_56(:,6), X_scaled_56(:,7), 'k*', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;

% 2. Highlight closest point pair in 10D space
h_pair = plot(X_scaled_56([pt1_idx, pt2_idx], 6), X_scaled_56([pt1_idx, pt2_idx], 7), 'ro', ...
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none');
plot(X_scaled_56([pt1_idx, pt2_idx], 6), X_scaled_56([pt1_idx, pt2_idx], 7), 'r--', 'LineWidth', 1.5);

title('Lower Skin Angle Subspace', 'Interpreter', 'latex', 'FontSize', title_size);
xlabel('$\phi_{\mathrm{lo}}$ (deg)', 'Interpreter', 'latex', 'FontSize', label_size); 
ylabel('$\psi_{\mathrm{lo}}$ (deg)', 'Interpreter', 'latex', 'FontSize', label_size);
xlim([10 40]); ylim([65 90]); grid on; grid minor; box on;
set(gca, 'FontSize', font_size);
legend([h_lhs, h_pair], {'LHS Points', sprintf('Closest Pair ($d_{\\mathrm{min}} = %.4f$)', d_min)}, ...
    'Interpreter', 'latex', 'Location', 'bestoutside', 'FontSize', font_size-4);

% --- Figure 5: Aggregate Repetition Filter Compliance ---
figure('Color', 'w', 'Name', 'Aggregate Filter', 'Position', [300 300 700 550]);
r_sig_3_selected = sum(X_fea_56(:, [3:5, 8:10]), 2);
histogram(r_sig_3_selected, 'BinWidth', 2, 'FaceColor', [0.3 0.7 0.4], 'EdgeColor', 'k'); hold on;
xline(94, 'r-', 'Threshold (r_{\Sigma,3} \geq 94)', 'LineWidth', 2.5, 'LabelOrientation', 'aligned', 'Interpreter', 'tex', 'FontSize', font_size);
title('Aggregate Repetition Filter ($r_{\Sigma,3}$)', 'Interpreter', 'latex', 'FontSize', title_size);
xlabel('Total Repetitions $r_{\Sigma,3}$', 'Interpreter', 'latex', 'FontSize', label_size); 
ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', label_size);
set(gca, 'FontSize', font_size);
grid on; box on; axis square;

% --- Figure 6: Upper Skin Allocation Ratio Filter Compliance ---
figure('Color', 'w', 'Name', 'Allocation Ratio Filter', 'Position', [350 350 700 550]);
alloc_ratio_selected = sum(X_fea_56(:, 3:5), 2) ./ r_sig_3_selected;
histogram(alloc_ratio_selected * 100, 'BinWidth', 1.5, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'k'); hold on;
xline(60, 'r-', 'Threshold (\geq 60%)', 'LineWidth', 2.5, 'LabelOrientation', 'aligned', 'Interpreter', 'tex', 'FontSize', font_size-3);
title('Upper Skin Allocation Ratio', 'Interpreter', 'latex', 'FontSize', title_size);
xlabel('Upper Skin Repeats Proportion (\%)', 'Interpreter', 'latex', 'FontSize', label_size); 
ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', label_size);
set(gca, 'FontSize', font_size);
grid on; box on; axis square;

% --- Figure 7: Discrete Repetition Profiles (Upper & Lower Skins) ---
figure('Color', 'w', 'Name', 'Repetition Profiles', 'Position', [400 400 750 550]);
zones = [1, 2, 3];

% Calculate mean values per zone
mean_up = mean(X_fea_56(:, 3:5), 1);
mean_lo = mean(X_fea_56(:, 8:10), 1);

plot(zones, X_fea_56(:, 3:5)', 'Color', [0.8 0.3 0.3 0.1], 'LineWidth', 1.0); hold on;
plot(zones, X_fea_56(:, 8:10)', 'Color', [0.3 0.3 0.8 0.1], 'LineWidth', 1.0);
h_up = plot(zones, mean_up, 'r-o', 'LineWidth', 3, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
h_lo = plot(zones, mean_lo, 'b-o', 'LineWidth', 3, 'MarkerFaceColor', 'b', 'MarkerSize', 8);

% Add numerical text labels with mean values for each zone
for i = 1:length(zones)
    % Upper Skin Mean (shifted upward)
    text(zones(i), mean_up(i) + 1.2, sprintf('%.2f', mean_up(i)), ...
        'Color', [0.8 0 0], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'FontSize', font_size - 2);
    
    % Lower Skin Mean (shifted downward)
    text(zones(i), mean_lo(i) - 1.2, sprintf('%.2f', mean_lo(i)), ...
        'Color', [0 0 0.8], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'FontSize', font_size - 2);
end

title('Discrete Repetition Profiles', 'Interpreter', 'latex', 'FontSize', title_size);
xlabel('Zone', 'Interpreter', 'latex', 'FontSize', label_size); 
ylabel('Discrete Repetitions', 'Interpreter', 'latex', 'FontSize', label_size);
xticks([1 2 3]); xticklabels({'Root', 'Mid', 'Tip'});
xlim([0.8 3.2]); ylim([2 32]); grid on; box on;
set(gca, 'FontSize', font_size);
legend([h_up, h_lo], {'Mean Upper Skin', 'Mean Lower Skin'}, ...
    'Interpreter', 'latex', 'Location', 'bestoutside', 'FontSize', font_size - 4);

% --- Figure 8: Optimized 10D Parallel Coordinates ---
figure('Color', 'w', 'Name', '10D Parallel Coordinates', 'Position', [450 450 900 500]);
x_axes = 1:10;
hold on;

% Plot individual normalized sample lines with transparency (Alpha = 0.15)
plot(x_axes, X_normalized_56', 'Color', [0.2 0.4 0.8 0.15], 'LineWidth', 1.2);

% Central tendency profiles (Mean and Median)
mean_profile   = mean(X_normalized_56, 1);
median_profile = median(X_normalized_56, 1);

h_mean = plot(x_axes, mean_profile, 'r-o', 'LineWidth', 3, 'MarkerSize', 7, 'MarkerFaceColor', 'r');
h_med  = plot(x_axes, median_profile, 'k--s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k');

xticks(1:10); xticklabels(var_names);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
title('10D Parallel Coordinates', 'Interpreter', 'latex', 'FontSize', title_size);
xlabel('Design Variables', 'Interpreter', 'latex', 'FontSize', label_size);
ylabel('Normalized Values [0, 1]', 'Interpreter', 'latex', 'FontSize', label_size);
xlim([0.8 10.2]); ylim([-0.02 1.02]);
grid on; box on;

legend([h_mean, h_med], {'Mean Profile', 'Median Profile'}, ...
    'Interpreter', 'latex', 'Location', 'bestoutside', 'FontSize', font_size - 2);

%% 5. DATA EXPORT TO TEXT FILE
filename = 'lhs_phase3_60cases2.txt';
fileID = fopen(filename, 'w');

fprintf(fileID, '========================================================================================================================\n');
fprintf(fileID, '  PHASE 3: 10D LHS INITIAL TRAINING SET (60 CASES: 4 P2 ANCHORS + 56 NEW ASYMMETRIC)\n');
fprintf(fileID, '========================================================================================================================\n');
fprintf(fileID, '%-6s %-10s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Pt', 'Type', 'Phi_up', 'Psi_up', 'r_R_up', 'r_M_up', 'r_T_up', 'Phi_lo', 'Psi_lo', 'r_R_lo', 'r_M_lo', 'r_T_lo');
fprintf(fileID, '------------------------------------------------------------------------------------------------------------------------\n');

for i = 1:60
    type_str = 'ASYM';
    if i <= 4, type_str = 'ANCHOR'; end
    fprintf(fileID, '%-6d %-10s %-8.1f %-8.1f %-8d %-8d %-8d %-8.1f %-8.1f %-8d %-8d %-8d\n', ...
        i, type_str, X_fea_60(i,1), X_fea_60(i,2), X_fea_60(i,3), X_fea_60(i,4), X_fea_60(i,5), ...
        X_fea_60(i,6), X_fea_60(i,7), X_fea_60(i,8), X_fea_60(i,9), X_fea_60(i,10));
end

fclose(fileID);
disp(['Data matrix exported successfully to: ', filename]);
