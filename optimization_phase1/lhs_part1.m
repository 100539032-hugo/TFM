% =========================================================================
% LATIN HYPERCUBE SAMPLING (LHS) FOR COMPOSITE LAYUP DESIGN
% Part 1: Dynamic Data Generation and Space-Filling Evaluation
% =========================================================================

close all
clear all
clc

% --- Sampling Setup ---
n_samples = 20;       % Number of sample points
lb = [0, 0];          % Lower bounds [Phi_min, Psi_min] in degrees
ub = [90, 90];        % Upper bounds [Phi_max, Psi_max] in degrees

% Generate sample points using Latin Hypercube Sampling
[X_scaled, X_normalized] = lhsdesign_modified(n_samples, lb, ub);

% -------------------------------------------------------------------------
% 1. MATHEMATICAL EVALUATION
% -------------------------------------------------------------------------

% A. Correlation Matrix (Checks linear dependence between Phi and Psi)
R = corrcoef(X_scaled(:,1), X_scaled(:,2));
disp('--- LHS Design Evaluation ---');
fprintf('Correlation between phi and psi: %.4f\n', R(1,2));

% B. Maximin Criterion (Evaluates minimum Euclidean distance between points)
% pdist calculates pairwise Euclidean distances across the sample space
distances = pdist(X_scaled);
min_distance = min(distances);
fprintf('Minimum distance between two points (Maximin): %.4f\n\n', min_distance);

% Identify the two closest points to highlight in the plot
dist_matrix = squareform(distances);
dist_matrix(dist_matrix == 0) = inf; % Ignore self-distance (diagonal = 0)
[~, linear_idx] = min(dist_matrix(:));
[point1, point2] = ind2sub(size(dist_matrix), linear_idx);

% -------------------------------------------------------------------------
% 2. PLOT 1: Design Space and Maximin Evaluation
% -------------------------------------------------------------------------

figure('Color', 'w')
plot(X_scaled(:,1), X_scaled(:,2), 'b*', 'MarkerSize', 8)
hold on

% Highlight the two closest points detected by the maximin criterion
plot(X_scaled([point1, point2], 1), X_scaled([point1, point2], 2), ...
    'ro', 'MarkerSize', 10, 'LineWidth', 1.5)

% Format title and axis labels using LaTeX rendering
title('Latin Hypercube Sampling - Experiment Part 1', 'Interpreter', 'latex', 'FontSize', 28)
xlabel('$\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('$\psi$ (deg)', 'Interpreter', 'latex', 'FontSize', 24)

% Set axis limits and major tick intervals
xlim([0 90])
ylim([0 90])
xticks(0:5:90)
yticks(0:5:90)

% Draw thin grid lines corresponding to LHS bin size (90 / 20 = 4.5 degrees)
bin_size = (ub(1) - lb(1)) / n_samples;
xline(0:bin_size:90, ':', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
yline(0:bin_size:90, ':', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');

% Axis visual formatting
ax = gca;
ax.FontSize = 20; 
ax.GridAlpha = 0.25;

box on
legend('LHS Points', 'Closest Points', 'Location', 'best', 'FontSize', 22)
hold off

% -------------------------------------------------------------------------
% 3. PLOT 2: Marginal Uniformity (Scatterhist)
% -------------------------------------------------------------------------

figure('Color', 'w')

[hAxes] = scatterhist(X_scaled(:,1), X_scaled(:,2), ...
    'NBins', [20 20], ...           % 20 bins matching LHS interval resolution
    'Location', 'SouthEast', ...
    'Direction', 'out', ...
    'Color', [0 0.4470 0.7410], ... % Accent blue color
    'Marker', '*', ...
    'MarkerSize', 9);

% Global title
sgtitle('LHS Marginal Uniformity', 'Interpreter', 'latex', 'FontSize', 28);

% Limits and ticks strictly set to 10-degree increments
xlim(hAxes(1), [0 90])
ylim(hAxes(1), [0 90])
xticks(hAxes(1), 0:10:90)
yticks(hAxes(1), 0:10:90)

% Axis labels
xlabel(hAxes(1), '$\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 22)
ylabel(hAxes(1), '$\psi$ (deg)', 'Interpreter', 'latex', 'FontSize', 22)

% Grid and font formatting for scatterhist main plot
grid(hAxes(1), 'on')
set(hAxes(1), ...
    'FontSize', 18, ...          % Tick numbers font size
    'GridAlpha', 0.25, ...       % Soft grid transparency
    'LineWidth', 0.5, ...        % Thin, crisp grid lines
    'GridLineStyle', ':', ...    % Fine dotted style
    'GridColor', [0.2 0.2 0.2]); % Neutral gray

% -------------------------------------------------------------------------
% 4. EXPORT ANGLE PAIRS TO TEXT FILE (.TXT)
% -------------------------------------------------------------------------

% Round to 1 decimal place (0.1 deg) based on ply manufacturing tolerances
X_rounded = round(X_scaled, 1);

filename = 'lhs_angle_pairs.txt';
fileID = fopen(filename, 'w');

% File header
fprintf(fileID, '=========================================================\n');
fprintf(fileID, '  LHS ANGLE PAIRS (PLY MANUFACTURING)\n');
fprintf(fileID, '  Tolerance / Precision: 0.1 degrees\n');
fprintf(fileID, '=========================================================\n');
fprintf(fileID, '%-8s %-15s %-15s\n', 'Point', 'Phi (deg)', 'Psi (deg)');
fprintf(fileID, '---------------------------------------------------------\n');

% Write each coordinate pair
for i = 1:size(X_rounded, 1)
    fprintf(fileID, '%-8d %-15.1f %-15.1f\n', i, X_rounded(i,1), X_rounded(i,2));
end

fclose(fileID);
disp(['Points successfully saved to: ', filename]);



