%% =========================================================
%% FINAL IEEE OVERLAID HISTOGRAM 
%% =========================================================

% IEEE standard font settings
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize', 10); 

% Create the figure (Sized for single-column IEEE width)
fig_hist = figure('Name', 'Monte Carlo Overlay', 'Color', 'w', ...
             'Units', 'inches', 'Position', [1, 1, 3.5, 2.8]); 
hold on; grid on; box on;

% 1. Plot Baseline (Red, slightly transparent)
h1 = histogram(freq_min_baseline, 'BinWidth', 0.001, ...
    'FaceColor', [0.85 0.32 0.09], 'FaceAlpha', 0.6, 'EdgeColor', 'k');

% 2. Plot Proposed (Blue, slightly transparent)
h2 = histogram(freq_min_proposed, 'BinWidth', 0.001, ...
    'FaceColor', [0 0.447 0.741], 'FaceAlpha', 0.7, 'EdgeColor', 'k');

% 3. Draw the Safety Limit Wall
x1 = xline(-0.008, 'k--', {'UFLS Safety','Limit'}, 'DisplayName','XLINE','LineWidth', 1.5);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';

% 4. Formatting
xlabel('Minimum Frequency Nadir, \Delta f_{min} (p.u.)');
ylabel('Number of Occurrences');
xlim([-0.015 0]); % Zooms in on the exact region of interest
legend([h1, h2], {'Conventional AGC', 'Proposed HO-CBF'}, ...
    'Location', 'west', 'FontSize', 8);

% 5. Export for IEEE
exportgraphics(fig_hist, 'MonteCarlo_Overlay.tif', 'Resolution', 600);
disp('Overlaid histogram saved as MonteCarlo_Overlay.tif!');

% %% Uncertainity analysis plot
% fig10 = figure('Name', 'Uncertainity', 'Color', 'w', ...
%              'Units', 'inches', 'Position', [1, 1, 6.16, 3]); % Adjust height as needed
% hold on; grid on; box on;
%     % Cyber Attack Window Patch
%     %fill([120 121.5 121.5 120], [-0.04 -0.04 0.01 0.01], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% % Plot the 3 Cases
%     histogram(freq_min_baseline, 15, 'FaceColor', 'r', 'FaceAlpha', 0.7); hold on; grid on;
%     xline(-0.02, 'k--', 'Safety Limit', 'LineWidth', 2);
%     title('Baseline PI-AGC (Case 2)'); xlabel('Min \Delta\omega (pu)'); ylabel('Occurrences');
%     xlim([-0.04 0]);
% 
% % Option A: Export as a Vector PDF (Highly Recommended for IEEE)
% % This keeps lines infinitely sharp and texts searchable.
% %exportgraphics(fig6, 'PBatt.pdf', 'ContentType', 'vector');
% 
% % Option B: Export as a Vector EPS (Some IEEE journals prefer this over PDF)
% %exportgraphics(fig6, 'PBatt.eps', 'ContentType', 'vector');
% 
% % Option C: Export as High-Resolution TIFF/PNG (If raster is strictly required)
% % IEEE requires 600 DPI for line-art raster images.
% exportgraphics(fig10, 'Baseline_PI_AGC.tif', 'Resolution', 600);
% 
% %
% fig11 = figure('Name', 'Uncertainity2', 'Color', 'w', ...
%              'Units', 'inches', 'Position', [1, 1, 6.16, 3]); % Adjust height as needed
% hold on; grid on; box on;
%     % Cyber Attack Window Patch
%     %fill([120 121.5 121.5 120], [-0.04 -0.04 0.01 0.01], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% % Plot the 3 Cases
%     histogram(freq_min_proposed, 15, 'FaceColor', 'b', 'FaceAlpha', 0.7); hold on; grid on;
%     xline(-0.02, 'k--', 'Safety Limit', 'LineWidth', 2);
%     title('Proposed CLF-HOCBF (Case 3)'); xlabel('Min \Delta\omega (pu)');
%     xlim([-0.04 0]);
% 
% % Option A: Export as a Vector PDF (Highly Recommended for IEEE)
% % This keeps lines infinitely sharp and texts searchable.
% %exportgraphics(fig6, 'PBatt.pdf', 'ContentType', 'vector');
% 
% % Option B: Export as a Vector EPS (Some IEEE journals prefer this over PDF)
% %exportgraphics(fig6, 'PBatt.eps', 'ContentType', 'vector');
% 
% % Option C: Export as High-Resolution TIFF/PNG (If raster is strictly required)
% % IEEE requires 600 DPI for line-art raster images.
% exportgraphics(fig11, 'Proposed.tif', 'Resolution', 600);
% 
% 
% % %% Statistical Plotting
% % figure('Position', [200, 200, 800, 400], 'Color', 'w');
% % subplot(1,2,1);
% % histogram(freq_min_baseline, 15, 'FaceColor', 'r', 'FaceAlpha', 0.7); hold on; grid on;
% % xline(-0.02, 'k--', 'Safety Limit', 'LineWidth', 2);
% % title('Baseline PI-AGC (Case 2)'); xlabel('Min \Delta\omega (pu)'); ylabel('Occurrences');
% % xlim([-0.04 0]);
% % 
% % subplot(1,2,2);
% % histogram(freq_min_proposed, 15, 'FaceColor', 'b', 'FaceAlpha', 0.7); hold on; grid on;
% % xline(-0.02, 'k--', 'Safety Limit', 'LineWidth', 2);
% % title('Proposed CLF-HOCBF (Case 3)'); xlabel('Min \Delta\omega (pu)');
% % xlim([-0.04 0]);
% % 
% % sgtitle('Monte-Carlo Cyber-Resilience Evaluation (100 Random DoS Profiles)', 'FontWeight', 'bold', 'FontSize', 14);
% % 
% % % Print metrics to console
% % prob_fail_base = sum(freq_min_baseline <= -0.02) / Nsim * 100;
% % prob_fail_prop = sum(freq_min_proposed <= -0.02) / Nsim * 100;
% % 
% % fprintf('\n--- Monte Carlo Results ---\n');
% % fprintf('Baseline Failure Rate: %.1f%%\n', prob_fail_base);
% % fprintf('Proposed Failure Rate: %.1f%%\n', prob_fail_prop);