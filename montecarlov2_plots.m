%% =========================================================
%% FINAL PLOTS
%% =========================================================

% Font settings
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize', 10); 

% Create the figure 
fig_hist = figure('Name', 'Monte Carlo Overlay', 'Color', 'w', ...
             'Units', 'inches', 'Position', [1, 1, 3.5, 2.8]); 
hold on; grid on; box on;

% 1. Plot Baseline 
h1 = histogram(freq_min_baseline, 'BinWidth', 0.001, ...
    'FaceColor', [0.85 0.32 0.09], 'FaceAlpha', 0.6, 'EdgeColor', 'k');

% 2. Plot Proposed 
h2 = histogram(freq_min_proposed, 'BinWidth', 0.001, ...
    'FaceColor', [0 0.447 0.741], 'FaceAlpha', 0.7, 'EdgeColor', 'k');

% 3. Draw the Safety Limit Wall
f_Safe_Limit = input('Enter the frequency safety limit in p.u.')
x1 = xline(-f_Safe_Limit, 'k--', {'UFLS Safety','Limit'}, 'DisplayName','XLINE','LineWidth', 1.5);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';

% 4. Formatting
xlabel('Minimum Frequency Nadir, \Delta f_{min} (p.u.)');
ylabel('Number of Occurrences');
legend([h1, h2], {'Conventional AGC', 'Proposed HO-CBF'}, ...
    'Location', 'west', 'FontSize', 8);

% 5. Export for IEEE
exportgraphics(fig_hist, 'MonteCarlo_Overlay.tif', 'Resolution', 1200);
disp('Overlaid histogram saved as MonteCarlo_Overlay.tif!');
% % prob_fail_prop = sum(freq_min_proposed <= -0.02) / Nsim * 100;
% % 
% % fprintf('\n--- Monte Carlo Results ---\n');
% % fprintf('Baseline Failure Rate: %.1f%%\n', prob_fail_base);
% % fprintf('Proposed Failure Rate: %.1f%%\n', prob_fail_prop);
