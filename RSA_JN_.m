% In the name of Allah
% Author: Jalaledin Noroozi
% Email: jalaledin.noroozi@gmail.com
% Code for PNAS article titled: "Frontotemporal Network Contribution to Occluded Face Processing"
% Tracking #: 2024-07457R

clc
clear
close all

%% Load PFC RSA Data
load('PFC_RSA.mat')

%% Perform Bootstrap Test for Significance on PFC RSA Data
ii = 0; sig_boot = [];
for i = [2, 5]
    ii = ii + 1;
    for itime = 1:size(mat_pfc_sp, 2)
        % Calculate the proportion of values ≤ 0 (null hypothesis) using bootstrap sampling
        sig_boot(ii, itime) = sum(mat_pfc_sp(:, itime, i) <= 0) / size(mat_pfc_sp, 1);
    end
end

%% Plot Bootstrapped PFC RSA with Confidence Intervals
color_g = [
    0,      0.4470, 0.7410;               
    0.8500, 0.3250, 0.0980;          
    0.9290, 0.6940, 0.1250;          
    0.4940, 0.1840, 0.5560;       
    0.4660, 0.6740, 0.1880;         
    0.3010, 0.7450, 0.9330;       
    0.6350, 0.0780, 0.1840;
    0, 0, 0];

hh = figure;
options.alpha = 0.3;
options.line_width = 1;
options.error = 'std';
options.smooth = 10;
options.x_axis = ti_win;

% Plot PFC RSA for intact and fully occluded conditions
options.color = color_g(2, :);
niceplot_r(mat_pfc_sp(:,:,2), options)
hold on
options.color = color_g(8, :);
niceplot_r(mat_pfc_sp(:,:,5), options)
hold on
axis tight

% Mark significance based on bootstrap test for each timepoint
sig_it = nan(1, size(mat_pfc_sp, 2));
for j = 1:size(mat_pfc_sp, 2)
    if sig_boot(1, j) < 0.05, sig_it(j) = -0.05; end
end
plot(ti_win(ti_win > 0 & ti_win < 500), sig_it(ti_win > 0 & ti_win < 500), '.', 'MarkerSize', 12, 'Color', [0 0.5 0.5]);

% Repeat for the fully occluded condition
sig_it = nan(1, size(mat_pfc_sp, 2));
for j = 1:size(mat_pfc_sp, 2)
    if sig_boot(2, j) < 0.05, sig_it(j) = -0.03; end
end
plot(ti_win(ti_win > 0 & ti_win < 500), sig_it(ti_win > 0 & ti_win < 500), '.', 'MarkerSize', 12, 'Color', [0.3 0.3 0.3]);

ylim([-0.06 0.15])
xlim([-50 500])
ylabel('Correlation')
xlabel('Time (ms)')
title('vlPFC OSNs')
set(gca, 'box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'fontsize', 12);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--')
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--')
line([350 350], ylim, 'Color', 'k', 'LineStyle', '--')

hh = findobj(gca, 'Type', 'line');
legend([hh(length(hh)), hh(length(hh)-1)], 'Intact', 'Full Occlusion', 'FontSize', 12, 'LineWidth', 1.2, 'box', 'off')

%% Bar Plot for Peak RSA in PFC during Specific Time Windows
color_fig = [
    0, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560];

st1 = 150; st2 = 200;
time_indices = find(ti_win > st1 & ti_win < st2);

for cond = 2:5
    IT_bar(cond-1, :) = mean(mat_pfc_sp(:, time_indices, cond), 2);
    it_m(cond-1) = mean(IT_bar(cond-1, :));
    it_e(cond-1) = std(IT_bar(cond-1, :)) / sqrt(size(IT_bar, 2));
end

figure; hold on
for k = 1:4
    bar(k, it_m(k), 'FaceColor', color_fig(k, :), 'EdgeColor', 'none', 'LineWidth', 1)
    errorbar(k, it_m(k), it_e(k), '.', 'Color', [0 0 0], 'LineWidth', 0.5)
end
xlim([0.5 4.5])
xlabel({'% Visible', 'Increasing Occlusion'})
ylabel('Correlation')
title('vlPFC OSNs')
set(gca, 'box', 'off', 'TickDir', 'out', 'LineWidth', 0.5, 'fontsize', 12)

%% Statistical Analysis with ANOVA and Post Hoc Test
[p, tbl, stats] = anova1(IT_bar');
fprintf('ANOVA F-Statistic: %f, p-value: %f\n', tbl{2,6})

% Tukey-Kramer post hoc test for multiple comparisons
[c, m, h, nms] = multcompare(stats, 'CType', 'tukey-kramer');
posthoc_results = table(nms(c(:,1)), nms(c(:,2)), c(:,3), c(:,5), c(:,6), 'VariableNames', {'Group1', 'Group2', 'Difference', 'LowerCI', 'UpperCI'});
disp(posthoc_results);

significant_pairs = c(c(:,6) < 0.05, 1:2);
disp('Pairs with significant differences:')
disp(nms(significant_pairs))

%% Load IT OSN Data and Perform Similar Processing as PFC
load('IT_RSA.mat')

% Baseline correction for IT data
ind_b = find(ti_win > -50 & ti_win < 0);
base_c = squeeze(mean(mean(mat_it_sp(:, ind_b, :), 1), 2));
base_rep = repmat(base_c, 1, size(mat_it_sp, 1), size(mat_it_sp, 2));
mat_it_sn_ = mat_it_sp - permute(base_rep, [2, 3, 1]);

ii = 0; sig_boot = [];
for i = [2, 5]
    ii = ii + 1;
    for itime = 1:size(mat_it_sn_, 2)
        sig_boot(ii, itime) = sum(mat_it_sn_(:, itime, i) <= 0) / size(mat_it_sn_, 1);
    end
end

% Plotting similar bootstrap results for IT RSA data
hh = figure;
options.color = color_g(6, :);
niceplot_r(mat_it_sn_(:,:,2), options)
hold on
options.color = color_g(8, :);
niceplot_r(mat_it_sn_(:,:,5), options)
hold on
axis tight
ylim([-0.09 0.15])
xlim([-50 500])
ylabel('Correlation')
xlabel('Time (ms)')
title('ITC OSNs')
set(gca, 'box', 'off', 'TickDir', 'out', 'LineWidth', 0.5, 'fontsize', 12);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--')
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--')
line([350 350], ylim, 'Color', 'k', 'LineStyle', '--')

lg = legend([hh(length(hh)), hh(length(hh)-1)], 'Intact', 'Full Occlusion');
lg.FontSize = 12;
lg.LineWidth = 1.2;
legend boxoff
