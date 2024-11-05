% Author: Jalaledin Noroozi
% Email: jalaledin.noroozi@gmail.com
% Code for PNAS article titled: "Frontotemporal Network Contribution to Occluded Face Processing"
% Tracking #: 2024-07457R

clc;
clear;

% Load data for neuron responses and slope tuning curves
load('4le_occ_ID_MeanGray_beta_value_overTime_allNeurons_100_1000.mat');
load('Slop_tuningCurve_IT_PFC.mat');

% Define color scheme
c = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250];
     [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330];
     [0.6350, 0.0780, 0.1840]; [0 0 0]];
str_model = {'Occluded', 'Intact', 'Around Zero'};
ti_win = -100:499;

% Set plot options
options.alpha = 0.3;
options.line_width = 1;
options.error = 'sem';
options.smooth = 10;
options.x_axis = ti_win;

% Find indices for occluded/intact neurons in IT and PFC based on coefficients
Occ_pfc = find(coef_pfc > 0);
Int_pfc = find(coef_pfc < 0);
Occ_it = find(coef_it > 0);
Int_it = find(coef_it < 0);
var_ = 4;

% Compute normalized beta values for IT and PFC occluded conditions
bv_it_occ = beta_values(:,:,var_) ./ beta_variance(:,:,var_);
bv_pfc_occ = beta_values_pfc(:,:,var_) ./ beta_variance_pfc(:,:,var_);

% Statistical analysis on peak and trough of responses
[~, lo] = min(smooth(nanmean(bv_it_occ,1), 30));
strOut = statReprt(bv_it_occ(:,223), 3, 'sem');

signrank(nanmean(bv_it_occ(:,80:99), 2), nanmean(bv_it_occ(:,210:229), 2));
[~, lo] = max(smooth(nanmean(bv_pfc_occ,1), 30));
strOut = statReprt(bv_pfc_occ(:,230), 3, 'sem');
signrank(nanmean(bv_pfc_occ(:,80:99), 2), nanmean(bv_pfc_occ(:,220:239), 2));

% Plot IT responses for intact and occluded conditions
figure;
options.color = c(1,:);
niceplot_r(bv_it_occ(Int_it, :), options);
hold on;
options.color = [0 0 0];
niceplot_r(bv_it_occ(Occ_it, :), options);

% Statistical significance markers for IT intact and occluded neurons
sig_it = [];
for j = 1:size(bv_it_occ, 2)
    sig_it(j) = ranksum(bv_it_occ(Int_it, j), nanmean(bv_it_occ(Int_it, 50:100), 2));
    if sig_it(j) < 0.05
        sig_it(j) = 1;
    else
        sig_it(j) = NaN;
    end
end
t_sig = find(ti_win > 0 & ti_win < 400);
ph = plot(ti_win(t_sig), sig_it(t_sig) * -6.4, 'Marker', '.', 'MarkerSize', 9, 'Color', c(1,:));

% Repeat for occluded neurons
sig_it = [];
for j = 1:size(bv_it_occ, 2)
    sig_it(j) = ranksum(bv_it_occ(Occ_it, j), nanmean(bv_it_occ(Occ_it, 50:100), 2));
    if sig_it(j) < 0.05
        sig_it(j) = 1;
    else
        sig_it(j) = NaN;
    end
end
ph = plot(ti_win(t_sig), sig_it(t_sig) * 3.7, 'Marker', '.', 'MarkerSize', 9, 'Color', c(8,:));

% Formatting IT plot
axis tight;
xlim([-100 500]);
ylim([-6.7 3.9]);
set(gca, 'box', 'off', 'TickDir', 'out', 'LineWidth', 0.5);
line([0 0], ylim, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
line(xlim, [0 0], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-');
line([350  350], ylim, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
xlabel('Time (ms)');
ylabel('Normalized Beta');
legend({'ISNs', 'OSNs'}, 'boxoff');
title('IT');

% Plot PFC responses
figure;
options.color = c(7,:);
niceplot_r(bv_pfc_occ(Int_pfc, :), options);
hold on;
options.color = [0 0 0];
niceplot_r(bv_pfc_occ(Occ_pfc, :), options);

% Statistical significance markers for PFC intact and occluded neurons
sig_pfc = [];
for j = 1:size(bv_pfc_occ, 2)
    sig_pfc(j) = ranksum(bv_pfc_occ(Int_pfc, j), nanmean(bv_pfc_occ(Int_pfc, 50:100), 2));
    if sig_pfc(j) < 0.05
        sig_pfc(j) = 1;
    else
        sig_pfc(j) = NaN;
    end
end
ph = plot(ti_win(t_sig), sig_pfc(t_sig) * -1.9, 'Marker', '.', 'MarkerSize', 9, 'Color', c(7,:));

sig_pfc = [];
for j = 1:size(bv_pfc_occ, 2)
    sig_pfc(j) = ranksum(bv_pfc_occ(Occ_pfc, j), nanmean(bv_pfc_occ(Occ_pfc, 50:100), 2));
    if sig_pfc(j) < 0.05
        sig_pfc(j) = 1;
    else
        sig_pfc(j) = NaN;
    end
end
ph = plot(ti_win(t_sig), sig_pfc(t_sig) * 2.2, 'Marker', '.', 'MarkerSize', 9, 'Color', c(8,:));

% Formatting PFC plot
axis tight;
xlim([-100 500]);
ylim([-2 2.4]);
set(gca, 'box', 'off', 'TickDir', 'out', 'LineWidth', 0.5);
line([0 0], ylim, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
line(xlim, [0 0], 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-');
line([350  350], ylim, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
xlabel('Time (ms)');
ylabel('Normalized Beta');
legend({'ISNs', 'OSNs'}, 'boxoff');
title('PFC');
