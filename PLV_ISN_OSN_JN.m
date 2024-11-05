% In the name of Allah
% Author: Jalaledin Noroozi
% Email: jalaledin.noroozi@gmail.com
% Code for PNAS article titled: "Frontotemporal Network Contribution to Occluded Face Processing"
% Tracking #: 2024-07457R

%% Load Data
load('E:\Jalal\PFC_IT\Data\Jenab\Main\Code\LFP\Data_plv\PLV_IT_PFC_all_shuf.mat')

%% Define Parameters and Compute PLV Significance
freq = [1:1:45 55:5:160];
ind_b = 1:100;
ind_stim = 401:1200;
tt = ind_stim - 500;

% Frequency band and time windows for baseline and stimulus
freq_b = [7 17];
st_t = 170;
en_t = 200;
base_ac = find(tt > st_t & tt < en_t);
base_ix = find(tt > -(en_t - st_t) & tt < 0);

% Initialize significance arrays
pp = []; p_value = []; hh = [];

for ss = 1:111
    var_b = squeeze(nanmean(plv_it_pfc_all(ss, :, freq >= freq_b(1) & freq <= freq_b(2), base_ix), 3));
    var_h = squeeze(nanmean(plv_it_pfc_all(ss, :, freq >= freq_b(1) & freq <= freq_b(2), base_ac), 3));
    
    for ii = 1:4
        % Perform Wilcoxon signed-rank and t-tests
        [pp(ss, ii), ~] = signrank(var_b(ii, :)', var_h(ii, :)', 'tail', 'left', 'alpha', 0.05);
        [hh(ss, ii), p_value(ss, ii)] = ttest(var_b(ii, :)', var_h(ii, :)', 'Tail', 'left', 'Alpha', 0.05);
    end
end

% Identify significant indices based on rank and t-tests
sig_rank = find(pp(:, 1) <= 0.05 & pp(:, 4) <= 0.05);
sig_ttest = find(p_value(:, 1) <= 0.05 & p_value(:, 4) <= 0.05);

%% Compute Base-Corrected Heatmap
fr_ix = find(freq > 3 & freq < 160);
val_b = squeeze(plv_shufC_all(:, :, fr_ix, base_ix));
val_bm = nanmean(val_b, 4);
val_sd = std(val_bm, [], 3);

plv_norm = (plv_shufC_all(:, :, fr_ix, :) - val_bm) ./ val_sd;

% Apply Gaussian smoothing
sigma_freq = 2;
sigma_time = 15;
plv_norm_smoothed = zeros(size(plv_norm));

for ss = 1:length(sig_rank)
    for condi = 1:4
        plv_norm_smoothed(ss, condi, :, :) = imgaussfilt(squeeze(plv_norm(sig_rank(ss), condi, :, :)), [sigma_freq, sigma_time]);
    end
end

% Average and plot smoothed PLV values
var_hh = squeeze(nanmean(plv_norm_smoothed, 1));
name_str = {'Intact', 'Low', 'High', 'Full'};

for ii = [1, 4]
    figure('Position', [ii * 380, 400, 350, 300])
    ax = subplot(1, 1, 1);
    h = pcolor(tt, freq(fr_ix), squeeze(var_hh(ii, :, :)));
    h.EdgeColor = 'none';
    colormap(jet); colorbar('Label', 'Norm. PLV (a.u.)');
    ax.YDir = 'normal'; ax.YScale = 'log';
    line([0 0], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
    line([350 350], ylim, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
    caxis([-4 4]); % Set color limits
    ax.YTick = [8 18 40 100 160];
    ax.XTick = [0 350 500];
    xlim([-40 500]);
    xlabel('Time (ms)');
    ylabel('Freq. (Hz)');
end

%% Plot PLV Difference Between Conditions Over Time
freq_b = [7 17];
c = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
options.alpha = 0.5;
options.line_width = 1;
options.error = 'sem';
options.smooth = 10;
options.x_axis = tt;

for ifr = 1
    figure('Position', [400 400 350 350])
    ind_freq = (freq >= freq_b(ifr, 1) & freq <= freq_b(ifr, 2));
    var_h1 = squeeze(nanmean(plv_norm(sig_rank, 1, ind_freq, :), 3));
    var_h4 = squeeze(nanmean(plv_norm(sig_rank, 4, ind_freq, :), 3));
    
    % Plotting condition curves
    options.color = c(1, :); niceplot_r(var_h1, options); hold on
    options.color = c(7, :); niceplot_r(var_h4, options);
    
    xlim([-40 500]);
    xticks([0 150 350 500]);
    line([0 0], ylim, 'Linewidth', 0.3, 'Color', 'k', 'LineStyle', '--');
    line([350 350], ylim, 'Linewidth', 0.3, 'Color', 'k', 'LineStyle', '--');
    xlabel('Time (ms)');
    ylabel('Phase Locking Value');
    legend({'Intact', 'Full Oc.'}, 'Location', 'northeast');
end
