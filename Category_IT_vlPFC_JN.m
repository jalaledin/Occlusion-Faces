%%%%% In the name of Allah
% Author: Jalaledin Noroozi
% Email: jalaledin.noroozi@gmail.com
% Code for PNAS article titled: "Frontotemporal Network Contribution to Occluded Face Processing"
% Tracking #: 2024-07457R

%% Load Data
clc; clear; close all;
load('data.mat');
addpath \Functions

% Load condition-specific data
[Res_z, cond_z] = all_conditions_zebel;
[Res_j, cond_j] = all_conditions_jenab;

% Define face categories: human and monkey
all_f = [1:27, 36:62, 71:97, 106:132, 141:167, 176:310];
f_ix = all_f(1:135); 
m_f = all_f(136:270);
cat_face = [f_ix; m_f]';

%% Calculate Mean Responses for Human and Monkey Faces
win_ = 30;  % Smoothing window for PSTH
psth  = @(x) ndass_smooth(1000 * nanmean(x, 1), win_);
t_step = 2;
win_size = 50;
t_st = 400;
t_end = 1100;
tim_tar1 = t_st:t_step:(t_end - win_size);
tim_tar2 = (t_st + win_size):t_step:t_end;

%% Identify Responsive Neurons in IT Cortex
start_ = 500;
end_ = 900;
stp = 10;
baseline_t = 500 - stp:500;
win_res = [start_:stp:end_ - stp; start_ + stp:stp:end_];
res_it_01 = zeros(size(it, 1), 1);

for ss = 1:size(it, 1)
    base_t = nanmean(it{ss}(:, baseline_t), 1);
    for tim_res = 1:size(win_res, 2)
        win_time = win_res(1, tim_res):win_res(2, tim_res);
        ac_t = nanmean(it{ss}(:, win_time), 1);
        [p, ~] = ranksum(ac_t', base_t');
        if p < 0.001 && mean(ac_t) > mean(base_t) + std(base_t) * 3
            res_it_01(ss) = 1;
            break;
        end
    end
end
it_responsive = find(res_it_01 == 1);

%% Identify Responsive Neurons in PFC Cortex
start_ = 500;
end_ = 850;
stp = 10;
baseline_t = 500 - stp:500;
win_res = [start_:stp:end_ - stp; start_ + stp:stp:end_];
res_pfc_01 = zeros(size(pfc, 1), 1);

for ss = 1:size(pfc, 1)
    base_t = smooth(nanmean(pfc{ss}(:, baseline_t), 1), 50);
    for tim_res = 1:size(win_res, 2)
        win_time = win_res(1, tim_res):win_res(2, tim_res);
        ac_t = smooth(nanmean(pfc{ss}(:, win_time), 1), 50);
        [p, ~] = ranksum(ac_t', base_t');
        if p < 0.01 && mean(ac_t) > mean(base_t) + std(base_t) * 2
            res_pfc_01(ss) = 1;
            break;
        end
    end
end
pfc_responsive = find(res_pfc_01 == 1);

%% Calculate F-statistic for Human vs. Monkey Categories
it_ = it(it_responsive);
pfc_ = pfc(pfc_responsive);
t_res = [];

for t1 = 1:length(tim_tar1)
    tim_ix_1 = tim_tar1(t1):tim_tar2(t1);
    for ne = 1:length(it_)
        IT = []; PFC = [];
        for s = 1:388
            val_h = nanmean(1000 * it_{ne}(s, tim_ix_1), 2);
            IT = [IT; val_h];
            if ne <= size(pfc_, 1)
                PFC = [PFC; nanmean(1000 * pfc_{ne}(s, tim_ix_1), 2)];
            end
        end
        grp_s = ones(388, 1);
        for ii = 1:size(cat_face, 2)
            grp_s(ismember(cond_j, cat_face(:, ii))) = ii + 1;
        end
        ind_keep = ismember(grp_s, 2:ii + 1);

        [~, tb] = anova1(grp_s(ind_keep), IT(ind_keep), 'off');
        t_res.it(t1, ne) = tb{2, 5};

        if ne <= size(pfc_, 1)
            [~, tb] = anova1(grp_s(ind_keep), PFC(ind_keep), 'off');
            t_res.pfc(t1, ne) = tb{2, 5};
        end
    end
end

%% Plot Results for IT and PFC Responses
color_g = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980];
IT_f = t_res.it;
PFC_f = t_res.pfc;
t = (round(tim_tar1 + tim_tar2) / 2) - 500;
hh = figure('Position', [600, 500, 400, 300]);

options.alpha = 0.5;
options.line_width = 2;
options.error = 'sem';
options.smooth = 5;
options.x_axis = t;

% Plot IT responses
options.color = color_g(1, :);
niceplot_r(IT_f', options);
hold on;
text(280, 0.68, strcat('N= ', num2str(size(IT_f, 2))), 'fontsize', 10, 'Color', options.color);

% Plot PFC responses
options.color = color_g(2, :);
niceplot_r(PFC_f', options);
text(280, 0.63, strcat('N= ', num2str(size(PFC_f, 2))), 'fontsize', 10, 'Color', options.color);
axis tight;
title('Human vs. Monkey');
xlabel('Time (ms)');
ylabel('F-statistic');
legend('IT', 'PFC', 'FontSize', 12);
legend boxoff;
