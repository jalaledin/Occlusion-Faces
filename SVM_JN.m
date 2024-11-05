clc; clear; close all;

% Add path for custom functions
addpath('Functions')

% Load data
load('E:\Jalal\PFC_IT\Data\Jenab\Main\Data2m\IT_PFC_Untrail_trailed.mat');
clearvars -except pfc it trail_pfc trail_it pfc_resp it_resp

% Parameter setup for stimulus conditions
[Res_z, cond_z] = all_conditions_zebel;
[Res_j, cond_j] = all_conditions_jenab;

no_fa_nt = find(Res_j(:,1) < 6 & Res_j(:,4) > 0)';
nt_ = reshape(no_fa_nt, 8, []);
on_ = 311:350;
on_fa_nt = reshape(on_, 8, []);
le_oc = [nt_, on_fa_nt];

% Reordering levels of occlusion
lev_occ = [];
for j = [1, 4, 3, 7, 2, 6, 5, 8]
    lev_occ = [lev_occ; le_oc(j, :)];
end

level_occ_nat = lev_occ';
level_occ_4l = [reshape(lev_occ(1:2, :), [], 1), reshape(lev_occ(3:4, :), [], 1), reshape(lev_occ(5:6, :), [], 1), reshape(lev_occ(7:8, :), [], 1)];
ind_ss = ismember(cond_j', reshape(level_occ_4l, 1, []));

stim_select = reshape(level_occ_4l, [], 1);
W_id = le_oc(:, [1, 2, 3, 9]);
M_id = le_oc(:, [4, 5, 7, 8]);

% Defining intensity levels of occlusion
int_ = lev_occ(1:2, :);
low_ = lev_occ(3:4, :);
hig_ = lev_occ(5:6, :);
full_ = lev_occ(7:8, :);

id_W = reshape([W_id], [], 1);
id_M = reshape([M_id], [], 1);
lab_wm = [id_W, id_M];

% Group occlusion levels for vlPFC and ITC
int_wm = [reshape(int_(:, [1:3, 9]), [], 1), reshape(int_(:, [4, 5, 7, 8]), [], 1)];
low_wm = [reshape(low_(:, [1:3, 9]), [], 1), reshape(low_(:, [4, 5, 7, 8]), [], 1)];
hig_wm = [reshape(hig_(:, [1:3, 9]), [], 1), reshape(hig_(:, [4, 5, 7, 8]), [], 1)];
ful_wm = [reshape(full_(:, [1:3, 9]), [], 1), reshape(full_(:, [4, 5, 7, 8]), [], 1)];

% Define time window parameters
t_step = 5;
win_size = 30;
t_st = 400;
t_end = 1100;
win_h = [t_st:t_step:t_end - win_size; t_st + win_size:t_step:t_end]';
t = nanmean(win_h, 2) - 500;

rate = 0.5;
rep = 1000;
numStim = 90;
boot = 20;
norm_type = 'spectral'; % Options: 'spectral', 'trace', 'frobenius'
out = [];
t_res_shuffle = [];

% Define response indices for vlPFC and ITC
IT_res_4l = [35, 38, 40, 45, 46, 62, 64, 106:119, 121:184];
PFC_res_4l = [24:29, 83:96, 98:121, 123:132, 135:137, 139:159, 161];

% Loop through each occlusion level group
for gi = 1:4
    switch gi
        case 1, ix_label = int_wm;  % Intact occlusion level
        case 2, ix_label = low_wm;  % Low occlusion level
        case 3, ix_label = hig_wm;  % High occlusion level
        case 4, ix_label = ful_wm;  % Full occlusion level
    end

    % Number of trials for IT and PFC
    num_trail_it = arrayfun(@(cell) sum(ismember(trail_it{cell}, ix_label(:, j))), 1:size(it_resp, 1));
    num_trail_pfc = arrayfun(@(cell) sum(ismember(trail_pfc{cell}, ix_label(:, j))), 1:size(pfc_resp, 1));

    clc; fprintf('Level SVM %d of %d\n', gi, 4);

    % Define time indices and stimulus count for IT and PFC
    numStim_it = 100;
    numStim_pfc = 24;
    tim_ix_1 = 570:700;
    tim_ix_2 = 500:850;

    % Collect IT and PFC responses
    IT = collect_responses(trail_it, it_resp, num_trail_it, ix_label, numStim_it, tim_ix_1);
    PFC = collect_responses(trail_pfc, pfc_resp, num_trail_pfc, ix_label, numStim_pfc, tim_ix_2);

    % Labeling groups
    grp_it = repelem((1:size(ix_label, 2))', numStim_it);
    grp_pfc = repelem((1:size(ix_label, 2))', numStim_pfc);

    % SVM decoding for IT and PFC
    t_res.it(gi) = gen_fx_get_svm_half_trail_boot(grp_it, IT, rate, rep, boot);
    t_res.pfc(gi) = gen_fx_get_svm_half_trail_boot(grp_pfc, PFC, rate, rep, boot);

    % Shuffle-based SVM for comparison
    t_res_shuffle.it(gi) = gen_fx_get_svm_half_trail_sh(grp_it, IT, rate, rep);
    t_res_shuffle.pfc(gi) = gen_fx_get_svm_half_trail_sh(grp_pfc, PFC, rate, rep);
end

% Data preparation for bar plot
data_it = cell2mat(arrayfun(@(x) x.pt, t_res.it, 'UniformOutput', false));
data_pfc = cell2mat(arrayfun(@(x) x.pt, t_res.pfc, 'UniformOutput', false));

% Bar plot with error bars for IT
figure;
bar(mean(data_it, 2));
hold on;
errorbar(mean(data_it, 2), std(data_it, [], 2), 'k.', 'LineWidth', 1);
errorbar(mean(data_it_shuffle, 2), std(data_it_shuffle, [], 2), 'r.', 'LineWidth', 1);
xlabel('Group Index'); ylabel('Decoding Accuracy'); title('IT Decoding Accuracy with Shuffle Data');
legend('Actual Data', 'Shuffle Data');
hold off;

% Bar plot with error bars for vlPFC
figure;
bar(mean(data_pfc, 2));
hold on;
errorbar(mean(data_pfc, 2), std(data_pfc, [], 2), 'k.', 'LineWidth', 1);
errorbar(mean(data_pfc_shuffle, 2), std(data_pfc_shuffle, [], 2), 'r.', 'LineWidth', 1);
xlabel('Group Index'); ylabel('Decoding Accuracy'); title('vlPFC Decoding Accuracy with Shuffle Data');
legend('Actual Data', 'Shuffle Data');
hold off;

% ANOVA for IT and PFC data
[p, tbl, stats] = anova1(data_pfc);
[p, tbl, stats] = anova1(data_it);
fprintf('ANOVA results: F = %f, p = %f\n', tbl{2, 6}, p);
[c, m, h, nms] = multcompare(stats, 'CType', 'tukey-kramer');

% Display significant pairs
significant_pairs = c(c(:, 6) < 0.05, 1:2);
significant_differences = nms(significant_pairs);
disp('Pairs with significant differences:'); disp(significant_differences);
