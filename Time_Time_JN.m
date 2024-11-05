% In the name of Allah
% Author: Jalaledin Noroozi
% Email: jalaledin.noroozi@gmail.com
% Code for PNAS article: "Frontotemporal Network Contribution to Occluded Face Processing"

%% Load Data
% Load low-occlusion and full-occlusion data for analysis
load('Time_Time_Intact_LowOc.mat')
low = t_res_perm;

load('Time_Time_Full_LowOc.mat')
full = t_res_perm;

% Compute mean response for IT and PFC in low and full occlusion
it_low = mean(low.it.pt, 3);
pfc_low = mean(low.pfc.pt, 3);
it_full = mean(full.it.pt, 3);
pfc_full = mean(full.pfc.pt, 3);

% Set time vector for analysis
ti_ = mean(tim_tar, 2) - 500;

%% Statistical Analysis for Specific Time Windows
% Perform statistical analysis for specific time windows on IT data

tst1 = 70; tend1 = 200;
ind_tim = (ti_ > tst1 & ti_ < tend1);
statReprt(reshape(it_low(ind_tim, ind_tim), [], 1), 3, 'std')

tst2 = 200; tend2 = 350;
ind_tim = (ti_ > tst2 & ti_ < tend2);
statReprt(reshape(it_low(ind_tim, ind_tim), [], 1), 3, 'std')

%% Wilcoxon Signed-Rank Test on IT Data
% Conduct Wilcoxon signed-rank test for two time windows on IT data

sig_t = (ti_ > 70 & ti_ < 200);
numbers1 = it_low(sig_t, sig_t);
numbers2 = it_full(sig_t, sig_t);
[p1, ~] = signrank(numbers1(:), numbers2(:), 'tail', 'both');

sig_t = (ti_ > 200 & ti_ < 350);
numbers1 = it_low(sig_t, sig_t);
numbers2 = it_full(sig_t, sig_t);
[p2, ~] = signrank(numbers1(:), numbers2(:), 'tail', 'both');

%% Compute Significance for Low Occlusion Data
% Calculate p-values for each time point in IT and PFC for low occlusion data

it_low_perm = mean(low.it.pt_perm, 4);
pfc_low_perm = mean(low.pfc.pt_perm, 4);

pval_it_low = arrayfun(@(i, j) mean(it_low_perm(i, j, :) >= it_low(i, j)), ...
                       1:size(it_low, 1), 1:size(it_low, 2));
pval_it_low(pval_it_low <= 0.001) = 1;
pval_it_low(pval_it_low > 0.001) = NaN;

%% Plot IT Low Occlusion Response
% Generate smoothed IT response image for low occlusion

figure;
imagesc(ti_, ti_, imgaussfilt(it_low, [2, 1]));
colormap(plasma);
colorbar('Label', 'Accuracy');
set(gca, 'YDir', 'normal', 'box', 'off', 'TickDir', 'out', 'FontSize', 12);
caxis([0.1 0.4]);
title('IT Intact-Low Oc.');
xlabel('Testing time (ms)');
ylabel('Training time (ms)');

% Add significance boundaries
significantPatch = pval_it_low < 0.01;
boundaries = bwboundaries(significantPatch);
hold on;
for k = 1:numel(boundaries)
    plot(ti_(boundaries{k}(:, 2)), ti_(boundaries{k}(:, 1)), 'w-', 'LineWidth', 0.25);
end
