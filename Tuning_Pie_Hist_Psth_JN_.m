% In the name of Allah
% Author: Jalaledin Noroozi
% Email: jalaledin.noroozi@gmail.com
% Code for PNAS article titled: "Frontotemporal Network Contribution to Occluded Face Processing"
% Tracking #: 2024-07457R

%% Load Data and Initialize
clc; clear; close all;
load('data.mat');  % Load data file
addpath('Functions');  % Add function directory to path

%% Define Parameters and Conditions
[Res_z, cond_z] = all_conditions_zebel();
[Res_j, cond_j] = all_conditions_jenab();

% Define condition indices for different occlusion levels
no_fa_nt = find(Res_j(:,1) < 6 & Res_j(:,4) > 0)';
nt_ = reshape(no_fa_nt, 8, []);
on_ = 311:350;
on_fa_nt = reshape(on_, 8, []);
le_oc = [nt_, on_fa_nt];

% Reorder occlusion levels
lev_occ = [];
for j = [1, 4, 3, 7, 2, 6, 5, 8]
    lev_occ = [lev_occ; le_oc(j,:)];
end

% Organize occlusion conditions into four levels
level_occ_nat = lev_occ';
level_occ_4l = [reshape(lev_occ(1:2,:),[],1), reshape(lev_occ(3:4,:),[],1), ...
                reshape(lev_occ(5:6,:),[],1), reshape(lev_occ(7:8,:),[],1)];
ind_ss = ismember(cond_j', reshape(level_occ_4l, 1, []));

clearvars -except pfc it level_occ_4l cond_j

%% Tuning Curve for Active Time for Best Neurons
% Define sets of best neurons for IT and PFC
IT_best_neuron = [35,38,40,45,46,62,64,106:119,121:184];
PFC_best_neuron = [24:29,83:96,98:121,123:132,135:137,139:159,161];

% Define normalization function
norm_max_min = @(x) (x - min(x)) ./ (max(x) - min(x));

% Define time windows for analysis
ActiveTime = 500:850;
EarlyTime = 500:650;
LateTime = 670:850;

% Analyze IT Best Neurons
for ne = 1:length(IT_best_neuron)
    for jj = 1:4
        cond_indices = ismember(cond_j', reshape(level_occ_4l(:,jj), 1, []));
        TunBestNeuronIT(ne,:,jj) = nanmean(it{IT_best_neuron(ne)}(cond_indices, ActiveTime), 2);
        TunBestNeuronIT_e(ne,:,jj) = nanmean(it{IT_best_neuron(ne)}(cond_indices, EarlyTime), 2);
        TunBestNeuronIT_l(ne,:,jj) = nanmean(it{IT_best_neuron(ne)}(cond_indices, LateTime), 2);
    end
    
    % Normalize and calculate slope coefficients
    resMeanIT(ne,:) = norm_max_min(nanmean(squeeze(TunBestNeuronIT(ne,:,:))));
    resMeanIT_e(ne,:) = norm_max_min(nanmean(squeeze(TunBestNeuronIT_e(ne,:,:))));
    resMeanIT_l(ne,:) = norm_max_min(nanmean(squeeze(TunBestNeuronIT_l(ne,:,:))));
    
    coef_best_it(ne) = polyfit(1:4, resMeanIT(ne,:), 1);
    coef_best_it_e(ne) = polyfit(1:4, resMeanIT_e(ne,:), 1);
    coef_best_it_l(ne) = polyfit(1:4, resMeanIT_l(ne,:), 1);
end

% Analyze PFC Best Neurons (similar structure as IT neurons)
for ne = 1:length(PFC_best_neuron)
    for jj = 1:4
        cond_indices = ismember(cond_j', reshape(level_occ_4l(:,jj), 1, []));
        TunBestNeuronPFC(ne,:,jj) = nanmean(pfc{PFC_best_neuron(ne)}(cond_indices, ActiveTime), 2);
        TunBestNeuronPFC_e(ne,:,jj) = nanmean(pfc{PFC_best_neuron(ne)}(cond_indices, EarlyTime), 2);
        TunBestNeuronPFC_l(ne,:,jj) = nanmean(pfc{PFC_best_neuron(ne)}(cond_indices, LateTime), 2);
    end
    
    % Normalize and calculate slope coefficients
    resMeanPFC(ne,:) = norm_max_min(nanmean(squeeze(TunBestNeuronPFC(ne,:,:))));
    resMeanPFC_e(ne,:) = norm_max_min(nanmean(squeeze(TunBestNeuronPFC_e(ne,:,:))));
    resMeanPFC_l(ne,:) = norm_max_min(nanmean(squeeze(TunBestNeuronPFC_l(ne,:,:))));
    
    coef_best_pfc(ne) = polyfit(1:4, resMeanPFC(ne,:), 1);
    coef_best_pfc_e(ne) = polyfit(1:4, resMeanPFC_e(ne,:), 1);
    coef_best_pfc_l(ne) = polyfit(1:4, resMeanPFC_l(ne,:), 1);
end

%% Tuning Curve for Active Time for All Neurons
% Repeat analysis for all neurons in both IT and PFC
for ne = 1:length(it)
    for jj = 1:4
        cond_indices = ismember(cond_j', reshape(level_occ_4l(:,jj), 1, []));
        TunAllIT(ne,:,jj) = nanmean(it{ne}(cond_indices, ActiveTime), 2);
        TunAllIT_e(ne,:,jj) = nanmean(it{ne}(cond_indices, EarlyTime), 2);
        TunAllIT_l(ne,:,jj) = nanmean(it{ne}(cond_indices, LateTime), 2);
    end
    
    resMeanAllIT(ne,:) = norm_max_min(nanmean(squeeze(TunAllIT(ne,:,:))));
    coef_it(ne) = glmfit(1:4, resMeanAllIT(ne,:), 'normal');
end

for ne = 1:length(pfc)
    for jj = 1:4
        cond_indices = ismember(cond_j', reshape(level_occ_4l(:,jj), 1, []));
        TunAllPFC(ne,:,jj) = nanmean(pfc{ne}(cond_indices, ActiveTime), 2);
        TunAllPFC_e(ne,:,jj) = nanmean(pfc{ne}(cond_indices, EarlyTime), 2);
        TunAllPFC_l(ne,:,jj) = nanmean(pfc{ne}(cond_indices, LateTime), 2);
    end
    
    resMeanAllPFC(ne,:) = norm_max_min(nanmean(squeeze(TunAllPFC(ne,:,:))));
    coef_pfc(ne) = glmfit(1:4, resMeanAllPFC(ne,:), 'normal');
end

%% Plot Histogram for Active Time Slopes
figure;
histogram(coef_it, 'BinWidth', 0.1, 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none');
title('IT'); xlabel('Slope'); ylabel('Number of neurons');
