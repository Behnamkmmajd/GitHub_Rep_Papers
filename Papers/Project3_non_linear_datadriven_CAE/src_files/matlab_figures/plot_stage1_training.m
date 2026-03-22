%% PLOT_TRAINING  Plot training loss curves (Stage-1 CAE and/or Stage-2 MLP).
%
%  Reads loss_history.json from each run directory and plots semilogy
%  convergence curves.  Handles CNN format (Train/Val), KTH format
%  (Total), and MLP format (Train/Val).
%
%  Usage:  Run this script in MATLAB from the matlab_figures/ folder.

clear; clc; close all;

%% ========================================================================
%  USER CONFIGURATION — Stage-1 (CAE)
%  ========================================================================

% model_name = 'CNN-ROM (Wang)';
% json_path  = fullfile('..', 'simple_cae', 'runs', ...
%              'Advection_Ch1_Nt240_B32_LD4_LR0.0001_E5000_Norm', ...
%              'loss_history.json');

model_name = 'ResNet-CAE';
json_path  = fullfile('..', 'resnet_cae', 'runs', ...
             'Advection_Ch1_Nt240_B8_NL6_LR1e-05_L21e-06_E300_Norm', ...
             'loss_history.json');

%% ========================================================================
%  USER CONFIGURATION — Stage-2 (MLP)
%  ========================================================================

% mlp_model_name = 'Latent MLP (CNN encoder)';
% mlp_json_path  = fullfile('..', 'latent_mlp', 'runs', ...
%                  'Advection_Ch1_CNN_MLP_HD64_LD4_LR0.001_E5000', ...
%                  'loss_history.json');

mlp_model_name = 'Latent MLP (ResNet-CAE)';
mlp_json_path  = fullfile('..', 'latent_mlp', 'runs', ...
                 'Advection_Ch1_KTH_MLP_HD64_LD1024_LR0.001_E5000', ...
                 'loss_history.json');

%% ========================================================================
%  PLOT — Stage-1 (CAE)
%  ========================================================================

fprintf('=== %s ===\n', model_name);
raw = fileread(json_path);
d   = jsondecode(raw);

figure('Position', [100 100 600 500]);
fields = fieldnames(d);
for k = 1:numel(fields)
    vals = d.(fields{k});
    semilogy(1:numel(vals), vals, 'LineWidth', 1.2); hold on;
end

xlabel('Epoch', 'FontSize', 20);
ylabel('Loss (log)', 'FontSize', 20);
title(model_name, 'FontSize', 20);
legend(fields, 'Location', 'northeast', 'FontSize', 20);
grid on;
set(gca, 'FontSize', 20);

%% ========================================================================
%  PLOT — Stage-2 (MLP)
%  ========================================================================

fprintf('=== %s ===\n', mlp_model_name);
raw_mlp = fileread(mlp_json_path);
d_mlp   = jsondecode(raw_mlp);

figure('Position', [750 100 600 500]);
fields_mlp = fieldnames(d_mlp);
for k = 1:numel(fields_mlp)
    vals = d_mlp.(fields_mlp{k});
    semilogy(1:numel(vals), vals, 'LineWidth', 1.2); hold on;
end

xlabel('Epoch', 'FontSize', 20);
ylabel('Loss (log)', 'FontSize', 20);
title(mlp_model_name, 'FontSize', 20);
legend(fields_mlp, 'Location', 'northeast', 'FontSize', 20);
grid on;
set(gca, 'FontSize', 20);
