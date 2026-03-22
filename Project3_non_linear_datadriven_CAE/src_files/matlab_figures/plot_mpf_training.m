%% PLOT_MPF_TRAINING  Plot training loss curves for MPF (Stage-1 CAE + Stage-2 MLP).
%
%  Reads loss_history.json from the MPF KTH run directory and the MPF MLP
%  run directory and plots semilogy convergence curves.
%
%  Usage:  Run this script in MATLAB from the matlab_figures/ folder.

clear; clc; close all;

%% ========================================================================
%  USER CONFIGURATION — Stage-1 (CAE)
%  ========================================================================

model_name = 'ResNet-CAE (MPF 3-ch)';
json_path  = fullfile('..', 'resnet_cae', 'runs', ...
             '2DMPF_GUV_Nt2424_B16_NL8_LR1e-05_L21e-06_E300', ...
             'loss_history.json');

%% ========================================================================
%  USER CONFIGURATION — Stage-2 (MLP)
%  ========================================================================

mlp_model_name = 'Latent MLP (MPF ResNet-CAE)';
mlp_json_path  = fullfile('..', 'latent_mlp', 'runs', ...
                 'MPF_KTH_MLP_HD128_LD128_LR0.001_E5000', ...
                 'loss_history.json');

%% ========================================================================
%  PLOT — Stage-1 Training Loss (CAE)
%  ========================================================================

fprintf('=== %s ===\n', model_name);
raw = fileread(json_path);
d   = jsondecode(raw);

fields = fieldnames(d);

figure('Position', [100 100 600 500]);

colors = lines(numel(fields));
for k = 1:numel(fields)
    vals   = d.(fields{k});
    epochs = 1:numel(vals);
    semilogy(epochs, vals, 'LineWidth', 1.5, 'Color', colors(k, :)); hold on;

    % Print summary
    [best_val, best_ep] = min(vals);
    fprintf('  %s: %d epochs | final=%.4e | best=%.4e (ep %d)\n', ...
        fields{k}, numel(vals), vals(end), best_val, best_ep);
end

xlabel('Epoch', 'FontSize', 20);
ylabel('Loss (log)', 'FontSize', 20);
title(model_name, 'FontSize', 20);
legend(fields, 'Location', 'northeast', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 20);

%% ========================================================================
%  PLOT — Stage-2 Training Loss (MLP)
%  ========================================================================

fprintf('=== %s ===\n', mlp_model_name);
raw_mlp = fileread(mlp_json_path);
d_mlp   = jsondecode(raw_mlp);

fields_mlp = fieldnames(d_mlp);

figure('Position', [750 100 600 500]);

colors_mlp = lines(numel(fields_mlp));
for k = 1:numel(fields_mlp)
    vals   = d_mlp.(fields_mlp{k});
    epochs = 1:numel(vals);
    semilogy(epochs, vals, 'LineWidth', 1.5, 'Color', colors_mlp(k, :)); hold on;

    % Print summary
    [best_val, best_ep] = min(vals);
    fprintf('  %s: %d epochs | final=%.4e | best=%.4e (ep %d)\n', ...
        fields_mlp{k}, numel(vals), vals(end), best_val, best_ep);
end

xlabel('Epoch', 'FontSize', 20);
ylabel('Loss (log)', 'FontSize', 20);
title(mlp_model_name, 'FontSize', 20);
legend(fields_mlp, 'Location', 'northeast', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 20);
