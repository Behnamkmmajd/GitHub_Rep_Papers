%% PLOT_PIPELINE_PREDICTION  Plot full-pipeline predictions vs ground truth.
%
%  Reads predictions.npz from full_pipeline/runs/ and plots:
%    Row 1: True u(x,y)
%    Row 2: Predicted u(x,y)
%    Row 3: |Error|
%
%  Usage:  Run this script in MATLAB from the matlab_figures/ folder.

clear; clc; close all;

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

% model_name = 'CNN-ROM Pipeline';
% npz_path   = fullfile('..', 'full_pipeline', 'runs', ...
%              'Advection_CNN_prediction', 'predictions.npz');

model_name = 'ResNet-CAE Pipeline';
npz_path   = fullfile('..', 'full_pipeline', 'runs', ...
             'Advection_KTH_prediction', 'predictions.npz');

% Channel to plot (1-based)
CH = 1;

%% ========================================================================
%  LOAD DATA
%  ========================================================================

fprintf('=== %s ===\n', model_name);
d = load_npz(npz_path);

U_true  = d.U_true;     % (N, C, H, W)
U_pred  = d.U_pred;     % (N, C, H, W)
X_test  = d.X_test;     % (N, 2)  columns: [mu, t]
n_samples = size(U_true, 1);

fprintf('  Snapshots: %d\n', n_samples);
fprintf('  Field size: %dx%d\n', size(U_true, 3), size(U_true, 4));

%% ========================================================================
%  PLOT
%  ========================================================================

figure('Position', [100 100 300*n_samples 900]);

for i = 1:n_samples
    orig = squeeze(U_true(i, CH, :, :));
    pred = squeeze(U_pred(i, CH, :, :));
    err  = abs(orig - pred);
    vmin = min(orig(:));
    vmax = max(orig(:));
    mse_val = mean((orig(:) - pred(:)).^2);

    mu_i = X_test(i, 1);
    t_i  = X_test(i, 2);

    % --- True ---
    subplot(3, n_samples, i);
    imagesc(orig');
    axis xy equal tight;
    clim([vmin vmax]);
    colormap(gca, 'parula');
    title(sprintf('$\\mu$=%.2f, t=%.1f', mu_i, t_i), ...
        'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'XTick', [], 'YTick', []);
    if i == 1
        ylabel('True', 'FontSize', 20, 'FontWeight', 'bold');
    end

    % --- Predicted ---
    subplot(3, n_samples, n_samples + i);
    imagesc(pred');
    axis xy equal tight;
    clim([vmin vmax]);
    colormap(gca, 'parula');
    title(sprintf('MSE = %.1e', mse_val), 'FontSize', 20);
    set(gca, 'XTick', [], 'YTick', []);
    if i == 1
        ylabel('Predicted', 'FontSize', 20, 'FontWeight', 'bold');
    end

    % --- Error ---
    subplot(3, n_samples, 2*n_samples + i);
    imagesc(err');
    axis xy equal tight;
    colorbar;
    colormap(gca, 'hot');
    title(sprintf('Max |err| = %.3f', max(err(:))), 'FontSize', 20);
    set(gca, 'XTick', [], 'YTick', []);
    if i == 1
        ylabel('|Error|', 'FontSize', 20, 'FontWeight', 'bold');
    end
end

sgtitle(sprintf('%s: $(\\mu, t) \\to$ MLP $\\to \\mathbf{z} \\to$ Decoder $\\to \\tilde{u}(x,y)$', model_name), ...
    'Interpreter', 'latex', 'FontSize', 20);
