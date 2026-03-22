%% PLOT_STAGE1_RECONSTRUCTION  Plot CAE reconstruction quality.
%
%  Reads reconstruction_test.npz produced by each project's reconstruct.py
%  and displays original / reconstructed / |error| grids.
%
%  Usage:  Run this script in MATLAB from the matlab_figures/ folder.

clear; clc; close all;

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

% model_name = 'CNN-ROM (Wang)';
% npz_path   = fullfile('..', 'simple_cae', 'runs', ...
%              'Advection_Ch1_Nt240_B32_LD4_LR0.0001_E5000_Norm', ...
%              'reconstruction_test.npz');

model_name = 'ResNet-CAE';
npz_path   = fullfile('..', 'resnet_cae', 'runs', ...
             'Advection_Ch1_Nt240_B8_NL6_LR1e-05_L21e-06_E300_Norm', ...
             'reconstruction_test.npz');

%% ========================================================================
%  PLOT
%  ========================================================================

fprintf('=== %s ===\n', model_name);
d = load_npz(npz_path);

originals       = d.originals;        % (N, C, H, W)
reconstructions = d.reconstructions;
n_samples  = size(originals, 1);
n_channels = size(originals, 2);

% Read source labels (val/test) if available
if isfield(d, 'sources')
    sources = d.sources;  % cell array of 'val'/'test'
else
    sources = repmat({'unseen'}, 1, n_samples);
end

for ch = 1:n_channels
    figure('Position', [100 100 300*n_samples 900]);

    for i = 1:n_samples
        orig  = squeeze(originals(i, ch, :, :));
        recon = squeeze(reconstructions(i, ch, :, :));
        err   = abs(orig - recon);
        vmin  = min(orig(:));
        vmax  = max(orig(:));

        % --- Original ---
        subplot(3, n_samples, i);
        imagesc(orig');
        axis xy equal tight;
        clim([vmin vmax]);
        colormap(gca, 'parula');
        title(sprintf('#%d (%s)', i, sources{i}), 'FontSize', 20);
        set(gca, 'XTick', [], 'YTick', []);
        if i == 1
            ylabel('Original', 'FontSize', 20, 'FontWeight', 'bold');
        end

        % --- Reconstructed ---
        subplot(3, n_samples, n_samples + i);
        imagesc(recon');
        axis xy equal tight;
        clim([vmin vmax]);
        colormap(gca, 'parula');
        mse_val = mean((orig(:) - recon(:)).^2);
        title(sprintf('MSE = %.1e', mse_val), 'FontSize', 20);
        set(gca, 'XTick', [], 'YTick', []);
        if i == 1
            ylabel('Reconstructed', 'FontSize', 20, 'FontWeight', 'bold');
        end

        % --- Error ---
        subplot(3, n_samples, 2*n_samples + i);
        imagesc(err');
        axis xy equal tight;
        colorbar;
        colormap(gca, 'hot');
        title(sprintf('Max err = %.3f', max(err(:))), 'FontSize', 20);
        set(gca, 'XTick', [], 'YTick', []);
        if i == 1
            ylabel('|Error|', 'FontSize', 20, 'FontWeight', 'bold');
        end
    end

    if n_channels > 1
        sgtitle(sprintf('%s — Channel %d', model_name, ch), ...
            'FontSize', 20, 'FontWeight', 'bold');
    else
        sgtitle(sprintf('%s — Reconstruction', model_name), ...
            'FontSize', 20, 'FontWeight', 'bold');
    end
end
