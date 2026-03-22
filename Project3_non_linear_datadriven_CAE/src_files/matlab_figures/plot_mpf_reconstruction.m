%% PLOT_MPF_RECONSTRUCTION  Plot MPF 3-channel CAE reconstruction quality.
%
%  Reads reconstruction_test.npz produced by reconstruct_mpf.py and
%  displays Original / Reconstructed / |Error| for each channel (G, U, V).
%
%  NPZ contents:
%    originals       (N, 3, 256, 512)   — normalized to [0, 1]
%    reconstructions (N, 3, 256, 512)
%    channel_names   {'G (gas)', 'U (vel-x)', 'V (vel-y)'}
%    sources         {'train', 'val', 'test', ...}
%    timesteps       [0, 50, 100, ...]
%    params          (N, 4)  — [DR, VR, Re, Bo]
%
%  Usage:  Run this script in MATLAB from the matlab_figures/ folder.

clear; clc; close all;

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

model_name = 'KTH ResNet-CAE (MPF)';
npz_path   = fullfile('..', 'resnet_cae', 'runs', ...
             '2DMPF_GUV_Nt2424_B16_NL8_LR1e-05_L21e-06_E300', ...
             'reconstruction_test.npz');

% Which samples to plot (1-based, into the NPZ arrays)
% Default: plot all
SAMPLE_INDICES = [];   % empty = plot all

% Maximum samples per figure (to avoid overcrowding)
MAX_PER_FIG = 5;

%% ========================================================================
%  LOAD DATA
%  ========================================================================

fprintf('=== %s ===\n', model_name);
d = load_npz(npz_path);

originals       = d.originals;         % (N, 3, 256, 512)
reconstructions = d.reconstructions;
channel_names   = d.channel_names;     % cell array
sources         = d.sources;           % cell array
timesteps       = d.timesteps;         % (N, 1) or (N,)
params          = d.params;            % (N, 4)

n_total    = size(originals, 1);
n_channels = size(originals, 2);

if isempty(SAMPLE_INDICES)
    SAMPLE_INDICES = 1:n_total;
end

fprintf('  Loaded: %d samples, %d channels\n', n_total, n_channels);

%% ========================================================================
%  PLOT — Per-channel figures
%  ========================================================================

for ch = 1:n_channels
    ch_name = channel_names{ch};
    idx = SAMPLE_INDICES;
    n_samples = numel(idx);

    % Split into pages if too many samples
    n_pages = ceil(n_samples / MAX_PER_FIG);

    for page = 1:n_pages
        i_start = (page - 1) * MAX_PER_FIG + 1;
        i_end   = min(page * MAX_PER_FIG, n_samples);
        page_idx = idx(i_start:i_end);
        n_s = numel(page_idx);

        figure('Position', [50 50 320*n_s 900]);

        for j = 1:n_s
            si   = page_idx(j);
            orig  = squeeze(originals(si, ch, :, :));
            recon = squeeze(reconstructions(si, ch, :, :));
            err   = abs(orig - recon);
            vmin  = min(orig(:));
            vmax  = max(orig(:));

            src = sources{si};
            t   = timesteps(si);
            p   = params(si, :);

            % --- Original ---
            subplot(3, n_s, j);
            imagesc(orig');
            axis xy equal tight;
            clim([vmin vmax]);
            colormap(gca, 'parula');
            title(sprintf('%s t=%d', src, t), 'FontSize', 16);
            set(gca, 'XTick', [], 'YTick', []);
            if j == 1
                ylabel('Original', 'FontSize', 18, 'FontWeight', 'bold');
            end

            % --- Reconstructed ---
            subplot(3, n_s, n_s + j);
            imagesc(recon');
            axis xy equal tight;
            clim([vmin vmax]);
            colormap(gca, 'parula');
            mse_val = mean((orig(:) - recon(:)).^2);
            title(sprintf('MSE = %.2e', mse_val), 'FontSize', 16);
            set(gca, 'XTick', [], 'YTick', []);
            if j == 1
                ylabel('Reconstructed', 'FontSize', 18, 'FontWeight', 'bold');
            end

            % --- Error ---
            subplot(3, n_s, 2*n_s + j);
            imagesc(err');
            axis xy equal tight;
            colorbar;
            colormap(gca, 'hot');
            title(sprintf('Max err = %.3f', max(err(:))), 'FontSize', 16);
            set(gca, 'XTick', [], 'YTick', []);
            if j == 1
                ylabel('|Error|', 'FontSize', 18, 'FontWeight', 'bold');
            end
        end

        if n_pages > 1
            sgtitle(sprintf('%s — %s (page %d/%d)', model_name, ch_name, ...
                page, n_pages), 'FontSize', 20, 'FontWeight', 'bold');
        else
            sgtitle(sprintf('%s — %s', model_name, ch_name), ...
                'FontSize', 20, 'FontWeight', 'bold');
        end
    end
end

%% ========================================================================
%  PRINT SUMMARY TABLE
%  ========================================================================

fprintf('\n  Per-sample metrics:\n');
fprintf('  %-6s %-6s %-4s %-20s %-12s %-12s\n', ...
    'Idx', 'Source', 't', 'Params', 'MSE', 'MAE');
fprintf('  %s\n', repmat('-', 1, 62));

for j = 1:n_total
    orig_j  = squeeze(originals(j, :, :, :));
    recon_j = squeeze(reconstructions(j, :, :, :));
    mse_j = mean((orig_j(:) - recon_j(:)).^2);
    mae_j = mean(abs(orig_j(:) - recon_j(:)));
    p = params(j, :);
    fprintf('  %-6d %-6s %-4d [%3d,%2d,%3d,%3d]  %12.4e %12.4e\n', ...
        j, sources{j}, timesteps(j), p(1), p(2), p(3), p(4), mse_j, mae_j);
end

fprintf('\n  Per-channel metrics (all samples):\n');
for ch = 1:n_channels
    o = originals(:, ch, :, :);
    r = reconstructions(:, ch, :, :);
    mse_ch = mean((o(:) - r(:)).^2);
    mae_ch = mean(abs(o(:) - r(:)));
    fprintf('  Ch%d (%s): MSE=%.4e, MAE=%.4e\n', ch, channel_names{ch}, ...
        mse_ch, mae_ch);
end
