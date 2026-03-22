%% PLOT_MPF_PIPELINE  Plot MPF full-pipeline predictions vs ground truth.
%
%  Reads pipeline_predictions.npz from full_pipeline/runs/MPF_KTH_prediction/
%  and displays True / Predicted / |Error| for each channel (G, U, V).
%
%  NPZ contents:
%    originals     (N, 3, 256, 512)  — normalized [0, 1]
%    predictions   (N, 3, 256, 512)
%    channel_names {'G (gas)', 'U (vel-x)', 'V (vel-y)'}
%    case_ids      (N,)
%    timesteps     (N,)
%    params        (N, 4)  — [DR, VR, Re, Bo]
%    mae, mse, linf  (N,)  — per-sample metrics
%
%  Usage:  Run this script in MATLAB from the matlab_figures/ folder.

clear; clc; close all;

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

model_name = 'ResNet-CAE + MLP (MPF)';
npz_path   = fullfile('..', 'full_pipeline', 'runs', ...
             'MPF_KTH_prediction', 'pipeline_predictions.npz');

% Which samples to plot (1-based). Empty = plot all.
SAMPLE_INDICES = [];

% Max samples per figure page (per channel)
MAX_PER_FIG = 5;

%% ========================================================================
%  LOAD DATA
%  ========================================================================

fprintf('=== %s ===\n', model_name);
d = load_npz(npz_path);

originals      = d.originals;       % (N, 3, 256, 512)
predictions    = d.predictions;
channel_names  = d.channel_names;   % cell array of 3 strings
case_ids       = d.case_ids;        % (N,)
timesteps      = d.timesteps;       % (N,)
params         = d.params;          % (N, 4)
mae_arr        = d.mae;             % (N,)
mse_arr        = d.mse;             % (N,)

n_total    = size(originals, 1);
n_channels = size(originals, 2);

if isempty(SAMPLE_INDICES)
    SAMPLE_INDICES = 1:n_total;
end

fprintf('  Loaded: %d samples, %d channels\n', n_total, n_channels);

% %% ========================================================================
% %  PLOT — Per-channel figures (True / Predicted / |Error|)
% %  ========================================================================
% 
% for ch = 1:n_channels
%     ch_name = channel_names{ch};
%     idx = SAMPLE_INDICES;
%     n_samples = numel(idx);
% 
%     n_pages = ceil(n_samples / MAX_PER_FIG);
% 
%     for page = 1:n_pages
%         i_start  = (page - 1) * MAX_PER_FIG + 1;
%         i_end    = min(page * MAX_PER_FIG, n_samples);
%         page_idx = idx(i_start:i_end);
%         n_s      = numel(page_idx);
% 
%         figure('Position', [50 50 320*n_s 900]);
% 
%         for j = 1:n_s
%             si   = page_idx(j);
%             orig = squeeze(originals(si, ch, :, :));
%             pred = squeeze(predictions(si, ch, :, :));
%             err  = abs(orig - pred);
%             vmin = min(orig(:));
%             vmax = max(orig(:));
% 
%             cid = case_ids(si);
%             t   = timesteps(si);
%             p   = params(si, :);
% 
%             % --- True ---
%             subplot(3, n_s, j);
%             imagesc(orig');
%             axis xy equal tight;
%             clim([vmin vmax]);
%             colormap(gca, 'parula');
%             title(sprintf('case %d, t=%d', cid, t), 'FontSize', 16);
%             set(gca, 'XTick', [], 'YTick', []);
%             if j == 1
%                 ylabel('True', 'FontSize', 18, 'FontWeight', 'bold');
%             end
% 
%             % --- Predicted ---
%             subplot(3, n_s, n_s + j);
%             imagesc(pred');
%             axis xy equal tight;
%             clim([vmin vmax]);
%             colormap(gca, 'parula');
%             mse_ch = mean((orig(:) - pred(:)).^2);
%             title(sprintf('MSE = %.2e', mse_ch), 'FontSize', 16);
%             set(gca, 'XTick', [], 'YTick', []);
%             if j == 1
%                 ylabel('Predicted', 'FontSize', 18, 'FontWeight', 'bold');
%             end
% 
%             % --- Error ---
%             subplot(3, n_s, 2*n_s + j);
%             imagesc(err');
%             axis xy equal tight;
%             colorbar;
%             colormap(gca, 'hot');
%             title(sprintf('Max err = %.3f', max(err(:))), 'FontSize', 16);
%             set(gca, 'XTick', [], 'YTick', []);
%             if j == 1
%                 ylabel('|Error|', 'FontSize', 18, 'FontWeight', 'bold');
%             end
%         end
% 
%         if n_pages > 1
%             sgtitle(sprintf('%s — %s (page %d/%d)', model_name, ch_name, ...
%                 page, n_pages), 'FontSize', 20, 'FontWeight', 'bold');
%         else
%             sgtitle(sprintf('%s — %s', model_name, ch_name), ...
%                 'FontSize', 20, 'FontWeight', 'bold');
%         end
%     end
% end

%% ========================================================================
%  COMBINED SUMMARY FIGURE — one row per case, columns = timesteps
%  ========================================================================

unique_cases = unique(case_ids, 'stable');
n_cases = numel(unique_cases);

for ch = 1:n_channels
    ch_name = channel_names{ch};

    for ci = 1:n_cases
        cid = unique_cases(ci);
        mask = (case_ids == cid);
        c_idx = find(mask);
        n_t = numel(c_idx);

        figure('Position', [50 50 300*n_t 900]);

        for j = 1:n_t
            si   = c_idx(j);
            orig = squeeze(originals(si, ch, :, :));
            pred = squeeze(predictions(si, ch, :, :));
            err  = abs(orig - pred);
            vmin = min(orig(:));
            vmax = max(orig(:));
            t    = timesteps(si);

            subplot(3, n_t, j);
            imagesc(orig');
            axis xy equal tight;
            clim([vmin vmax]);
            colormap(gca, 'parula');
            title(sprintf('t = %d', t), 'FontSize', 16);
            set(gca, 'XTick', [], 'YTick', []);
            if j == 1, ylabel('True', 'FontSize', 18, 'FontWeight', 'bold'); end

            subplot(3, n_t, n_t + j);
            imagesc(pred');
            axis xy equal tight;
            clim([vmin vmax]);
            colormap(gca, 'parula');
            mse_ch = mean((orig(:) - pred(:)).^2);
            title(sprintf('MSE = %.2e', mse_ch), 'FontSize', 16);
            set(gca, 'XTick', [], 'YTick', []);
            if j == 1, ylabel('Predicted', 'FontSize', 18, 'FontWeight', 'bold'); end

            subplot(3, n_t, 2*n_t + j);
            imagesc(err');
            axis xy equal tight;
            colorbar;
            colormap(gca, 'hot');
            title(sprintf('Max = %.3f', max(err(:))), 'FontSize', 16);
            set(gca, 'XTick', [], 'YTick', []);
            if j == 1, ylabel('|Error|', 'FontSize', 18, 'FontWeight', 'bold'); end
        end

        p = params(c_idx(1), :);
        sgtitle(sprintf('%s — %s\nCase %d  [DR=%d, VR=%d, Re=%d, Bo=%d]', ...
            model_name, ch_name, cid, p(1), p(2), p(3), p(4)), ...
            'FontSize', 18, 'FontWeight', 'bold');
    end
end

%% ========================================================================
%  PRINT SUMMARY TABLE
%  ========================================================================

fprintf('\n  Per-sample metrics (normalized space):\n');
fprintf('  %-4s %-6s %-4s %-20s %-12s %-12s\n', ...
    '#', 'Case', 't', 'Params', 'MSE', 'MAE');
fprintf('  %s\n', repmat('-', 1, 55));

for j = 1:n_total
    cid = case_ids(j);
    t   = timesteps(j);
    p   = params(j, :);
    orig_j = squeeze(originals(j, :, :, :));
    pred_j = squeeze(predictions(j, :, :, :));
    mse_j  = mean((orig_j(:) - pred_j(:)).^2);
    mae_j  = mean(abs(orig_j(:) - pred_j(:)));
    fprintf('  %-4d %-6d %-4d [%3d,%2d,%3d,%3d]  %12.4e %12.4e\n', ...
        j, cid, t, p(1), p(2), p(3), p(4), mse_j, mae_j);
end

fprintf('\n  Per-channel metrics (all samples):\n');
for ch = 1:n_channels
    o = originals(:, ch, :, :);
    r = predictions(:, ch, :, :);
    mse_ch = mean((o(:) - r(:)).^2);
    mae_ch = mean(abs(o(:) - r(:)));
    fprintf('  Ch%d (%s): MSE=%.4e, MAE=%.4e\n', ch, channel_names{ch}, ...
        mse_ch, mae_ch);
end

fprintf('\n  Overall: MSE=%.4e, MAE=%.4e\n', ...
    mean(mse_arr), mean(mae_arr));
