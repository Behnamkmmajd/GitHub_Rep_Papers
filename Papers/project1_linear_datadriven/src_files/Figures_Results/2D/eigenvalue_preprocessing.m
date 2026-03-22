clc; clear; close all

%% =====================  FIGURE STYLE SETTINGS  =====================
%  Modify these values to change the appearance of all plots uniformly.
%  ---- Colors (reference: Figure 3 – spectral compressibility) ----
col_E   = [0.0  0.0  0.0 ];   % Eulerian      — black
col_L   = [0.0  0.0  0.69];   % Lagrangian    — blue
col_ME  = [0.0  0.61 0.2 ];   % Moving-Euler. — green
col_DNS = [0.8  0.0  0.0 ];   % DNS reference — red
%  ---- Line widths ----
lw      = 3.0;                 % data lines
lw_ax   = 2.5;                 % axes border / ticks
%  ---- Font ----
fontName    = 'Times New Roman';
fontSize    = 38;
lgdFontSize = 34;
%  ---- Axes preset (pass to set(gca, ...)) ----
ax_style = {'FontName', fontName, 'FontSize', fontSize, ...
            'TickDir', 'in', 'LineWidth', lw_ax, 'Box', 'on', ...
            'XMinorTick', 'on', 'YMinorTick', 'on'};
%  ---- Eigenvalue-specific colors ----
col_raw   = col_DNS;               % red    — problematic (mean0 / raw)
col_fix   = col_L;                 % blue   — fixed (mean1 / der1)
col_UC    = [0.4  0.4  0.4 ];      % grey   — unit circle
mk_raw    = 'o';
mk_fix    = 'o';
ms_raw    = 10;
ms_fix    = 10;
%  ====================================================================

%% ========================================================================
%  Figure: Eigenvalue spectra — preprocessing effects on SPDMD stability
%
%  Panel (a): Lagrangian position data (LP, C1)
%             mean0 (raw) vs der1 (FD derivative)
%             Shows catastrophic |mu| >> 1  →  |mu| ≈ 1 after derivative
%
%  Panel (b): Moving-Eulerian gas fraction (MEG, C3)
%             mean0 vs mean1 (mean subtraction)
%             Worst ME gas case: |mu|_max = 1.17, log10(Vand) = 34
%
%  Panel (c): Moving-Eulerian velocity (MEV, C4)
%             mean0 vs mean1 (mean subtraction)
%             Worst ME velocity case: 329/490 modes outside UC, ADMM failed
%  ========================================================================

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
DMD_path  = fullfile(scriptDir, '..', '..', 'Codes', 'DMD', 'Results');
DMD_path  = char(java.io.File(DMD_path).getCanonicalPath());

make_dmd = @(obs, dim, s, e, DR, VR, Re, Bo, var) ...
    fullfile(DMD_path, sprintf( ...
        'Results_DMD_%s_%snew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
        obs, dim, s, e, DR, VR, Re, Bo, var));

%% ---- Read eigenvalue data -----------------------------------------------
% (a) LP position: Re=20, Bo=20, 2D, 300 snapshots
[mu_LP_m0,  ~] = read_eigenvalues(make_dmd('LP','2D',1,300,10,10,20,20,'mean0'));
[mu_LP_d1,  ~] = read_eigenvalues(make_dmd('LP','2D',1,300,10,10,20,20,'der1'));

% (b) ME gas: Re=100, Bo=100, 2D, 500 snapshots  (worst ME gas: |mu|_max=1.17)
[mu_LEG_m0, ~] = read_eigenvalues(make_dmd('LEGi','2D',1,500,10,10,100,100,'mean0'));
[mu_LEG_m1, ~] = read_eigenvalues(make_dmd('LEGi','2D',1,500,10,10,100,100,'mean1'));

% (c) ME velocity: Re=168, Bo=2, 3D, 490 snapshots  (worst ME vel: ADMM failed)
[mu_LEV_m0, ~] = read_eigenvalues(make_dmd('LEUVWi','3D',1,490,10,10,168,2,'mean0'));
[mu_LEV_m1, ~] = read_eigenvalues(make_dmd('LEUVWi','3D',1,490,10,10,168,2,'mean1'));

%% ---- Figure setup -------------------------------------------------------
theta = linspace(0, 2*pi, 500);
uc_x  = cos(theta);
uc_y  = sin(theta);

% ======================================================================
% (a) LP position: mean0 vs derivative
% ======================================================================
figure('Units','centimeters','Position',[2 2 22 14]);
axes('Position', [0.11 0.17 0.84 0.75]); hold on; grid on;
set(gca, 'GridLineWidth', 0.5, 'GridAlpha', 0.15);

% Rotate panel by 90 degrees so the long axis lies horizontally.
plot(uc_y, uc_x, '-', 'Color', col_UC, 'LineWidth', 2, ...
     'HandleVisibility','off');

% Raw eigenvalues (mean0) — stable = faint, unstable = prominent
mu_raw = mu_LP_m0;
stable   = abs(mu_raw) <= 1.001;
unstable = abs(mu_raw) >  1.001;

scatter(imag(mu_raw(stable)),   real(mu_raw(stable)),   ms_raw^2/2, ...
    'MarkerEdgeColor', col_raw, 'MarkerEdgeAlpha', 0.3, ...
    'MarkerFaceColor', col_raw, 'MarkerFaceAlpha', 0.15, ...
    'HandleVisibility','off');
h_raw = scatter(imag(mu_raw(unstable)), real(mu_raw(unstable)), ms_raw^2, ...
    'MarkerEdgeColor', col_raw, 'MarkerFaceColor', col_raw, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);

% Derivative eigenvalues (der1) — all on/near UC
h_fix = scatter(imag(mu_LP_d1), real(mu_LP_d1), ms_fix^2, ...
    'MarkerEdgeColor', col_fix, 'MarkerFaceColor', col_fix, ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);

axis equal;
xlim([-8.5 8.5]); ylim([-1.5 3.0]);
xlabel('Im({\it\mu})','FontName',fontName,'FontSize',fontSize);
ylabel('Re({\it\mu})','FontName',fontName,'FontSize',fontSize);
set(gca, ax_style{:});
legend([h_raw, h_fix], ...
    {sprintf('$\\mathbf{P}\\;(|\\mu|_{\\mathrm{max}}=%.1f)$', max(abs(mu_LP_m0))), ...
     sprintf('$\\dot{\\mathbf{P}}\\;(|\\mu|_{\\mathrm{max}}=%.3f)$', max(abs(mu_LP_d1)))}, ...
    'Interpreter', 'latex', ...
    'FontSize', lgdFontSize, ...
    'Location', 'northeast');
%title('(a) LP position, {\itC}_1','FontName',fontName,'FontSize',fontSize,'FontWeight','normal');

% ======================================================================
% (b) ME gas fraction: mean0 vs mean1
% ======================================================================
figure('Units','centimeters','Position',[26 2 16 16]);
axes('Position', [0.14 0.14 0.80 0.78]); hold on; grid on;
set(gca, 'GridLineWidth', 0.5, 'GridAlpha', 0.15);

plot(uc_x, uc_y, '-', 'Color', col_UC, 'LineWidth', 2, ...
     'HandleVisibility','off');

% mean0 eigenvalues
mu_raw = mu_LEG_m0;
stable   = abs(mu_raw) <= 1.001;
unstable = abs(mu_raw) >  1.001;

scatter(real(mu_raw(stable)),   imag(mu_raw(stable)),   ms_raw^2/2, ...
    'MarkerEdgeColor', col_raw, 'MarkerEdgeAlpha', 0.3, ...
    'MarkerFaceColor', col_raw, 'MarkerFaceAlpha', 0.15, ...
    'HandleVisibility','off');
h_raw = scatter(real(mu_raw(unstable)), imag(mu_raw(unstable)), ms_raw^2, ...
    'MarkerEdgeColor', col_raw, 'MarkerFaceColor', col_raw, ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.9);

% mean1 eigenvalues
h_fix = scatter(real(mu_LEG_m1), imag(mu_LEG_m1), ms_fix^2, ...
    'MarkerEdgeColor', col_fix, 'MarkerFaceColor', col_fix, ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);

axis equal;
lim = 1.25;
xlim([-lim lim]); ylim([-lim lim]);
xlabel('Re({\it\mu})','FontName',fontName,'FontSize',fontSize);
set(gca, ax_style{:});
legend([h_raw, h_fix], ...
    {sprintf('No mean sub.  (|{\\it\\mu}|_{max} = %.2f)', max(abs(mu_LEG_m0))), ...
     sprintf('Mean sub.  (|{\\it\\mu}|_{max} \\approx 1)', max(abs(mu_LEG_m1)))}, ...
    'FontSize',lgdFontSize,'FontName',fontName,'Location','southeast');
title('(b) ME gas fraction, {\itC}_3', ...
    'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');

% % -- Zoom inset for panel (b): top-right corner --
% inset_w = pw * 0.42;   inset_h = ph * 0.42;
% ax_inset_b = axes('Position', ...
%     [pos_b(1) + pos_b(3) - inset_w - 0.005, ...
%      pos_b(2) + pos_b(4) - inset_h - 0.005, ...
%      inset_w, inset_h]);
% hold on; grid on;
% plot(uc_x, uc_y, '-', 'Color', col_UC, 'LineWidth', 1.5, ...
%      'HandleVisibility','off');
% scatter(real(mu_raw(unstable)), imag(mu_raw(unstable)), ms_raw^2*2, ...
%     'MarkerEdgeColor', col_raw, 'MarkerFaceColor', col_raw, ...
%     'MarkerFaceAlpha', 0.8);
% % Show mean1 modes in the zoomed region
% mask_m1 = abs(mu_LEG_m1) > 0.90;
% scatter(real(mu_LEG_m1(mask_m1)), imag(mu_LEG_m1(mask_m1)), ms_fix^2*2, ...
%     'MarkerEdgeColor', col_fix, 'MarkerFaceColor', col_fix, ...
%     'MarkerFaceAlpha', 0.5);
% axis equal;
% xlim([0.85 1.20]); ylim([-0.55 0.55]);
% set(ax_inset_b, 'FontName','Times New Roman','FontSize',11, ...
%     'TickDir','in','LineWidth',1.5,'Box','on');

% ======================================================================
% (c) ME velocity: mean0 vs mean1
% ======================================================================
figure('Units','centimeters','Position',[44 2 16 16]);
axes('Position', [0.14 0.14 0.80 0.78]); hold on; grid on;
set(gca, 'GridLineWidth', 0.5, 'GridAlpha', 0.15);

plot(uc_x, uc_y, '-', 'Color', col_UC, 'LineWidth', 2, ...
     'HandleVisibility','off');

% mean0 eigenvalues
mu_raw = mu_LEV_m0;
stable   = abs(mu_raw) <= 1.001;
unstable = abs(mu_raw) >  1.001;

scatter(real(mu_raw(stable)),   imag(mu_raw(stable)),   ms_raw^2/2, ...
    'MarkerEdgeColor', col_raw, 'MarkerEdgeAlpha', 0.3, ...
    'MarkerFaceColor', col_raw, 'MarkerFaceAlpha', 0.15, ...
    'HandleVisibility','off');
h_raw = scatter(real(mu_raw(unstable)), imag(mu_raw(unstable)), ms_raw^2, ...
    'MarkerEdgeColor', col_raw, 'MarkerFaceColor', col_raw, ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.9);

% mean1 eigenvalues
h_fix = scatter(real(mu_LEV_m1), imag(mu_LEV_m1), ms_fix^2, ...
    'MarkerEdgeColor', col_fix, 'MarkerFaceColor', col_fix, ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.7);

axis equal;
lim = 1.20;
xlim([-lim lim]); ylim([-lim lim]);
xlabel('Re({\it\mu})','FontName',fontName,'FontSize',fontSize);
set(gca, ax_style{:});
legend([h_raw, h_fix], ...
    {sprintf('No mean sub.  (|{\\it\\mu}|_{max} = %.2f)', max(abs(mu_LEV_m0))), ...
     sprintf('Mean sub.  (|{\\it\\mu}|_{max} \\approx 1)', max(abs(mu_LEV_m1)))}, ...
    'FontSize',lgdFontSize,'FontName',fontName,'Location','southeast');
title('(c) ME velocity, {\itC}_4', ...
    'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');

% % -- Zoom inset for panel (c): top-right corner --
% ax_inset_c = axes('Position', ...
%     [pos_c(1) + pos_c(3) - inset_w - 0.005, ...
%      pos_c(2) + pos_c(4) - inset_h - 0.005, ...
%      inset_w, inset_h]);
% hold on; grid on;
% plot(uc_x, uc_y, '-', 'Color', col_UC, 'LineWidth', 1.5, ...
%      'HandleVisibility','off');
% scatter(real(mu_raw(unstable)), imag(mu_raw(unstable)), ms_raw^2*2, ...
%     'MarkerEdgeColor', col_raw, 'MarkerFaceColor', col_raw, ...
%     'MarkerFaceAlpha', 0.8);
% mask_m1 = abs(mu_LEV_m1) > 0.90;
% scatter(real(mu_LEV_m1(mask_m1)), imag(mu_LEV_m1(mask_m1)), ms_fix^2*2, ...
%     'MarkerEdgeColor', col_fix, 'MarkerFaceColor', col_fix, ...
%     'MarkerFaceAlpha', 0.5);
% axis equal;
% xlim([0.95 1.15]); ylim([-0.25 0.25]);
% set(ax_inset_c, 'FontName','Times New Roman','FontSize',11, ...
%     'TickDir','in','LineWidth',1.5,'Box','on');

%% ---- Print summary to console ----
fprintf('\n===== Preprocessing Eigenvalue Summary =====\n');
fprintf('(a) LP position, C1 (Re=20, Bo=20):\n');
fprintf('    Raw:        n=%d, |mu|_max = %.4f, n(|mu|>1) = %d\n', ...
    numel(mu_LP_m0), max(abs(mu_LP_m0)), sum(abs(mu_LP_m0)>1.001));
fprintf('    Derivative: n=%d, |mu|_max = %.6f, n(|mu|>1) = %d\n', ...
    numel(mu_LP_d1), max(abs(mu_LP_d1)), sum(abs(mu_LP_d1)>1.001));
fprintf('(b) ME gas fraction, C3 (Re=100, Bo=100):\n');
fprintf('    mean0: n=%d, |mu|_max = %.4f, n(|mu|>1) = %d, log10(Vand) = %.1f\n', ...
    numel(mu_LEG_m0), max(abs(mu_LEG_m0)), sum(abs(mu_LEG_m0)>1.001), ...
    (numel(mu_LEG_m0)-1)*log10(max(abs(mu_LEG_m0))));
fprintf('    mean1: n=%d, |mu|_max = %.6f, n(|mu|>1) = %d\n', ...
    numel(mu_LEG_m1), max(abs(mu_LEG_m1)), sum(abs(mu_LEG_m1)>1.001));
fprintf('(c) ME velocity, C4 (Re=168, Bo=2):\n');
fprintf('    mean0: n=%d, |mu|_max = %.4f, n(|mu|>1) = %d, log10(Vand) = %.1f\n', ...
    numel(mu_LEV_m0), max(abs(mu_LEV_m0)), sum(abs(mu_LEV_m0)>1.001), ...
    (numel(mu_LEV_m0)-1)*log10(max(abs(mu_LEV_m0))));
fprintf('    mean1: n=%d, |mu|_max = %.6f, n(|mu|>1) = %d\n', ...
    numel(mu_LEV_m1), max(abs(mu_LEV_m1)), sum(abs(mu_LEV_m1)>1.001));
fprintf('=============================================\n');


%% ========================================================================
%  Local function: read eigenvalues from DMD HDF5
%  ========================================================================
function [mu, mu_proj] = read_eigenvalues(h5file)
    % Read original (pre-projection) eigenvalues if available
    try
        re_orig = h5read(h5file, '/Raw/EigenV_Orig_Real');
        im_orig = h5read(h5file, '/Raw/EigenV_Orig_Imag');
        mu = re_orig + 1i*im_orig;
    catch
        re = h5read(h5file, '/Raw/EigenV_Real');
        im = h5read(h5file, '/Raw/EigenV_Imag');
        mu = re + 1i*im;
    end
    % Also read projected eigenvalues
    try
        re_p = h5read(h5file, '/Raw/EigenV_Real');
        im_p = h5read(h5file, '/Raw/EigenV_Imag');
        mu_proj = re_p + 1i*im_p;
    catch
        mu_proj = mu;
    end
end
