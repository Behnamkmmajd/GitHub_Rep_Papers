clc; clear; close all

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
POD_path  = fullfile(scriptDir, '..', '..', 'Codes', 'POD', 'Results');
POD_path  = char(java.io.File(POD_path).getCanonicalPath());

%% Helper: build POD filename
make_pod = @(obs, s, e, DR, VR, Re, Bo, var) ...
    fullfile(POD_path, sprintf( ...
        'Results_POD_%s_2Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
        obs, s, e, DR, VR, Re, Bo, var));

DMD_path  = fullfile(scriptDir, '..', '..', 'Codes', 'DMD', 'Results');
DMD_path  = char(java.io.File(DMD_path).getCanonicalPath());

%% Helper: build DMD filename
make_dmd = @(obs, s, e, DR, VR, Re, Bo, var) ...
    fullfile(DMD_path, sprintf( ...
        'Results_DMD_%s_2Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
        obs, s, e, DR, VR, Re, Bo, var));

%% Read eigenvalues (stored descending in /Raw/Sigma)
grp = '/Raw/Sigma';

% --- Re = 20, Bo = 20 (300 snapshots) ---
EG_20  = h5read(make_pod('EG',    1,300,10,10,20,20,'mean1'), grp);
EV_20  = h5read(make_pod('EUV',   1,300,10,10,20,20,'mean1'), grp);
LG_20  = h5read(make_pod('LG',    1,300,10,10,20,20,'mean1'), grp);
LV_20  = h5read(make_pod('LUV',   1,300,10,10,20,20,'mean1'), grp);
LEG_20 = h5read(make_pod('LEGi',  1,300,10,10,20,20,'mean1'), grp);
LEV_20 = h5read(make_pod('LEUVi', 1,300,10,10,20,20,'mean1'), grp);
LP_20  = h5read(make_pod('LP',    1,300,10,10,20,20,'der1'),  grp);

% --- Re = 100, Bo = 100 (500 snapshots) ---
EG_100  = h5read(make_pod('EG',    1,500,10,10,100,100,'mean1'), grp);
EV_100  = h5read(make_pod('EUV',   1,500,10,10,100,100,'mean1'), grp);
LG_100  = h5read(make_pod('LG',    1,500,10,10,100,100,'mean1'), grp);
LV_100  = h5read(make_pod('LUV',   1,500,10,10,100,100,'mean1'), grp);
LEG_100 = h5read(make_pod('LEGi',  1,500,10,10,100,100,'mean1'), grp);
LEV_100 = h5read(make_pod('LEUVi', 1,500,10,10,100,100,'mean1'), grp);
LP_100  = h5read(make_pod('LP',    1,500,10,10,100,100,'der1'),  grp);

%% Read reconstruction errors and mean0 Sigma for fair-error profiles
grp_re = '/Reconstruction_error';

% mean1 reconstruction error vectors
re_EG_20   = h5read(make_pod('EG',    1,300,10,10, 20, 20,'mean1'), grp_re);
re_EV_20   = h5read(make_pod('EUV',   1,300,10,10, 20, 20,'mean1'), grp_re);
re_LG_20   = h5read(make_pod('LG',    1,300,10,10, 20, 20,'mean1'), grp_re);
re_LV_20   = h5read(make_pod('LUV',   1,300,10,10, 20, 20,'mean1'), grp_re);
re_LEG_20  = h5read(make_pod('LEGi',  1,300,10,10, 20, 20,'mean1'), grp_re);
re_LEV_20  = h5read(make_pod('LEUVi', 1,300,10,10, 20, 20,'mean1'), grp_re);
re_LP_20   = h5read(make_pod('LP',    1,300,10,10, 20, 20,'der1'),  grp_re);

re_EG_100  = h5read(make_pod('EG',    1,500,10,10,100,100,'mean1'), grp_re);
re_EV_100  = h5read(make_pod('EUV',   1,500,10,10,100,100,'mean1'), grp_re);
re_LG_100  = h5read(make_pod('LG',    1,500,10,10,100,100,'mean1'), grp_re);
re_LV_100  = h5read(make_pod('LUV',   1,500,10,10,100,100,'mean1'), grp_re);
re_LEG_100 = h5read(make_pod('LEGi',  1,500,10,10,100,100,'mean1'), grp_re);
re_LEV_100 = h5read(make_pod('LEUVi', 1,500,10,10,100,100,'mean1'), grp_re);
re_LP_100  = h5read(make_pod('LP',    1,500,10,10,100,100,'der1'),  grp_re);

% mean0 Sigma (for alpha = sum|Sigma_m1| / sum|Sigma_m0|)
s0_EG_20   = h5read(make_pod('EG',    1,300,10,10, 20, 20,'mean0'), grp);
s0_EV_20   = h5read(make_pod('EUV',   1,300,10,10, 20, 20,'mean0'), grp);
s0_LG_20   = h5read(make_pod('LG',    1,300,10,10, 20, 20,'mean0'), grp);
s0_LV_20   = h5read(make_pod('LUV',   1,300,10,10, 20, 20,'mean0'), grp);
s0_LEG_20  = h5read(make_pod('LEGi',  1,300,10,10, 20, 20,'mean0'), grp);
s0_LEV_20  = h5read(make_pod('LEUVi', 1,300,10,10, 20, 20,'mean0'), grp);

s0_EG_100  = h5read(make_pod('EG',    1,500,10,10,100,100,'mean0'), grp);
s0_EV_100  = h5read(make_pod('EUV',   1,500,10,10,100,100,'mean0'), grp);
s0_LG_100  = h5read(make_pod('LG',    1,500,10,10,100,100,'mean0'), grp);
s0_LV_100  = h5read(make_pod('LUV',   1,500,10,10,100,100,'mean0'), grp);
s0_LEG_100 = h5read(make_pod('LEGi',  1,500,10,10,100,100,'mean0'), grp);
s0_LEV_100 = h5read(make_pod('LEUVi', 1,500,10,10,100,100,'mean0'), grp);

%% Read DMD reconstruction errors and SVD_S for fair-error profiles
grp_dmd_rv  = '/Standard_Decomposition/Recons_vec';
grp_dmd_svd = '/Standard_Decomposition/SVD_S';

% --- Re = 20, Bo = 20 — mean1 Recons_vec ---
rv_EG_20_m1  = h5read(make_dmd('EG',    1,300,10,10, 20, 20,'mean1'), grp_dmd_rv);
rv_EV_20_m1  = h5read(make_dmd('EUV',   1,300,10,10, 20, 20,'mean1'), grp_dmd_rv);
rv_LG_20_m1  = h5read(make_dmd('LG',    1,300,10,10, 20, 20,'mean1'), grp_dmd_rv);
rv_LV_20_m1  = h5read(make_dmd('LUV',   1,300,10,10, 20, 20,'mean1'), grp_dmd_rv);
rv_LEG_20_m1 = h5read(make_dmd('LEGi',  1,300,10,10, 20, 20,'mean1'), grp_dmd_rv);
rv_LEV_20_m1 = h5read(make_dmd('LEUVi', 1,300,10,10, 20, 20,'mean1'), grp_dmd_rv);
rv_LP_20     = h5read(make_dmd('LP',    1,300,10,10, 20, 20,'der1'),  grp_dmd_rv);

% --- Re = 100, Bo = 100 — mean1 Recons_vec ---
rv_EG_100_m1  = h5read(make_dmd('EG',    1,500,10,10,100,100,'mean1'), grp_dmd_rv);
rv_EV_100_m1  = h5read(make_dmd('EUV',   1,500,10,10,100,100,'mean1'), grp_dmd_rv);
rv_LG_100_m1  = h5read(make_dmd('LG',    1,500,10,10,100,100,'mean1'), grp_dmd_rv);
rv_LV_100_m1  = h5read(make_dmd('LUV',   1,500,10,10,100,100,'mean1'), grp_dmd_rv);
rv_LEG_100_m1 = h5read(make_dmd('LEGi',  1,500,10,10,100,100,'mean1'), grp_dmd_rv);
rv_LEV_100_m1 = h5read(make_dmd('LEUVi', 1,500,10,10,100,100,'mean1'), grp_dmd_rv);
rv_LP_100     = h5read(make_dmd('LP',    1,500,10,10,100,100,'der1'),  grp_dmd_rv);

% --- Re = 20, Bo = 20 — SVD_S for alpha (||S_m1||^2 / ||S_m0||^2) ---
ss0_EG_20   = h5read(make_dmd('EG',    1,300,10,10, 20, 20,'mean0'), grp_dmd_svd);
ss0_EV_20   = h5read(make_dmd('EUV',   1,300,10,10, 20, 20,'mean0'), grp_dmd_svd);
ss0_LG_20   = h5read(make_dmd('LG',    1,300,10,10, 20, 20,'mean0'), grp_dmd_svd);
ss0_LV_20   = h5read(make_dmd('LUV',   1,300,10,10, 20, 20,'mean0'), grp_dmd_svd);
ss0_LEG_20  = h5read(make_dmd('LEGi',  1,300,10,10, 20, 20,'mean0'), grp_dmd_svd);
ss0_LEV_20  = h5read(make_dmd('LEUVi', 1,300,10,10, 20, 20,'mean0'), grp_dmd_svd);
ss1_EG_20   = h5read(make_dmd('EG',    1,300,10,10, 20, 20,'mean1'), grp_dmd_svd);
ss1_EV_20   = h5read(make_dmd('EUV',   1,300,10,10, 20, 20,'mean1'), grp_dmd_svd);
ss1_LG_20   = h5read(make_dmd('LG',    1,300,10,10, 20, 20,'mean1'), grp_dmd_svd);
ss1_LV_20   = h5read(make_dmd('LUV',   1,300,10,10, 20, 20,'mean1'), grp_dmd_svd);
ss1_LEG_20  = h5read(make_dmd('LEGi',  1,300,10,10, 20, 20,'mean1'), grp_dmd_svd);
ss1_LEV_20  = h5read(make_dmd('LEUVi', 1,300,10,10, 20, 20,'mean1'), grp_dmd_svd);

% --- Re = 100, Bo = 100 — SVD_S for alpha ---
ss0_EG_100  = h5read(make_dmd('EG',    1,500,10,10,100,100,'mean0'), grp_dmd_svd);
ss0_EV_100  = h5read(make_dmd('EUV',   1,500,10,10,100,100,'mean0'), grp_dmd_svd);
ss0_LG_100  = h5read(make_dmd('LG',    1,500,10,10,100,100,'mean0'), grp_dmd_svd);
ss0_LV_100  = h5read(make_dmd('LUV',   1,500,10,10,100,100,'mean0'), grp_dmd_svd);
ss0_LEG_100 = h5read(make_dmd('LEGi',  1,500,10,10,100,100,'mean0'), grp_dmd_svd);
ss0_LEV_100 = h5read(make_dmd('LEUVi', 1,500,10,10,100,100,'mean0'), grp_dmd_svd);
ss1_EG_100  = h5read(make_dmd('EG',    1,500,10,10,100,100,'mean1'), grp_dmd_svd);
ss1_EV_100  = h5read(make_dmd('EUV',   1,500,10,10,100,100,'mean1'), grp_dmd_svd);
ss1_LG_100  = h5read(make_dmd('LG',    1,500,10,10,100,100,'mean1'), grp_dmd_svd);
ss1_LV_100  = h5read(make_dmd('LUV',   1,500,10,10,100,100,'mean1'), grp_dmd_svd);
ss1_LEG_100 = h5read(make_dmd('LEGi',  1,500,10,10,100,100,'mean1'), grp_dmd_svd);
ss1_LEV_100 = h5read(make_dmd('LEUVi', 1,500,10,10,100,100,'mean1'), grp_dmd_svd);

%% =====================  FIGURE STYLE SETTINGS  =====================
%  Modify these values to change appearance across all panels
col_E   = [0.0  0.0  0.0 ];   % Eulerian        = black
col_L   = [0.0  0.0  0.92];   % Lagrangian      = blue
col_ME  = [0.0  0.85 0.2 ];   % Moving-Eulerian = green
col_DNS = [0.8  0.0  0.0 ];   % DNS reference   = red
lw      = 3.0;                 % line width (data)
lw_DNS  = 3.5;                 % line width (DNS)
lw_ax   = 2.5;                 % axis line width
ls_POD  = '-';                 % line style POD
ls_DMD  = '--';                % line style DMD
ls_DNS  = '--';                % line style DNS
fontName    = 'Times New Roman';
fontSize    = 24;
lgdFontSize = 22;
ax_style = {'FontName',fontName,'FontSize',fontSize, ...
            'TickDir','in','LineWidth',lw_ax,'Box','on', ...
            'XMinorTick','on','YMinorTick','on'};
lighter  = @(c) c*0.4 + 0.6;  % pastel shade for Re 20 bars
nmax     = 50;                 % modes shown on x-axis
%  ====================================================================

%% Metric functions (kept for console summary)
cum    = @(e) cumsum(abs(e)) ./ sum(abs(e)) * 100;
rt     = @(e,t) find(cumsum(abs(e))./sum(abs(e)) >= t, 1, 'first');

%% Fair reconstruction error profiles (alpha-corrected; LP der1 is uncorrected)
fe_EG_20   = re_EG_20   * (sum(abs(EG_20))   / sum(abs(s0_EG_20)));
fe_EV_20   = re_EV_20   * (sum(abs(EV_20))   / sum(abs(s0_EV_20)));
fe_LG_20   = re_LG_20   * (sum(abs(LG_20))   / sum(abs(s0_LG_20)));
fe_LV_20   = re_LV_20   * (sum(abs(LV_20))   / sum(abs(s0_LV_20)));
fe_LEG_20  = re_LEG_20  * (sum(abs(LEG_20))  / sum(abs(s0_LEG_20)));
fe_LEV_20  = re_LEV_20  * (sum(abs(LEV_20))  / sum(abs(s0_LEV_20)));
fe_LP_20   = re_LP_20;   % der1: no mean correction

fe_EG_100  = re_EG_100  * (sum(abs(EG_100))  / sum(abs(s0_EG_100)));
fe_EV_100  = re_EV_100  * (sum(abs(EV_100))  / sum(abs(s0_EV_100)));
fe_LG_100  = re_LG_100  * (sum(abs(LG_100))  / sum(abs(s0_LG_100)));
fe_LV_100  = re_LV_100  * (sum(abs(LV_100))  / sum(abs(s0_LV_100)));
fe_LEG_100 = re_LEG_100 * (sum(abs(LEG_100)) / sum(abs(s0_LEG_100)));
fe_LEV_100 = re_LEV_100 * (sum(abs(LEV_100)) / sum(abs(s0_LEV_100)));
fe_LP_100  = re_LP_100;  % der1: no mean correction

%% DMD fair reconstruction error profiles (alpha_DMD = ||SVD_S_m1||^2 / ||SVD_S_m0||^2)
% SVD_S is stored as a full square diagonal matrix; use diag() to get singular values.
% dmd_err_vec extracts error-vs-rank from Recons_vec (see local function at bottom)
dfe_EG_20   = dmd_err_vec(rv_EG_20_m1,  nmax) * (sum(diag(ss1_EG_20).^2)   / sum(diag(ss0_EG_20).^2));
dfe_EV_20   = dmd_err_vec(rv_EV_20_m1,  nmax) * (sum(diag(ss1_EV_20).^2)   / sum(diag(ss0_EV_20).^2));
dfe_LG_20   = dmd_err_vec(rv_LG_20_m1,  nmax) * (sum(diag(ss1_LG_20).^2)   / sum(diag(ss0_LG_20).^2));
dfe_LV_20   = dmd_err_vec(rv_LV_20_m1,  nmax) * (sum(diag(ss1_LV_20).^2)   / sum(diag(ss0_LV_20).^2));
dfe_LEG_20  = dmd_err_vec(rv_LEG_20_m1, nmax) * (sum(diag(ss1_LEG_20).^2)  / sum(diag(ss0_LEG_20).^2));
dfe_LEV_20  = dmd_err_vec(rv_LEV_20_m1, nmax) * (sum(diag(ss1_LEV_20).^2)  / sum(diag(ss0_LEV_20).^2));
dfe_LP_20   = dmd_err_vec(rv_LP_20,     nmax);   % der1: no mean correction

dfe_EG_100  = dmd_err_vec(rv_EG_100_m1,  nmax) * (sum(diag(ss1_EG_100).^2)  / sum(diag(ss0_EG_100).^2));
dfe_EV_100  = dmd_err_vec(rv_EV_100_m1,  nmax) * (sum(diag(ss1_EV_100).^2)  / sum(diag(ss0_EV_100).^2));
dfe_LG_100  = dmd_err_vec(rv_LG_100_m1,  nmax) * (sum(diag(ss1_LG_100).^2)  / sum(diag(ss0_LG_100).^2));
dfe_LV_100  = dmd_err_vec(rv_LV_100_m1,  nmax) * (sum(diag(ss1_LV_100).^2)  / sum(diag(ss0_LV_100).^2));
dfe_LEG_100 = dmd_err_vec(rv_LEG_100_m1, nmax) * (sum(diag(ss1_LEG_100).^2) / sum(diag(ss0_LEG_100).^2));
dfe_LEV_100 = dmd_err_vec(rv_LEV_100_m1, nmax) * (sum(diag(ss1_LEV_100).^2) / sum(diag(ss0_LEV_100).^2));
dfe_LP_100  = dmd_err_vec(rv_LP_100,     nmax);  % der1: no mean correction

%% ========================================================================
%  Figure: 2 rows × 2 columns
%  Row 1: DMD reconstruction error (fair, %)
%  Row 2: POD reconstruction error (fair, %)
%  Col 1: Gas volume fraction (EG / LG / MEG)
%  Col 2: Velocity field + LP (EV / LV / MEV / LP)
%  Re 20 = dashed,  Re 100 = solid
%  ========================================================================
fig = figure('Units','centimeters','Position',[0 0 50 45]);
ax = gobjects(4,1);

% ---- (a) DMD — Gas VF -------------------------------------------------------
ax(1) = subplot(2,2,1); hold on; grid on;
set(gca, 'GridLineWidth', 0.5, 'MinorGridLineWidth', 0.3, 'GridAlpha', 0.3);
% Re 100 — solid
plot(1:nmax, dfe_EG_100(1:nmax),  '-', 'Color',col_E, 'LineWidth',lw);
plot(1:nmax, dfe_LG_100(1:nmax),  '-', 'Color',col_L, 'LineWidth',lw);
plot(1:nmax, dfe_LEG_100(1:nmax), '-', 'Color',col_ME,'LineWidth',lw);
% Re 20 — dashed
plot(1:nmax, dfe_EG_20(1:nmax),  '--', 'Color',col_E, 'LineWidth',lw);
plot(1:nmax, dfe_LG_20(1:nmax),  '--', 'Color',col_L, 'LineWidth',lw);
plot(1:nmax, dfe_LEG_20(1:nmax), '--', 'Color',col_ME,'LineWidth',lw);
set(gca,'YScale','log', ax_style{:}); xlim([1 nmax]);
ylabel('$R_{\rm err,DMD}(r)\%$', 'Interpreter', 'latex', 'FontSize', fontSize);
hE  = plot(NaN,NaN,'-', 'Color',col_E, 'LineWidth',lw);
hL  = plot(NaN,NaN,'-', 'Color',col_L, 'LineWidth',lw);
hME = plot(NaN,NaN,'-', 'Color',col_ME,'LineWidth',lw);
hS  = plot(NaN,NaN,'-k', 'LineWidth',lw);
hD  = plot(NaN,NaN,'--k','LineWidth',lw);
% legend([hE hL hME hS hD], {'E','L','ME','Re = 100','Re = 20'}, ...
%        'FontSize',lgdFontSize,'FontName',fontName,'Location','northeast');

% ---- (b) DMD — Velocity + LP ------------------------------------------------
ax(2) = subplot(2,2,2); hold on; grid on;
set(gca, 'GridLineWidth', 0.5, 'MinorGridLineWidth', 0.3, 'GridAlpha', 0.3);
% Re 100 — solid
plot(1:nmax, dfe_EV_100(1:nmax),  '-',  'Color',col_E, 'LineWidth',lw);
plot(1:nmax, dfe_LV_100(1:nmax),  '-',  'Color',col_L, 'LineWidth',lw);
plot(1:nmax, dfe_LEV_100(1:nmax), '-',  'Color',col_ME,'LineWidth',lw);
plot(1:nmax, dfe_LP_100(1:nmax),  '-',  'Color',col_L, 'LineWidth',lw, ...
     'Marker','o','MarkerSize',8,'MarkerFaceColor',col_L,'MarkerIndices',1:5:nmax);
% Re 20 — dashed
plot(1:nmax, dfe_EV_20(1:nmax),  '--', 'Color',col_E, 'LineWidth',lw);
plot(1:nmax, dfe_LV_20(1:nmax),  '--', 'Color',col_L, 'LineWidth',lw);
plot(1:nmax, dfe_LEV_20(1:nmax), '--', 'Color',col_ME,'LineWidth',lw);
plot(1:nmax, dfe_LP_20(1:nmax),  '--', 'Color',col_L, 'LineWidth',lw, ...
     'Marker','o','MarkerSize',8,'MarkerFaceColor',col_L,'MarkerIndices',1:5:nmax);
set(gca,'YScale','log', ax_style{:}); xlim([1 nmax]);
%ylabel('$R_{\rm err,DMD}(r)$', 'Interpreter', 'latex', 'FontSize', fontSize);
hE  = plot(NaN,NaN,'-',  'Color',col_E, 'LineWidth',lw);
hL  = plot(NaN,NaN,'-',  'Color',col_L, 'LineWidth',lw);
hME = plot(NaN,NaN,'-',  'Color',col_ME,'LineWidth',lw);
hLP = plot(NaN,NaN,'-o', 'Color',col_L, 'LineWidth',lw, ...
           'MarkerSize',8,'MarkerFaceColor',col_L);
hS  = plot(NaN,NaN,'-k', 'LineWidth',lw);
hD  = plot(NaN,NaN,'--k','LineWidth',lw);
hE_leg  = hE;
hL_leg  = hL;
hME_leg = hME;
hLP_leg = hLP;
hS_leg  = hS;
hD_leg  = hD;
% legend([hE hL hME hLP hS hD], {'E','L','ME','LP','Re = 100','Re = 20'}, ...
%        'FontSize',lgdFontSize,'FontName',fontName,'Location','northeast');

% ---- (c) POD — Gas VF -------------------------------------------------------
ax(3) = subplot(2,2,3); hold on; grid on;
set(gca, 'GridLineWidth', 0.5, 'MinorGridLineWidth', 0.3, 'GridAlpha', 0.3);
% Re 100 — solid
plot(1:nmax, fe_EG_100(1:nmax),  '-', 'Color',col_E, 'LineWidth',lw);
plot(1:nmax, fe_LG_100(1:nmax),  '-', 'Color',col_L, 'LineWidth',lw);
plot(1:nmax, fe_LEG_100(1:nmax), '-', 'Color',col_ME,'LineWidth',lw);
% Re 20 — dashed
plot(1:nmax, fe_EG_20(1:nmax),  '--', 'Color',col_E, 'LineWidth',lw);
plot(1:nmax, fe_LG_20(1:nmax),  '--', 'Color',col_L, 'LineWidth',lw);
plot(1:nmax, fe_LEG_20(1:nmax), '--', 'Color',col_ME,'LineWidth',lw);
set(gca,'YScale','log', ax_style{:}); xlim([1 nmax]);
ylabel('$R_{\rm err,POD}(r)\%$', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$Number\ of\ modes,\ r$', 'Interpreter', 'latex', 'FontSize', fontSize);

% ---- (d) POD — Velocity + LP ------------------------------------------------
ax(4) = subplot(2,2,4); hold on; grid on;
set(gca, 'GridLineWidth', 0.5, 'MinorGridLineWidth', 0.3, 'GridAlpha', 0.3);
% Re 100 — solid
plot(1:nmax, fe_EV_100(1:nmax),  '-',  'Color',col_E, 'LineWidth',lw);
plot(1:nmax, fe_LV_100(1:nmax),  '-',  'Color',col_L, 'LineWidth',lw);
plot(1:nmax, fe_LEV_100(1:nmax), '-',  'Color',col_ME,'LineWidth',lw);
plot(1:nmax, fe_LP_100(1:nmax),  '-',  'Color',col_L, 'LineWidth',lw, ...
     'Marker','o','MarkerSize',8,'MarkerFaceColor',col_L,'MarkerIndices',1:5:nmax);
% Re 20 — dashed
plot(1:nmax, fe_EV_20(1:nmax),  '--', 'Color',col_E, 'LineWidth',lw);
plot(1:nmax, fe_LV_20(1:nmax),  '--', 'Color',col_L, 'LineWidth',lw);
plot(1:nmax, fe_LEV_20(1:nmax), '--', 'Color',col_ME,'LineWidth',lw);
plot(1:nmax, fe_LP_20(1:nmax),  '--', 'Color',col_L, 'LineWidth',lw, ...
     'Marker','o','MarkerSize',8,'MarkerFaceColor',col_L,'MarkerIndices',1:5:nmax);
set(gca,'YScale','log', ax_style{:}); xlim([1 nmax]);
%ylabel('$R_{\rm err,POD}(r)$', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$Number\ of\ modes,\ r$', 'Interpreter', 'latex', 'FontSize', fontSize);

% Leave room for a single shared legend below all panels.
for k = 1:numel(ax)
    pos = get(ax(k), 'Position');
    pos(2) = pos(2) + 0.05;
    pos(4) = pos(4) - 0.03;
    set(ax(k), 'Position', pos);
end

% Shared legend for the full figure.
lgd = legend(ax(2), [hE_leg hL_leg hME_leg hLP_leg], ...
    {'Eulerian', 'Lagrangian', 'Moving-Frame', '$\dot{\mathbf{P}}$'}, ...
    'Interpreter', 'latex', ...
    'Orientation', 'horizontal', ...
    'NumColumns', 4, ...
    'FontSize', lgdFontSize, ...
    'Box', 'off');




























%% ========================================================================
%  Console analysis
%  Compares POD vs DMD and Eulerian vs Lagrangian vs Moving-Eulerian
%  for gas volume fraction (GVF) and velocity field (VEL).
%  Metrics evaluated at representative ranks r = 5, 10, 20, 50.
%  Fair reconstruction error (%) — alpha-corrected as in the plots.
%  ========================================================================
probe_ranks = [5 10 20 50];

% ---- pack data for easy iteration ----------------------------------------
%   rows:  {label, POD_re20,  POD_re100,  DMD_re20,   DMD_re100}
configs = { ...
    'GVF  Eulerian (E)',         fe_EG_20,  fe_EG_100,  dfe_EG_20,  dfe_EG_100;  ...
    'GVF  Lagrangian (L)',       fe_LG_20,  fe_LG_100,  dfe_LG_20,  dfe_LG_100;  ...
    'GVF  Moving-Eul (ME)',      fe_LEG_20, fe_LEG_100, dfe_LEG_20, dfe_LEG_100; ...
    'VEL  Eulerian (E)',         fe_EV_20,  fe_EV_100,  dfe_EV_20,  dfe_EV_100;  ...
    'VEL  Lagrangian (L)',       fe_LV_20,  fe_LV_100,  dfe_LV_20,  dfe_LV_100;  ...
    'VEL  Moving-Eul (ME)',      fe_LEV_20, fe_LEV_100, dfe_LEV_20, dfe_LEV_100; ...
    'VEL  Lagr.Pos. (LP)',       fe_LP_20,  fe_LP_100,  dfe_LP_20,  dfe_LP_100;  ...
};
nC = size(configs, 1);

% ---- header ---------------------------------------------------------------
w = 24;   % label column width
sep = repmat('=', 1, w + 4 + numel(probe_ranks)*22 + 2);
fprintf('\n%s\n', sep);
fprintf('RECONSTRUCTION ERROR (fair, %%)  at r = %s\n', ...
    strjoin(arrayfun(@num2str, probe_ranks,'UniformOutput',false), ' / '));
fprintf('%s\n', sep);
hdr = sprintf('%-*s |  Re', w, 'Config');
for rr = probe_ranks
    hdr = [hdr sprintf('    POD r=%-2d   DMD r=%-2d', rr, rr)]; %#ok<AGROW>
end
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-',1,numel(hdr)));

for k = 1:nC
    lbl    = configs{k,1};
    pod20  = configs{k,2};  pod100 = configs{k,3};
    dmd20  = configs{k,4};  dmd100 = configs{k,5};

    for re_tag = {'20','100'}
        tag = re_tag{1};
        if strcmp(tag,'20');  pod = pod20;  dmd = dmd20;
        else;                 pod = pod100; dmd = dmd100; end
        line = sprintf('%-*s | %3s', w, lbl, tag);
        for rr = probe_ranks
            pe = get_err(pod, rr);
            de = get_err(dmd, rr);
            line = [line sprintf('  %9.3f  %9.3f', pe, de)]; %#ok<AGROW>
        end
        fprintf('%s\n', line);
    end
    fprintf('%s\n', repmat('-',1,numel(hdr)));
end

% ---- Best frame per method / field / Re ------------------------------------
fprintf('\nBEST FRAME (lowest fair error at r = 10)\n');
fprintf('%s\n', repmat('-',1,60));
fields   = {'GVF', 'VEL'};
methods  = {'POD', 'DMD'};
res_tags = {'Re=20', 'Re=100'};

% indices into configs{} for GVF rows (1:3) and VEL rows (4:6, excl LP 7)
field_rows = {1:3, 4:6};
frame_lbls = {'E','L','ME'};

for fi = 1:2
    rows = field_rows{fi};
    for mi = 1:2
        col_off = (mi-1)*2 + 2;   % col 2=pod20,3=pod100, 4=dmd20,5=dmd100
        for ri = 1:2
            col = col_off + (ri-1);
            errs = arrayfun(@(k) get_err(configs{k,col}, 10), rows);
            [best_err, bi] = min(errs);
            [worst_err, wi] = max(errs);
            fprintf('  %s | %s | %-6s → best: %s (%.3f%%), worst: %s (%.3f%%)\n', ...
                fields{fi}, methods{mi}, res_tags{ri}, ...
                frame_lbls{bi}, best_err, frame_lbls{wi}, worst_err);
        end
    end
end

% ---- POD vs DMD comparison (average over frames, r=10) --------------------
fprintf('\nPOD vs DMD AVERAGE ERROR (over E/L/ME frames, r = 10)\n');
fprintf('%s\n', repmat('-',1,60));
for fi = 1:2
    rows = field_rows{fi};
    for ri = 1:2
        pod_col = 1 + ri;    % 2 or 3
        dmd_col = 3 + ri;    % 4 or 5
        pod_avg = mean(arrayfun(@(k) get_err(configs{k,pod_col},10), rows));
        dmd_avg = mean(arrayfun(@(k) get_err(configs{k,dmd_col},10), rows));
        winner = 'POD'; if dmd_avg < pod_avg; winner = 'DMD'; end
        fprintf('  %s | %-6s → POD avg=%.3f%%  DMD avg=%.3f%%  → %s lower\n', ...
            fields{fi}, res_tags{ri}, pod_avg, dmd_avg, winner);
    end
end
fprintf('%s\n', repmat('=',1,60));

export_data = struct();
export_data.rank = (1:nmax)';

export_data.pod.gvf.re20.eulerian = fe_EG_20(1:nmax);
export_data.pod.gvf.re20.lagrangian = fe_LG_20(1:nmax);
export_data.pod.gvf.re20.moving_frame = fe_LEG_20(1:nmax);
export_data.pod.gvf.re100.eulerian = fe_EG_100(1:nmax);
export_data.pod.gvf.re100.lagrangian = fe_LG_100(1:nmax);
export_data.pod.gvf.re100.moving_frame = fe_LEG_100(1:nmax);

export_data.pod.velocity.re20.eulerian = fe_EV_20(1:nmax);
export_data.pod.velocity.re20.lagrangian = fe_LV_20(1:nmax);
export_data.pod.velocity.re20.moving_frame = fe_LEV_20(1:nmax);
export_data.pod.velocity.re20.lagrangian_position = fe_LP_20(1:nmax);
export_data.pod.velocity.re100.eulerian = fe_EV_100(1:nmax);
export_data.pod.velocity.re100.lagrangian = fe_LV_100(1:nmax);
export_data.pod.velocity.re100.moving_frame = fe_LEV_100(1:nmax);
export_data.pod.velocity.re100.lagrangian_position = fe_LP_100(1:nmax);

export_data.dmd.gvf.re20.eulerian = dfe_EG_20(1:nmax);
export_data.dmd.gvf.re20.lagrangian = dfe_LG_20(1:nmax);
export_data.dmd.gvf.re20.moving_frame = dfe_LEG_20(1:nmax);
export_data.dmd.gvf.re100.eulerian = dfe_EG_100(1:nmax);
export_data.dmd.gvf.re100.lagrangian = dfe_LG_100(1:nmax);
export_data.dmd.gvf.re100.moving_frame = dfe_LEG_100(1:nmax);

export_data.dmd.velocity.re20.eulerian = dfe_EV_20(1:nmax);
export_data.dmd.velocity.re20.lagrangian = dfe_LV_20(1:nmax);
export_data.dmd.velocity.re20.moving_frame = dfe_LEV_20(1:nmax);
export_data.dmd.velocity.re20.lagrangian_position = dfe_LP_20(1:nmax);
export_data.dmd.velocity.re100.eulerian = dfe_EV_100(1:nmax);
export_data.dmd.velocity.re100.lagrangian = dfe_LV_100(1:nmax);
export_data.dmd.velocity.re100.moving_frame = dfe_LEV_100(1:nmax);
export_data.dmd.velocity.re100.lagrangian_position = dfe_LP_100(1:nmax);

metadata = struct();
metadata.nmax = nmax;
metadata.note = 'Arrays exported from Reconstruction_error.m. These are the plotted fair reconstruction-error profiles used in Figure 4.';
metadata.re20 = struct('re', 20, 'bo', 20, 'start_num', 1, 'end_num', 300);
metadata.re40 = struct('re', 40, 'bo', 40, 'start_num', 1, 'end_num', 400);
metadata.re100 = struct('re', 100, 'bo', 100, 'start_num', 1, 'end_num', 500);

save_reconstruction_error_export(export_data, metadata);

% ---- helper ---------------------------------------------------------------
function v = get_err(vec, r)
    if r <= numel(vec) && ~isnan(vec(r))
        v = vec(r);
    else
        v = NaN;
    end
end

%% ========================================================================
%  Local functions
%  ========================================================================
function Hn = compute_Hnorm(eig)
    e = abs(eig);
    p = e / sum(e);
    p = p(p > 0);
    H = -sum(p .* log(p));
    Hn = H / log(numel(p));
end

function draw_grouped_bars(names, colors, vals_20, vals_100, ~, ax_style, is_int, ylim_range, fontName, lgdFontSize)
    n  = numel(names);
    bw = 0.31;
    for i = 1:n
        % Re 20 bar — same colour, will get hatch overlay
        bh20 = bar(i - bw/2, vals_20(i), bw, ...
            'FaceColor', colors(i,:), 'EdgeColor','k','LineWidth',1.5);
        % Re 100 bar — solid fill
        bar(i + bw/2, vals_100(i), bw, ...
            'FaceColor', colors(i,:), 'EdgeColor','k','LineWidth',1.5);
        % Add diagonal hatch to Re 20 bar
        add_hatch(i - bw/2, bw, vals_20(i));
        % Value labels
        if is_int
            s20  = sprintf('%d',   vals_20(i));
            s100 = sprintf('%d',   vals_100(i));
        else
            s20  = sprintf('%.2f', vals_20(i));
            s100 = sprintf('%.2f', vals_100(i));
        end
        offset = max([vals_20(:); vals_100(:)]) * 0.02;
        text(i - bw/2-0.1, vals_20(i)  + offset, s20,  ...
             'HorizontalAlignment','center','VerticalAlignment','bottom', ...
             'FontSize',lgdFontSize,'FontName',fontName,'FontWeight','bold');
        text(i + bw/2+0.1, vals_100(i) + offset, s100, ...
             'HorizontalAlignment','center','VerticalAlignment','bottom', ...
             'FontSize',lgdFontSize,'FontName',fontName,'FontWeight','bold');
    end
    % Axis formatting
    labs = cellfun(@(s) ['{\it' s '}'], names, 'UniformOutput', false);
    set(gca, 'XTick', 1:n, 'XTickLabel', labs, ax_style{:});
    xlim([0.4  n + 0.6]);
    ylim(ylim_range);
    % Legend: solid grey = Re 100, grey with white lines = Re 20
    % We use two simple patch handles for the legend
    p_solid  = patch(NaN,NaN,[0.5 0.5 0.5],'EdgeColor','k','LineWidth',1.5);
    p_hatch  = patch(NaN,NaN,[0.5 0.5 0.5],'EdgeColor','k','LineWidth',1.5, ...
                     'LineStyle','--');
    % legend([p_hatch p_solid], {'Re = 20 (hatched)','Re = 100 (solid)'}, ...
    %        'FontSize',lgdFontSize,'FontName',fontName,'Location','northeast');
end

function add_hatch(xc, w, h)
    % Draw diagonal white lines inside a bar centred at xc, width w, height h.
    % Uses parametric clipping so it works regardless of axis aspect ratio.
    if isnan(xc) || isnan(h) || h <= 0
        return;
    end
    xl = xc - w/2;   % left edge
    xr = xc + w/2;   % right edge
    yb = 0;           % bottom
    yt = h;           % top

    shift   = w;      % horizontal shift from bottom to top (one bar-width)
    nlines  = 8;
    span    = w + shift;

    for k = 0:nlines
        % Line endpoints before clipping
        x_bot = xl + k * span / nlines;     % bottom point
        x_top = x_bot - shift;              % top point (shifted left)

        % Parametric form: P(t) = (1-t)*(x_bot,yb) + t*(x_top,yt),  t in [0,1]
        dx = x_top - x_bot;
        dy = yt - yb;

        t0 = 0;  t1 = 1;

        % Clip against vertical edges  x = xl  and  x = xr
        if dx ~= 0
            t_xl = (xl - x_bot) / dx;
            t_xr = (xr - x_bot) / dx;
            if dx < 0
                t0 = max(t0, t_xr);
                t1 = min(t1, t_xl);
            else
                t0 = max(t0, t_xl);
                t1 = min(t1, t_xr);
            end
        elseif x_bot < xl || x_bot > xr
            continue;
        end

        if t0 >= t1; continue; end

        px = [x_bot + t0*dx,  x_bot + t1*dx];
        py = [yb    + t0*dy,  yb    + t1*dy];

        plot(px, py, 'w-', 'LineWidth', 1.8, 'HandleVisibility','off');
    end
end

function err_vec = dmd_err_vec(rv, nmax)
    % Extract a reconstruction-error-vs-rank vector from DMD Recons_vec.
    % h5read transposes HDF5 dimensions, so rv has shape (4, N) in MATLAB:
    %   rv(1,:) = gamma_idx,  rv(2,:) = rank,  rv(3,:) = gamma,  rv(4,:) = error%
    % DMD evaluates at discrete gamma thresholds (even rank steps); interpolate
    % linearly to fill integer ranks 1:nmax.
    ranks  = double(rv(2,:));
    errors = double(rv(4,:));
    % Sort by rank ascending so interp1 receives monotonically increasing x
    [r_sorted, si] = sort(ranks);
    e_sorted = errors(si);
    % Interpolate to integer ranks 1:nmax; clamp extrapolation with nearest
    err_vec = interp1(r_sorted, e_sorted, (1:nmax)', 'linear', 'extrap');
end
