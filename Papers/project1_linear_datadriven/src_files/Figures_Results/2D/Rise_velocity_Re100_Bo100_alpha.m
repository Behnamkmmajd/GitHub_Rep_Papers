clc;
clear;
close all

%% =====================  FIGURE STYLE SETTINGS  =====================
%  Modify these values to change the appearance of all plots uniformly.
%  ---- Colors (reference: Figure 3 – spectral compressibility) ----
col_E   = [0.0  0.0  0.0 ];   % Eulerian      — black
col_L   = [0.0  0.0  0.92];   % Lagrangian    — blue
col_ME  = [0.0  0.85 0.2 ];   % Moving-Euler. — green
col_DNS = [0.8  0.0  0.0 ];   % DNS reference — red
%  ---- Line styles ----
ls_POD  = '-';                 % POD  — solid
ls_DMD  = '--';                % DMD  — dashed
ls_DNS  = '--';                % DNS  — dashed
%  ---- Line widths ----
lw      = 3.5;                 % data lines
lw_DNS  = 3.5;                 % DNS (slightly thicker)
lw_ax   = 4.0;                 % axes border / ticks
%  ---- Font ----
fontName    = 'Times New Roman';
fontSize    = 38;
lgdFontSize = 34;
%  ---- Axes preset (pass to set(gca, ...)) ----
ax_style = {'FontName', fontName, 'FontSize', fontSize, ...
            'TickDir', 'in', 'LineWidth', lw_ax, 'Box', 'on', ...
            'XMinorTick', 'on', 'YMinorTick', 'on'};
%  ====================================================================

%% Case Setup
D_bub    = 0.63;       % Bubble diameter (m)
startNum = 1;          % Starting snapshot number
endNum   = 500;        % Ending snapshot number
DR       = 10;         % Density ratio (liquid/gas)
VR       = 10;         % Viscosity ratio (liquid/gas)
Re       = 100;        % Reynolds number
Bo       = 100;        % Bond number
nt       = 498;        % Total number of snapshots
dt       = 0.013107046;    % Time step (s)

%% Scaling Parameters (only for time non-dimensionalization)
v_s_D = sqrt(9.81 * D_bub);   % gravitational velocity scale [m/s]
t_s_D = D_bub / v_s_D;         % gravitational time scale [s]

%% Reconstruction ranks per configuration
%              POD   DMD
rank_EG   = [  10,    10 ];
rank_EUV  = [  10,    10 ];
rank_LEG  = [  10,    10 ];
rank_LEUV = [  10,    10 ];
rank_LG   = [  10,    10 ];
rank_LUV  = [  10,    10 ];
rank_LP   = [  10,    10 ];

%% HDF5 Dataset Group Paths
grp_EG    = sprintf('/xPOD_rank%d/', rank_EG(1));
grp_EUV   = sprintf('/xPOD_rank%d/', rank_EUV(1));
grp_LEG   = sprintf('/xPOD_rank%d/', rank_LEG(1));
grp_LEUV  = sprintf('/xPOD_rank%d/', rank_LEUV(1));
grp_LG    = sprintf('/xPOD_rank%d/', rank_LG(1));
grp_LUV   = sprintf('/xPOD_rank%d/', rank_LUV(1));
grp_LP    = sprintf('/xPOD_rank%d/', rank_LP(1));

grp_dmd_EG    = sprintf('/xDMD_rank%d/', rank_EG(2));
grp_dmd_EUV   = sprintf('/xDMD_rank%d/', rank_EUV(2));
grp_dmd_LEG   = sprintf('/xDMD_rank%d/', rank_LEG(2));
grp_dmd_LEUV  = sprintf('/xDMD_rank%d/', rank_LEUV(2));
grp_dmd_LG    = sprintf('/xDMD_rank%d/', rank_LG(2));
grp_dmd_LUV   = sprintf('/xDMD_rank%d/', rank_LUV(2));
grp_dmd_LP    = sprintf('/xDMD_rank%d/', rank_LP(2));

grp_EXACT = '/x_stacked/';

%% Resolve Results Folders
scriptDir = fileparts(mfilename('fullpath'));

POD_path = fullfile(scriptDir, '..', '..', 'Codes', 'POD', 'Results');
POD_path = char(java.io.File(POD_path).getCanonicalPath());

DMD_path = fullfile(scriptDir, '..', '..', 'Codes', 'DMD', 'Results');
DMD_path = char(java.io.File(DMD_path).getCanonicalPath());

%% Build HDF5 File Paths
%  ---- Variant per configuration & method: 'mean0', 'mean1', or 'der1' ----
%  Customize below to switch any configuration's variant independently
%  for POD and DMD.
%               POD        DMD
var_EG   = {'mean1',   'mean1'};
var_EUV  = {'mean1',   'mean1'};
var_LEG  = {'mean1',   'mean1'};
var_LEUV = {'mean1',   'mean1'};
var_LG   = {'mean1',   'mean1'};
var_LUV  = {'mean1',   'mean1'};
var_LP   = {'der1',   'der1'};

base_fmt = @(method, variant) sprintf( ...
    'Results_%s_2Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
    method, startNum, endNum, DR, VR, Re, Bo, variant);

% Eulerian Gas + Velocity
EG_file_POD  = fullfile(POD_path, base_fmt('POD_EG',    var_EG{1}));
EV_file_POD  = fullfile(POD_path, base_fmt('POD_EUV',   var_EUV{1}));

% Semi-Lagrangian Gas + Velocity
LEG_file_POD = fullfile(POD_path, base_fmt('POD_LEGi',  var_LEG{1}));
LEV_file_POD = fullfile(POD_path, base_fmt('POD_LEUVi', var_LEUV{1}));

% Lagrangian Gas + Velocity + Positions (scattered)
LG_file_POD  = fullfile(POD_path, base_fmt('POD_LG',    var_LG{1}));
LUV_file_POD = fullfile(POD_path, base_fmt('POD_LUV',   var_LUV{1}));
LP_file_POD  = fullfile(POD_path, base_fmt('POD_LP',    var_LP{1}));

% DMD: Eulerian Gas + Velocity
EG_file_DMD  = fullfile(DMD_path, base_fmt('DMD_EG',    var_EG{2}));
EV_file_DMD  = fullfile(DMD_path, base_fmt('DMD_EUV',   var_EUV{2}));

% DMD: Semi-Lagrangian Gas + Velocity
LEG_file_DMD = fullfile(DMD_path, base_fmt('DMD_LEGi',  var_LEG{2}));
LEV_file_DMD = fullfile(DMD_path, base_fmt('DMD_LEUVi', var_LEUV{2}));

% DMD: Lagrangian Gas + Velocity + Positions (scattered)
LG_file_DMD  = fullfile(DMD_path, base_fmt('DMD_LG',    var_LG{2}));
LUV_file_DMD = fullfile(DMD_path, base_fmt('DMD_LUV',   var_LUV{2}));
LP_file_DMD  = fullfile(DMD_path, base_fmt('DMD_LP',    var_LP{2}));

%% Read Grid Dimensions
% Eulerian grid (fixed)
x_E  = h5read(EG_file_POD, '/x');
y_E  = h5read(EG_file_POD, '/y');
nx_E = size(x_E, 1);
ny_E = size(y_E, 1);
[Xg_E, Yg_E] = meshgrid(x_E, y_E);

% SL grid dimensions
x_LE  = h5read(LEG_file_POD, '/x');
y_LE  = h5read(LEG_file_POD, '/y');
nx_LE = size(x_LE, 1);
ny_LE = size(y_LE, 1);
clear x_LE y_LE

% Load semi-Lagrangian grid parameters (replaces per-frame LEP storage)
P_SL = load(fullfile('grid_generator_case_parameters', 'Re100Bo100.mat'));

%% Preallocate
time_nd = (1:nt)' * dt / t_s_D;       % nondimensional time

% Centroid POD: [1] E-POD, [2] ME-POD, [3] L-POD, [4] DNS
% Mean POD:     [5] E-POD, [6] ME-POD, [7] L-POD, [8] DNS
% Centroid DMD: [9] E-DMD, [10] ME-DMD, [11] L-DMD
% Mean DMD:     [12] E-DMD, [13] ME-DMD, [14] L-DMD
yc     = nan(nt, 7);
v_rise = nan(nt, 14);

%% Setup Figures
% --- Two separate figures: (a) centroid, (b) mean ---

% (a) Centroid Rise Velocity
fig_centroid = figure('Units', 'centimeters', 'Position', [0 0 31 29]);
ax1 = axes(fig_centroid);
hold(ax1, 'on');
box(ax1, 'on');

h_EPOD_c  = plot(ax1, NaN, NaN, 'Color', col_E,   'LineStyle', ls_POD, 'LineWidth', lw);
h_LPOD_c  = plot(ax1, NaN, NaN, 'Color', col_L,   'LineStyle', ls_POD, 'LineWidth', lw);
h_SLPOD_c = plot(ax1, NaN, NaN, 'Color', col_ME,  'LineStyle', ls_POD, 'LineWidth', lw);
h_DNS_c   = plot(ax1, NaN, NaN, 'Color', col_DNS,  'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
h_EDMD_c  = plot(ax1, NaN, NaN, 'Color', col_E,   'LineStyle', ls_DMD, 'LineWidth', lw);
h_LDMD_c  = plot(ax1, NaN, NaN, 'Color', col_L,   'LineStyle', ls_DMD, 'LineWidth', lw);
h_SLDMD_c = plot(ax1, NaN, NaN, 'Color', col_ME,  'LineStyle', ls_DMD, 'LineWidth', lw);
uistack(h_DNS_c, 'top');

xlim(ax1, [0 26]);
ylim(ax1, [0.0 0.3])
xlabel(ax1, '{\it\tau}');
%ylabel(ax1, '{\itV}');
set(ax1, ax_style{:});

lgd1 = legend([h_EPOD_c, h_LPOD_c, h_SLPOD_c, h_DNS_c, h_EDMD_c, h_LDMD_c, h_SLDMD_c], ...
    {'E-POD', 'L-POD', 'M-POD', 'DNS', 'E-DMD', 'L-DMD', 'M-DMD'}, ...
    'NumColumns', 4, 'Orientation', 'horizontal', 'FontSize', lgdFontSize, ...
    'FontName', fontName, Location='southeast');
lgd1.ItemTokenSize = [30, 5];

% (b) Mean Rise Velocity
fig_mean = figure('Units', 'centimeters', 'Position', [32 0 31 29]);
ax2 = axes(fig_mean);
hold(ax2, 'on');
box(ax2, 'on');

h_EPOD_m  = plot(ax2, NaN, NaN, 'Color', col_E,   'LineStyle', ls_POD, 'LineWidth', lw);
h_LPOD_m  = plot(ax2, NaN, NaN, 'Color', col_L,   'LineStyle', ls_POD, 'LineWidth', lw);
h_SLPOD_m = plot(ax2, NaN, NaN, 'Color', col_ME,  'LineStyle', ls_POD, 'LineWidth', lw);
h_DNS_m   = plot(ax2, NaN, NaN, 'Color', col_DNS,  'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
h_EDMD_m  = plot(ax2, NaN, NaN, 'Color', col_E,   'LineStyle', ls_DMD, 'LineWidth', lw);
h_LDMD_m  = plot(ax2, NaN, NaN, 'Color', col_L,   'LineStyle', ls_DMD, 'LineWidth', lw);
h_SLDMD_m = plot(ax2, NaN, NaN, 'Color', col_ME,  'LineStyle', ls_DMD, 'LineWidth', lw);
uistack(h_DNS_m, 'top');

xlim(ax2, [0 26]);
ylim(ax2, [0.0 0.3])
xlabel(ax2, '{\it\tau}');
%ylabel(ax2, '{\itV}');
set(ax2, ax_style{:});

lgd2 = legend([h_EPOD_m, h_LPOD_m, h_SLPOD_m, h_DNS_m, h_EDMD_m, h_LDMD_m, h_SLDMD_m], ...
    {'E-POD', 'L-POD', 'M-POD', 'DNS', 'E-DMD', 'L-DMD', 'M-DMD'}, ...
    'NumColumns', 4, 'Orientation', 'horizontal', 'FontSize', lgdFontSize, ...
    'FontName', fontName, Location='southeast');
lgd2.ItemTokenSize = [30, 5];

drawnow

%% Time Loop — Alpha-Weighted Approach (no mask)
for i = startNum:nt

    file_name = sprintf('Stacked_%d', i);
    t_k   = time_nd(i);
    t_km1 = time_nd(max(i - 1, startNum));

    %______________________________________________________________________
    % DNS — Eulerian
    %______________________________________________________________________
    G_dnsE = h5read(EG_file_POD, [grp_EXACT, file_name]);
    V_dnsE = h5read(EV_file_POD, [grp_EXACT, file_name]);
    G_grid = reshape(G_dnsE, [ny_E, nx_E]);
    V_grid = reshape(V_dnsE((ny_E*nx_E+1):end), [ny_E, nx_E]);

    if i == startNum
        [yc(i,4), v_rise(i,4)] = rise_velocity_centroid_alpha(G_grid, Xg_E, Yg_E, NaN, t_k, t_k);
    else
        [yc(i,4), v_rise(i,4)] = rise_velocity_centroid_alpha(G_grid, Xg_E, Yg_E, yc(i-1,4), t_k, t_km1);
    end
    v_rise(i, 8) = rise_velocity_mean_alpha(G_grid, Xg_E, V_grid);

    %______________________________________________________________________
    % E-POD
    %______________________________________________________________________
    G_epod = h5read(EG_file_POD, [grp_EG, file_name]);
    V_epod = h5read(EV_file_POD, [grp_EUV, file_name]);
    G_grid = reshape(G_epod, [ny_E, nx_E]);
    V_grid = reshape(V_epod((ny_E*nx_E+1):end), [ny_E, nx_E]);

    if i == startNum
        [yc(i,1), v_rise(i,1)] = rise_velocity_centroid_alpha(G_grid, Xg_E, Yg_E, NaN, t_k, t_k);
    else
        [yc(i,1), v_rise(i,1)] = rise_velocity_centroid_alpha(G_grid, Xg_E, Yg_E, yc(i-1,1), t_k, t_km1);
    end
    v_rise(i, 5) = rise_velocity_mean_alpha(G_grid, Xg_E, V_grid);

    %______________________________________________________________________
    % E-DMD
    %______________________________________________________________________
    G_edmd = h5read(EG_file_DMD, [grp_dmd_EG, file_name]);
    V_edmd = h5read(EV_file_DMD, [grp_dmd_EUV, file_name]);
    G_grid = reshape(G_edmd, [ny_E, nx_E]);
    V_grid = reshape(V_edmd((ny_E*nx_E+1):end), [ny_E, nx_E]);

    if i == startNum
        [yc(i,5), v_rise(i,9)] = rise_velocity_centroid_alpha(G_grid, Xg_E, Yg_E, NaN, t_k, t_k);
    else
        [yc(i,5), v_rise(i,9)] = rise_velocity_centroid_alpha(G_grid, Xg_E, Yg_E, yc(i-1,5), t_k, t_km1);
    end
    v_rise(i, 12) = rise_velocity_mean_alpha(G_grid, Xg_E, V_grid);

    %______________________________________________________________________
    % SL-POD
    %______________________________________________________________________

    % SL positions (generated via grid_generator_2D instead of LEP file)
    [x_sl, y_sl] = grid_generator_2D(i, P_SL);
    [Xg_LE, Yg_LE] = meshgrid(x_sl, y_sl);

    % SL G and V
    G_slpod = h5read(LEG_file_POD, [grp_LEG, file_name]);
    V_slpod = h5read(LEV_file_POD, [grp_LEUV, file_name]);
    G_grid  = reshape(G_slpod, [ny_LE, nx_LE]);
    V_grid  = reshape(V_slpod((ny_LE*nx_LE)+1:end), [ny_LE, nx_LE]);

    if i == startNum
        [yc(i,2), v_rise(i,2)] = rise_velocity_centroid_alpha(G_grid, Xg_LE, Yg_LE, NaN, t_k, t_k);
    else
        [yc(i,2), v_rise(i,2)] = rise_velocity_centroid_alpha(G_grid, Xg_LE, Yg_LE, yc(i-1,2), t_k, t_km1);
    end
    v_rise(i, 6) = rise_velocity_mean_alpha(G_grid, Xg_LE, V_grid);

    %______________________________________________________________________
    % SL-DMD (ME-DMD)
    %______________________________________________________________________

    % SL positions (generated via grid_generator_2D)
    [x_sl_d, y_sl_d] = grid_generator_2D(i, P_SL);
    [Xg_LE_d, Yg_LE_d] = meshgrid(x_sl_d, y_sl_d);

    % SL G and V (DMD)
    G_sldmd = h5read(LEG_file_DMD, [grp_dmd_LEG, file_name]);
    V_sldmd = h5read(LEV_file_DMD, [grp_dmd_LEUV, file_name]);
    G_grid  = reshape(G_sldmd, [ny_LE, nx_LE]);
    V_grid  = reshape(V_sldmd((ny_LE*nx_LE)+1:end), [ny_LE, nx_LE]);

    if i == startNum
        [yc(i,6), v_rise(i,10)] = rise_velocity_centroid_alpha(G_grid, Xg_LE_d, Yg_LE_d, NaN, t_k, t_k);
    else
        [yc(i,6), v_rise(i,10)] = rise_velocity_centroid_alpha(G_grid, Xg_LE_d, Yg_LE_d, yc(i-1,6), t_k, t_km1);
    end
    v_rise(i, 13) = rise_velocity_mean_alpha(G_grid, Xg_LE_d, V_grid);

    %______________________________________________________________________
    % L-POD (scattered Lagrangian)
    %______________________________________________________________________
    LG_lpod  = h5read(LG_file_POD,  [grp_LG,  file_name]);
    LUV_lpod = h5read(LUV_file_POD, [grp_LUV, file_name]);
    LXY_lpod = h5read(LP_file_POD,  [grp_LP,  file_name]);
    nL = numel(LG_lpod);
    Lx = LXY_lpod(1:nL);  Ly = LXY_lpod(nL+1:2*nL);
    Lv = LUV_lpod(nL+1:2*nL);

    if i == startNum
        [yc(i,3), v_rise(i,3)] = rise_velocity_centroid_alpha(LG_lpod, Lx, Ly, NaN, t_k, t_k);
    else
        [yc(i,3), v_rise(i,3)] = rise_velocity_centroid_alpha(LG_lpod, Lx, Ly, yc(i-1,3), t_k, t_km1);
    end
    v_rise(i, 7) = rise_velocity_mean_alpha(LG_lpod, Lx, Lv);

    %______________________________________________________________________
    % L-DMD (scattered Lagrangian)
    %______________________________________________________________________
    LG_ldmd  = h5read(LG_file_DMD,  [grp_dmd_LG,  file_name]);
    LUV_ldmd = h5read(LUV_file_DMD, [grp_dmd_LUV, file_name]);
    LXY_ldmd = h5read(LP_file_DMD,  [grp_dmd_LP,  file_name]);
    nL_d = numel(LG_ldmd);
    Lx_d = LXY_ldmd(1:nL_d);  Ly_d = LXY_ldmd(nL_d+1:2*nL_d);
    Lv_d = LUV_ldmd(nL_d+1:2*nL_d);

    if i == startNum
        [yc(i,7), v_rise(i,11)] = rise_velocity_centroid_alpha(LG_ldmd, Lx_d, Ly_d, NaN, t_k, t_k);
    else
        [yc(i,7), v_rise(i,11)] = rise_velocity_centroid_alpha(LG_ldmd, Lx_d, Ly_d, yc(i-1,7), t_k, t_km1);
    end
    v_rise(i, 14) = rise_velocity_mean_alpha(LG_ldmd, Lx_d, Lv_d);

    %______________________________________________________________________
    % Update Plots
    %______________________________________________________________________
    t_plot = time_nd(startNum:i);

    % Figure 1: Centroid
    set(h_EPOD_c,  'XData', t_plot, 'YData', v_rise(startNum:i, 1));
    set(h_EDMD_c,  'XData', t_plot, 'YData', v_rise(startNum:i, 9));
    set(h_SLPOD_c, 'XData', t_plot, 'YData', v_rise(startNum:i, 2));
    set(h_SLDMD_c, 'XData', t_plot, 'YData', v_rise(startNum:i, 10));
    set(h_LPOD_c,  'XData', t_plot, 'YData', v_rise(startNum:i, 3));
    set(h_LDMD_c,  'XData', t_plot, 'YData', v_rise(startNum:i, 11));
    set(h_DNS_c,   'XData', t_plot, 'YData', v_rise(startNum:i, 4));

    % Figure 2: Mean
    set(h_EPOD_m,  'XData', t_plot, 'YData', v_rise(startNum:i, 5));
    set(h_EDMD_m,  'XData', t_plot, 'YData', v_rise(startNum:i, 12));
    set(h_SLPOD_m, 'XData', t_plot, 'YData', v_rise(startNum:i, 6));
    set(h_SLDMD_m, 'XData', t_plot, 'YData', v_rise(startNum:i, 13));
    set(h_LPOD_m,  'XData', t_plot, 'YData', v_rise(startNum:i, 7));
    set(h_LDMD_m,  'XData', t_plot, 'YData', v_rise(startNum:i, 14));
    set(h_DNS_m,   'XData', t_plot, 'YData', v_rise(startNum:i, 8));

    drawnow
    disp(file_name)

end

export_data = struct();
export_data.centroid = struct( ...
    'e_pod', v_rise(:,1), ...
    'm_pod', v_rise(:,2), ...
    'l_pod', v_rise(:,3), ...
    'dns', v_rise(:,4), ...
    'e_dmd', v_rise(:,9), ...
    'm_dmd', v_rise(:,10), ...
    'l_dmd', v_rise(:,11));
export_data.mean = struct( ...
    'e_pod', v_rise(:,5), ...
    'm_pod', v_rise(:,6), ...
    'l_pod', v_rise(:,7), ...
    'dns', v_rise(:,8), ...
    'e_dmd', v_rise(:,12), ...
    'm_dmd', v_rise(:,13), ...
    'l_dmd', v_rise(:,14));
export_data.centroid_position = struct( ...
    'e_pod', yc(:,1), ...
    'm_pod', yc(:,2), ...
    'l_pod', yc(:,3), ...
    'dns', yc(:,4), ...
    'e_dmd', yc(:,5), ...
    'm_dmd', yc(:,6), ...
    'l_dmd', yc(:,7));
export_data.raw = struct('v_rise', v_rise, 'yc', yc);

export_metadata = struct( ...
    'start_num', startNum, ...
    'end_num', endNum, ...
    'nt', nt, ...
    'dt', dt, ...
    'D_bub', D_bub, ...
    'v_s_D', v_s_D, ...
    't_s_D', t_s_D, ...
    'rank_EG', rank_EG, ...
    'rank_EUV', rank_EUV, ...
    'rank_LEG', rank_LEG, ...
    'rank_LEUV', rank_LEUV, ...
    'rank_LG', rank_LG, ...
    'rank_LUV', rank_LUV, ...
    'rank_LP', rank_LP);
export_metadata.var_EG = var_EG;
export_metadata.var_EUV = var_EUV;
export_metadata.var_LEG = var_LEG;
export_metadata.var_LEUV = var_LEUV;
export_metadata.var_LG = var_LG;
export_metadata.var_LUV = var_LUV;
export_metadata.var_LP = var_LP;
export_metadata.centroid_definition = 'Rise velocity from time derivative of alpha-weighted bubble centroid position.';
export_metadata.mean_definition = 'Alpha-weighted mean rise velocity computed directly from the velocity field.';

save_rise_velocity_export(Re, Bo, time_nd, export_data, export_metadata);
