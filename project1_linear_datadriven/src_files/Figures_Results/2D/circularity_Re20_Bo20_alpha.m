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
lw_ax   = 4;                 % axes border / ticks
%  ---- Font ----
fontName    = 'Times New Roman';
fontSize    = 44;
lgdFontSize = 40;
%  ---- Axes preset (pass to set(gca, ...)) ----
ax_style = {'FontName', fontName, 'FontSize', fontSize, ...
            'TickDir', 'in', 'LineWidth', lw_ax, 'Box', 'on', ...
            'XMinorTick', 'on', 'YMinorTick', 'on'};
%  ====================================================================

%% Case Setup
D_bub    = 0.3;       % Bubble diameter (m)
startNum = 1;          % Starting snapshot number
endNum   = 300;        % Ending snapshot number
DR       = 10;         % Density ratio (liquid/gas)
VR       = 10;         % Viscosity ratio (liquid/gas)
Re       = 20;        % Reynolds number
Bo       = 20;        % Bond number
nt       = 299;        % Total number of snapshots
dt       = 0.013430661;    % Time step (s)

%% Scaling Parameters (only for time non-dimensionalization)
v_s_D = sqrt(9.81 * D_bub);   % gravitational velocity scale [m/s]
t_s_D = D_bub / v_s_D;         % gravitational time scale [s]

%% Reconstruction ranks per configuration
%              POD   DMD
rank_EG   = [  10,    10 ];
rank_LEG  = [  10,    10 ];
rank_LG   = [  10,    10 ];
rank_LP   = [   9,     9 ];

%% HDF5 Dataset Group Paths
grp_EG    = sprintf('/xPOD_rank%d/', rank_EG(1));
grp_LEG   = sprintf('/xPOD_rank%d/', rank_LEG(1));
grp_LG    = sprintf('/xPOD_rank%d/', rank_LG(1));
grp_LP    = sprintf('/xPOD_rank%d/', rank_LP(1));

grp_dmd_EG    = sprintf('/xDMD_rank%d/', rank_EG(2));
grp_dmd_LEG   = sprintf('/xDMD_rank%d/', rank_LEG(2));
grp_dmd_LG    = sprintf('/xDMD_rank%d/', rank_LG(2));
grp_dmd_LP    = sprintf('/xDMD_rank%d/', rank_LP(2));

grp_EXACT = '/x_stacked/';

%% Resolve Results Folders
scriptDir = fileparts(mfilename('fullpath'));

POD_path = fullfile(scriptDir, '..', '..', 'Codes', 'POD', 'Results');
POD_path = char(java.io.File(POD_path).getCanonicalPath());

DMD_path = fullfile(scriptDir, '..', '..', 'Codes', 'DMD', 'Results');
DMD_path = char(java.io.File(DMD_path).getCanonicalPath());

%% Build HDF5 File Paths
%               POD        DMD
var_EG   = {'mean1',   'mean1'};
var_LEG  = {'mean1',   'mean1'};
var_LG   = {'mean1',   'mean1'};
var_LP   = {'der1',   'der1'};

base_fmt = @(method, variant) sprintf( ...
    'Results_%s_2Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
    method, startNum, endNum, DR, VR, Re, Bo, variant);

% Eulerian Gas
EG_file_POD  = fullfile(POD_path, base_fmt('POD_EG',   var_EG{1}));
EG_file_DMD  = fullfile(DMD_path, base_fmt('DMD_EG',   var_EG{2}));

% Semi-Lagrangian Gas
LEG_file_POD = fullfile(POD_path, base_fmt('POD_LEGi', var_LEG{1}));
LEG_file_DMD = fullfile(DMD_path, base_fmt('DMD_LEGi', var_LEG{2}));

% Lagrangian Gas + Positions (scattered)
LG_file_POD  = fullfile(POD_path, base_fmt('POD_LG',   var_LG{1}));
LP_file_POD  = fullfile(POD_path, base_fmt('POD_LP',   var_LP{1}));
LG_file_DMD  = fullfile(DMD_path, base_fmt('DMD_LG',   var_LG{2}));
LP_file_DMD  = fullfile(DMD_path, base_fmt('DMD_LP',   var_LP{2}));

%% Read Grid Dimensions
% Eulerian grid (fixed)
x_E  = h5read(EG_file_POD, '/x');
y_E  = h5read(EG_file_POD, '/y');
nx_E = size(x_E, 1);
ny_E = size(y_E, 1);

% SL grid dimensions
x_LE  = h5read(LEG_file_POD, '/x');
y_LE  = h5read(LEG_file_POD, '/y');
nx_LE = size(x_LE, 1);
ny_LE = size(y_LE, 1);
clear x_LE y_LE

% Load semi-Lagrangian grid parameters
P_SL = load(fullfile('grid_generator_case_parameters', 'Re20Bo20.mat'));

%% Preallocate
time_nd = (1:nt)' * dt / t_s_D;       % nondimensional time

level_circ = 0.5;  % iso-level used to define bubble interface

% Circularity: [1] E-POD, [2] ME-POD, [3] DNS, [4] E-DMD, [5] ME-DMD, [6] L-POD, [7] L-DMD
circ = nan(nt, 7);

%% Setup Figure
fig = figure('Units', 'centimeters', 'Position', [0 0 31 29]);

ax1 = axes;
hold(ax1, 'on');
box(ax1, 'on');

h_circ_EPOD  = plot(ax1, NaN, NaN, 'Color', col_E,   'LineStyle', ls_POD, 'LineWidth', lw);
h_circ_LPOD  = plot(ax1, NaN, NaN, 'Color', col_L,   'LineStyle', ls_POD, 'LineWidth', lw);
h_circ_MEPOD = plot(ax1, NaN, NaN, 'Color', col_ME,  'LineStyle', ls_POD, 'LineWidth', lw);
h_circ_DNS   = plot(ax1, NaN, NaN, 'Color', col_DNS,  'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
h_circ_EDMD  = plot(ax1, NaN, NaN, 'Color', col_E,   'LineStyle', ls_DMD, 'LineWidth', lw);
h_circ_LDMD  = plot(ax1, NaN, NaN, 'Color', col_L,   'LineStyle', ls_DMD, 'LineWidth', lw);
h_circ_MEDMD = plot(ax1, NaN, NaN, 'Color', col_ME,  'LineStyle', ls_DMD, 'LineWidth', lw);
uistack(h_circ_DNS, 'top');

xlim(ax1, [0 26]);
ylim(ax1, [0 1.05])
xlabel(ax1, '{\it\tau}');
ylabel(ax1, 'Circularity');
set(ax1, ax_style{:});

lgd = legend([h_circ_EPOD, h_circ_LPOD, h_circ_MEPOD, h_circ_DNS, h_circ_EDMD, h_circ_LDMD, h_circ_MEDMD], ...
    {'E-POD', 'L-POD', 'M-POD', 'DNS', 'E-DMD', 'L-DMD', 'M-DMD'}, ...
    'NumColumns', 4, 'Orientation', 'horizontal', 'FontSize', lgdFontSize, ...
    'FontName', fontName, Location='southeast');
lgd.ItemTokenSize = [30, 5];

drawnow

%% Time Loop
for i = startNum:nt

    file_name = sprintf('Stacked_%d', i);

    %______________________________________________________________________
    % DNS — Eulerian
    %______________________________________________________________________
    G_dnsE = h5read(EG_file_POD, [grp_EXACT, file_name]);
    G_grid = reshape(G_dnsE, [ny_E, nx_E]);
    circ(i, 3) = circularity_alpha_grid(G_grid, x_E, y_E, level_circ);

    %______________________________________________________________________
    % E-POD
    %______________________________________________________________________
    G_epod = h5read(EG_file_POD, [grp_EG, file_name]);
    G_grid = reshape(G_epod, [ny_E, nx_E]);
    circ(i, 1) = circularity_alpha_grid(G_grid, x_E, y_E, level_circ);

    %______________________________________________________________________
    % E-DMD
    %______________________________________________________________________
    G_edmd = h5read(EG_file_DMD, [grp_dmd_EG, file_name]);
    G_grid = reshape(G_edmd, [ny_E, nx_E]);
    circ(i, 4) = circularity_alpha_grid(G_grid, x_E, y_E, level_circ);

    %______________________________________________________________________
    % SL-POD (ME-POD)
    %______________________________________________________________________
    [x_sl, y_sl] = grid_generator_2D(i, P_SL);
    G_slpod = h5read(LEG_file_POD, [grp_LEG, file_name]);
    G_grid  = reshape(G_slpod, [ny_LE, nx_LE]);
    circ(i, 2) = circularity_alpha_grid(G_grid, x_sl, y_sl, level_circ);

    %______________________________________________________________________
    % SL-DMD (ME-DMD)
    %______________________________________________________________________
    [x_sl_d, y_sl_d] = grid_generator_2D(i, P_SL);
    G_sldmd = h5read(LEG_file_DMD, [grp_dmd_LEG, file_name]);
    G_grid  = reshape(G_sldmd, [ny_LE, nx_LE]);
    circ(i, 5) = circularity_alpha_grid(G_grid, x_sl_d, y_sl_d, level_circ);

    %______________________________________________________________________
    % L-POD (scattered Lagrangian)
    %______________________________________________________________________
    LG_lpod  = h5read(LG_file_POD, [grp_LG, file_name]);
    LXY_lpod = h5read(LP_file_POD, [grp_LP, file_name]);
    nL = numel(LG_lpod);
    Lx = LXY_lpod(1:nL);  Ly = LXY_lpod(nL+1:2*nL);
    circ(i, 6) = circularity_alpha_scattered(LG_lpod, Lx, Ly, level_circ);

    %______________________________________________________________________
    % L-DMD (scattered Lagrangian)
    %______________________________________________________________________
    LG_ldmd  = h5read(LG_file_DMD, [grp_dmd_LG, file_name]);
    LXY_ldmd = h5read(LP_file_DMD, [grp_dmd_LP, file_name]);
    nL_d = numel(LG_ldmd);
    Lx_d = LXY_ldmd(1:nL_d);  Ly_d = LXY_ldmd(nL_d+1:2*nL_d);
    circ(i, 7) = circularity_alpha_scattered(LG_ldmd, Lx_d, Ly_d, level_circ);

    %______________________________________________________________________
    % Update Plot
    %______________________________________________________________________
    t_plot = time_nd(startNum:i);

    set(h_circ_EPOD,  'XData', t_plot, 'YData', circ(startNum:i, 1));
    set(h_circ_EDMD,  'XData', t_plot, 'YData', circ(startNum:i, 4));
    set(h_circ_MEPOD, 'XData', t_plot, 'YData', circ(startNum:i, 2));
    set(h_circ_MEDMD, 'XData', t_plot, 'YData', circ(startNum:i, 5));
    set(h_circ_LPOD,  'XData', t_plot, 'YData', circ(startNum:i, 6));
    set(h_circ_LDMD,  'XData', t_plot, 'YData', circ(startNum:i, 7));
    set(h_circ_DNS,   'XData', t_plot, 'YData', circ(startNum:i, 3));

    drawnow
    disp(file_name)

end

export_metadata = struct( ...
    'start_num', startNum, ...
    'end_num', endNum, ...
    'nt', nt, ...
    'dt', dt, ...
    'D_bub', D_bub, ...
    'v_s_D', v_s_D, ...
    't_s_D', t_s_D, ...
    'level_circ', level_circ, ...
    'rank_EG', rank_EG, ...
    'rank_LEG', rank_LEG, ...
    'rank_LG', rank_LG, ...
    'rank_LP', rank_LP);
export_metadata.var_EG = var_EG;
export_metadata.var_LEG = var_LEG;
export_metadata.var_LG = var_LG;
export_metadata.var_LP = var_LP;

save_circularity_export(Re, Bo, time_nd, circ, export_metadata);
