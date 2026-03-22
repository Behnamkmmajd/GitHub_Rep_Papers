clc;
clear;
close all

%% =====================  FIGURE STYLE SETTINGS  =====================
%  Modify these values to change appearance across all panels
col_E   = [0.0  0.0  0.0 ];   % Eulerian        = black
col_L   = [0.0  0.0  0.92];   % Lagrangian      = blue
col_ME  = [0.0  0.85 0.2 ];   % Moving-Eulerian = green
col_DNS = [0.8  0.0  0.0 ];   % DNS reference   = red
lw      = 3.5;                 % line width (data)
lw_DNS  = 3.5;                 % line width (DNS)
lw_ax   = 4.0;                 % axis line width
ls_POD  = '-';                 % line style POD
ls_DMD  = '--';                % line style DMD
ls_DNS  = '--';                % line style DNS
fontName    = 'Times New Roman';
fontSize    = 38;
lgdFontSize = 34;
insetFontSize = 18;
ax_style = {'FontName',fontName,'FontSize',fontSize, ...
            'TickDir','in','LineWidth',lw_ax,'Box','on', ...
            'XMinorTick','on','YMinorTick','on'};
%  ====================================================================

%% Case Setup
D_bub    = 0.027;      % Bubble diameter (m)
startNum = 1;          % Starting snapshot number
endNum   = 490;        % Ending snapshot number
DR       = 10;         % Density ratio (liquid/gas)
VR       = 10;         % Viscosity ratio (liquid/gas)
Re       = 168;        % Reynolds number
Bo       = 2;          % Bond number
nt       = 489;        % Total number of snapshots
dt       = 0.005;      % Time step (s)

%% Scaling Parameters (only for time non-dimensionalization)
v_s_D = sqrt(9.81 * D_bub);   % gravitational velocity scale [m/s]
t_s_D = D_bub / v_s_D;         % gravitational time scale [s]

%% Reconstruction ranks
%              POD   DMD
rank_LEG  = [  10,    10 ];
rank_LEUV = [  10,    10 ];

%% HDF5 Dataset Group Paths
grp_LEG   = sprintf('/xPOD_rank%d/', rank_LEG(1));
grp_LEUV  = sprintf('/xPOD_rank%d/', rank_LEUV(1));

grp_dmd_LEG   = sprintf('/xDMD_rank%d/', rank_LEG(2));
grp_dmd_LEUV  = sprintf('/xDMD_rank%d/', rank_LEUV(2));

grp_EXACT = '/x_stacked/';

%% Resolve Results Folders
thisFile = mfilename('fullpath');
if isempty(thisFile)
    scriptDir = pwd;
else
    scriptDir = fileparts(thisFile);
end

POD_path = fullfile(scriptDir, '..', '..', 'Codes', 'POD', 'Results');
POD_path = char(java.io.File(POD_path).getCanonicalPath());

DMD_path = fullfile(scriptDir, '..', '..', 'Codes', 'DMD', 'Results');
DMD_path = char(java.io.File(DMD_path).getCanonicalPath());

%% Build HDF5 File Paths
%  ---- Variant per configuration & method: 'mean0', 'mean1', or 'der1' ----
%               POD        DMD
var_LEG  = {'mean1',   'mean1'};
var_LEUV = {'mean1',   'mean1'};

base_fmt = @(method, variant) sprintf( ...
    'Results_%s_3Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
    method, startNum, endNum, DR, VR, Re, Bo, variant);


% Semi-Lagrangian Gas + Velocity (POD)
LEG_file_POD = fullfile(POD_path, base_fmt('POD_LEGi2',   var_LEG{1}));
LEV_file_POD = fullfile(POD_path, base_fmt('POD_LEUVWi2', var_LEUV{1}));

% Semi-Lagrangian Gas + Velocity (DMD)
LEG_file_DMD = fullfile(DMD_path, base_fmt('DMD_LEGi',   var_LEG{2}));
LEV_file_DMD = fullfile(DMD_path, base_fmt('DMD_LEUVWi', var_LEUV{2}));

%% Read Grid Dimensions
% SL grid dimensions (from HDF5 metadata)
x_LE  = h5read(LEG_file_POD, '/x');
y_LE  = h5read(LEG_file_POD, '/y');
z_LE  = h5read(LEG_file_POD, '/z');
nx_LE = size(x_LE, 1);
ny_LE = size(y_LE, 1);
nz_LE = size(z_LE, 1);
nPts  = nx_LE * ny_LE * nz_LE;
clear x_LE y_LE z_LE

% Load semi-Lagrangian grid parameters
P_SL = load(fullfile('grid_generator_case_parameters', 'Re168Bo2.mat'));

%% Preallocate
time_nd = (1:nt)' * dt / t_s_D;       % nondimensional time

yc_dns   = nan(nt, 1);
yc_pod   = nan(nt, 1);
yc_dmd   = nan(nt, 1);
v_cent_dns = nan(nt, 1);
v_cent_pod = nan(nt, 1);
v_cent_dmd = nan(nt, 1);
v_mean_dns = nan(nt, 1);
v_mean_pod = nan(nt, 1);
v_mean_dmd = nan(nt, 1);

%% Setup Figures
fig_size = [31 29];

fig1 = figure('Units','centimeters','Position',[0 0 31 29]);

% (a) Centroid Rise Velocity
ax1 = axes('Parent', fig1);
hold(ax1, 'on');
box(ax1, 'on');

h_DNS_c   = plot(ax1, NaN, NaN, 'Color', col_DNS,  'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
h_SLPOD_c = plot(ax1, NaN, NaN, 'Color', col_ME, 'LineStyle', ls_POD, 'LineWidth', lw);
h_SLDMD_c = plot(ax1, NaN, NaN, 'Color', col_ME, 'LineStyle', ls_DMD, 'LineWidth', lw);
uistack(h_DNS_c, 'top');

xlim(ax1, [0 50]);
ylim(ax1, [0.0 1.45])
xlabel(ax1, '{\it\tau}');
%ylabel(ax1, '{\itV}_{rise}');
set(ax1, ax_style{:});

legend_handles_1 = [h_SLPOD_c, h_DNS_c, h_SLDMD_c];
legend_labels_1 = {'M-POD', 'DNS', 'M-DMD'};
lgd1 = legend(legend_handles_1, legend_labels_1, ...
    'NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', lgdFontSize, ...
    'FontName', fontName, Location='southeast',NumColumns=1);
lgd1.ItemTokenSize = [30, 5];

fig2 = figure('Units','centimeters','Position',[33 0 31 29]);
% (b) Mean Rise Velocity
ax2 = axes('Parent', fig2);
hold(ax2, 'on');
box(ax2, 'on');

h_DNS_m   = plot(ax2, NaN, NaN, 'Color', col_DNS,  'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
h_SLPOD_m = plot(ax2, NaN, NaN, 'Color', col_ME, 'LineStyle', ls_POD, 'LineWidth', lw);
h_SLDMD_m = plot(ax2, NaN, NaN, 'Color', col_ME, 'LineStyle', ls_DMD, 'LineWidth', lw);
uistack(h_DNS_m, 'top');

xlim(ax2, [0 50]);
ylim(ax2, [0.0 1.45])
xlabel(ax2, '{\it\tau}');
%ylabel(ax2, '{\itV}_{rise}');
set(ax2, ax_style{:});

legend_handles_2 = [h_SLPOD_m, h_DNS_m, h_SLDMD_m];
legend_labels_2 = {'M-POD', 'DNS', 'M-DMD'};
lgd2 = legend(legend_handles_2, legend_labels_2, ...
    'NumColumns', 3, 'Orientation', 'horizontal', 'FontSize', lgdFontSize, ...
    'FontName', fontName, Location='southeast',NumColumns=1);
lgd2.ItemTokenSize = [30, 5];

% Inset for the late-time local mismatch in the average-based rise velocity.
% ax2_inset = axes('Parent', fig2, 'Position', [0.19 0.58 0.30 0.24]);
% hold(ax2_inset, 'on');
% box(ax2_inset, 'on');
% 
% h_DNS_m_inset = plot(ax2_inset, NaN, NaN, 'Color', col_DNS, 'LineStyle', ls_DNS, 'LineWidth', 2.2);
% h_SLPOD_m_inset = plot(ax2_inset, NaN, NaN, 'Color', col_ME, 'LineStyle', ls_POD, 'LineWidth', 2.0);
% h_SLDMD_m_inset = plot(ax2_inset, NaN, NaN, 'Color', col_ME, 'LineStyle', ls_DMD, 'LineWidth', 2.0);
% uistack(h_DNS_m_inset, 'top');
% 
% xlim(ax2_inset, [45.84 46.61]);
% ylim(ax2_inset, [0.90 1.14]);
% set(ax2_inset, 'FontName', fontName, 'FontSize', insetFontSize, ...
%     'TickDir', 'in', 'LineWidth', 2.0, 'Box', 'on', ...
%     'XMinorTick', 'on', 'YMinorTick', 'on');

drawnow

%% Time Loop — Alpha-Weighted Approach
for i = startNum:nt

    file_name = sprintf('Stacked_%d', i);
    t_k   = time_nd(i);
    t_km1 = time_nd(max(i - 1, startNum));

    % SL positions (generated via grid_generator_3D)
    [x_sl, y_sl, z_sl] = grid_generator_3D(i, P_SL);
    [~, Yg_LE, ~] = meshgrid(x_sl, y_sl, z_sl);

    %______________________________________________________________________
    % DNS exact — ME frame
    %______________________________________________________________________
    G_dns = h5read(LEG_file_POD, [grp_EXACT, file_name]);
    V_dns = h5read(LEV_file_POD, [grp_EXACT, file_name]);
    G_grid = permute(reshape(G_dns, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    V_grid = permute(reshape(V_dns(nPts+1:2*nPts), [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);

    if i == startNum
        [yc_dns(i), v_cent_dns(i)] = rise_velocity_centroid_alpha_3D(G_grid, Yg_LE, NaN, t_k, t_k);
    else
        [yc_dns(i), v_cent_dns(i)] = rise_velocity_centroid_alpha_3D(G_grid, Yg_LE, yc_dns(i-1), t_k, t_km1);
    end
    v_mean_dns(i) = rise_velocity_mean_alpha_3D(G_grid, V_grid);

    %______________________________________________________________________
    % ME-POD
    %______________________________________________________________________
    G_slpod = h5read(LEG_file_POD, [grp_LEG, file_name]);
    V_slpod = h5read(LEV_file_POD, [grp_LEUV, file_name]);
    G_grid  = permute(reshape(G_slpod, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    V_grid  = permute(reshape(V_slpod(nPts+1:2*nPts), [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);

    if i == startNum
        [yc_pod(i), v_cent_pod(i)] = rise_velocity_centroid_alpha_3D(G_grid, Yg_LE, NaN, t_k, t_k);
    else
        [yc_pod(i), v_cent_pod(i)] = rise_velocity_centroid_alpha_3D(G_grid, Yg_LE, yc_pod(i-1), t_k, t_km1);
    end
    v_mean_pod(i) = rise_velocity_mean_alpha_3D(G_grid, V_grid);
    if i == 1
        v_mean_pod(i) = 0;
    end
    %______________________________________________________________________
    % ME-DMD
    %______________________________________________________________________
    G_sldmd = h5read(LEG_file_DMD, [grp_dmd_LEG, file_name]);
    V_sldmd = h5read(LEV_file_DMD, [grp_dmd_LEUV, file_name]);
    G_grid  = permute(reshape(G_sldmd, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    V_grid  = permute(reshape(V_sldmd(nPts+1:2*nPts), [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);

    if i == startNum
        [yc_dmd(i), v_cent_dmd(i)] = rise_velocity_centroid_alpha_3D(G_grid, Yg_LE, NaN, t_k, t_k);
    else
        [yc_dmd(i), v_cent_dmd(i)] = rise_velocity_centroid_alpha_3D(G_grid, Yg_LE, yc_dmd(i-1), t_k, t_km1);
    end
    v_mean_dmd(i) = rise_velocity_mean_alpha_3D(G_grid, V_grid);

    %______________________________________________________________________
    % Update Plots
    %______________________________________________________________________
    t_plot = time_nd(startNum:i);

    % Figure 1: Centroid
    set(h_SLPOD_c, 'XData', t_plot, 'YData', v_cent_pod(startNum:i));
    set(h_SLDMD_c, 'XData', t_plot, 'YData', v_cent_dmd(startNum:i));
    set(h_DNS_c,   'XData', t_plot, 'YData', v_cent_dns(startNum:i));

    % Figure 2: Mean
    set(h_SLPOD_m, 'XData', t_plot, 'YData', v_mean_pod(startNum:i));
    set(h_SLDMD_m, 'XData', t_plot, 'YData', v_mean_dmd(startNum:i));
    set(h_DNS_m,   'XData', t_plot, 'YData', v_mean_dns(startNum:i));
    % set(h_SLPOD_m_inset, 'XData', t_plot, 'YData', v_mean_pod(startNum:i));
    % set(h_SLDMD_m_inset, 'XData', t_plot, 'YData', v_mean_dmd(startNum:i));
    % set(h_DNS_m_inset,   'XData', t_plot, 'YData', v_mean_dns(startNum:i));

    drawnow
    disp(file_name)

end

% print(fig1, fullfile(scriptDir, 'Figure7a1.eps'), '-depsc2', '-painters');
% print(fig2, fullfile(scriptDir, 'Figure7b1.eps'), '-depsc2', '-painters');

export_data = struct();
export_data.centroid = struct('dns', v_cent_dns, 'm_pod', v_cent_pod, 'm_dmd', v_cent_dmd);
export_data.mean = struct('dns', v_mean_dns, 'm_pod', v_mean_pod, 'm_dmd', v_mean_dmd);
export_data.centroid_position = struct('dns', yc_dns, 'm_pod', yc_pod, 'm_dmd', yc_dmd);
export_data.raw = struct( ...
    'yc_dns', yc_dns, 'yc_pod', yc_pod, 'yc_dmd', yc_dmd, ...
    'v_cent_dns', v_cent_dns, 'v_cent_pod', v_cent_pod, 'v_cent_dmd', v_cent_dmd, ...
    'v_mean_dns', v_mean_dns, 'v_mean_pod', v_mean_pod, 'v_mean_dmd', v_mean_dmd);

export_metadata = struct( ...
    'start_num', startNum, 'end_num', endNum, 'nt', nt, 'dt', dt, ...
    'D_bub', D_bub, 'v_s_D', v_s_D, 't_s_D', t_s_D, ...
    'rank_LEG', rank_LEG, 'rank_LEUV', rank_LEUV);
export_metadata.var_LEG = var_LEG;
export_metadata.var_LEUV = var_LEUV;
export_metadata.centroid_definition = 'Rise velocity from time derivative of the alpha-weighted 3D centroid coordinate in the rise direction.';
export_metadata.mean_definition = 'Alpha-weighted mean rise velocity in the rise direction computed directly from the 3D velocity field.';

save_rise_velocity_export_3D(Re, Bo, time_nd, export_data, export_metadata);

% fprintf('--- EPS export complete: Figure7a1-Figure7b1 ---\n');
