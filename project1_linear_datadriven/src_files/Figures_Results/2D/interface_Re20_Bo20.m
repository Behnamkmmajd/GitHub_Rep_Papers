clc;
clear;
close all

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
fontSize    = 42;
%  ---- Axes preset (pass to set(gca, ...)) ----
ax_style = {'FontName', fontName, 'FontSize', fontSize, ...
            'TickDir', 'in', 'LineWidth', lw_ax, 'Box', 'on', ...
            'XMinorTick', 'on', 'YMinorTick', 'on'};
%  ====================================================================

%% Case Setup
startNum = 1;          % Starting snapshot number
endNum   = 300;        % Ending snapshot number
DR       = 10;         % Density ratio (liquid/gas)
VR       = 10;         % Viscosity ratio (liquid/gas)
Re       = 20;         % Reynolds number
Bo       = 20;         % Bond number
nt       = 299;        % Total number of snapshots
threshold = 0.17;      % Gas volume fraction threshold for masking

%% Reconstruction ranks per configuration
%              POD   DMD
% rank_EG  = [  9,    9 ];
% rank_LEG = [  9,    9 ];
% rank_LG  = [  10,    10 ];
% rank_LP  = [  9,    9 ];

%___mean1
rank_EG  = [  10,    10 ];
rank_LEG = [  10,    10 ];
rank_LG  = [  10,    10 ];
rank_LP  = [  9,    9 ];



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

%% Plot scaling
vertical_shift = 0.027;
markerSize = 50;

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
%              POD        DMD
var_EG  = {'mean1',   'mean1'};
var_LEG = {'mean1',   'mean1'};
var_LG  = {'mean1',   'mean1'};
var_LP  = {'der1',   'der1'};

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

%% Read Grids
x_E  = h5read(EG_file_POD, '/x');
y_E  = h5read(EG_file_POD, '/y');
nx_E = size(x_E, 1);
ny_E = size(y_E, 1);
[Xg_E, Yg_E] = meshgrid(x_E, y_E);

x_LE  = h5read(LEG_file_POD, '/x');
y_LE  = h5read(LEG_file_POD, '/y');
nx_LE = size(x_LE, 1);
ny_LE = size(y_LE, 1);
clear x_LE y_LE

% Load semi-Lagrangian grid parameters (replaces per-frame LEP storage)
P_SL = load(fullfile('grid_generator_case_parameters', 'Re20Bo20.mat'));

fig_dmd = figure('Units', 'centimeters', 'Position', [0 0 25 35]);
fig_pod = figure('Units', 'centimeters', 'Position', [34 0 25 35]);

%% Time loop
tIndices = 1:44:nt;
snapshot_data = cell(numel(tIndices), 1);
snapshot_counter = 0;
for i = tIndices
    snapshot_counter = snapshot_counter + 1;

    file_name = sprintf('Stacked_%d', i);
    yOffset = i * vertical_shift;

    %______________________________________________________________________
    % DNS exact (Eulerian frame)
    %______________________________________________________________________
    G_E_exact = h5read(EG_file_POD, [grp_EXACT, file_name]);
    G_E_exact = reshape(G_E_exact, [ny_E, nx_E]);
    G_E_exact(:, 1) = G_E_exact(:, 3);
    G_E_exact(:, 2) = G_E_exact(:, 3);

    %======================================================================
    %  FIGURE 1 — DMD
    %======================================================================
    figure(fig_dmd);

    %______________________________________________________________________
    % Eulerian DMD
    %______________________________________________________________________
    subplot(1, 3, 1);
    EG_dmd = h5read(EG_file_DMD, [grp_dmd_EG, file_name]);

    G_grid = reshape(EG_dmd, [ny_E, nx_E]);
    G_grid(:, 1) = G_grid(:, 2);
    G_edmd_grid = G_grid;

    [row, col] = find(threshold < G_grid);
    Xg_masked = Xg_E(sub2ind(size(Xg_E), row, col));
    Yg_masked = Yg_E(sub2ind(size(Yg_E), row, col));
    G_vals    = G_grid(sub2ind(size(G_grid), row, col));

    scatter( Xg_masked, Yg_masked + yOffset, markerSize, G_vals, 'filled'); hold on
    scatter(-Xg_masked, Yg_masked + yOffset, markerSize, G_vals, 'filled');

    contour( Xg_E, Yg_E + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 3);
    contour(-Xg_E, Yg_E + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 3);

    contour( Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);
    contour(-Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);

    colormap(flipud(gray)); clim([0 1]); daspect([1 1 1])
    xlim([-1 1]); ylim([2 14]); box on
    ylabel('{\itZ}'); xlabel('{\itR}'); title('E-DMD')
    set(gca, ax_style{:});

    %______________________________________________________________________
    % Lagrangian DMD (scattered — no interpolation)
    %______________________________________________________________________
    subplot(1, 3, 2);

    % Positions (intentionally from POD file, matching Rise_velocity)
    Po_L_d = h5read(LP_file_DMD, [grp_dmd_LP, file_name]);
    nL_d   = length(Po_L_d) / 2;
    x_Ld   = Po_L_d(1:nL_d);
    y_Ld   = Po_L_d(nL_d+1:end);

    % Gas (from DMD file)
    G_ldmd = h5read(LG_file_DMD, [grp_dmd_LG, file_name]);
    ldmd_boundary_x = [];
    ldmd_boundary_y = [];

    % Scatter on actual positions (no interpolation)
    mask = G_ldmd > threshold;
    scatter( x_Ld(mask), y_Ld(mask) + yOffset, markerSize, G_ldmd(mask), 'filled'); hold on
    scatter(-x_Ld(mask), y_Ld(mask) + yOffset, markerSize, G_ldmd(mask), 'filled');

    % Approximate 0.5 iso-contour via boundary of alpha > 0.5 points
    bub = G_ldmd > 0.5;
    x_bub = x_Ld(bub);  y_bub = y_Ld(bub);
    if numel(x_bub) > 2
        k = boundary(x_bub, y_bub, 0.5);
        ldmd_boundary_x = x_bub(k);
        ldmd_boundary_y = y_bub(k);
        plot( x_bub(k), y_bub(k) + yOffset, '-y', 'LineWidth', 2);
        plot(-x_bub(k), y_bub(k) + yOffset, '-y', 'LineWidth', 2);
    end

    % DNS exact overlay
    contour( Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);
    contour(-Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);

    colormap(flipud(gray)); clim([0 1]); daspect([1 1 1])
    xlim([-1 1]); ylim([2 14]); box on
    xlabel('{\itR}'); title('L-DMD')
    set(gca, ax_style{:});

    %______________________________________________________________________
    % Semi-Lagrangian DMD (structured deformed grid — pcolor)
    %______________________________________________________________________
    subplot(1, 3, 3);

    % SL positions (generated via grid_generator_2D instead of LEP file)
    [x_sl_d, y_sl_d] = grid_generator_2D(i, P_SL);
    [Xg_LE_d, Yg_LE_d] = meshgrid(x_sl_d, y_sl_d);

    % SL gas (from DMD file)
    G_sldmd = h5read(LEG_file_DMD, [grp_dmd_LEG, file_name]);
    G_grid  = reshape(G_sldmd, [ny_LE, nx_LE]);
    G_grid(:, 1) = G_grid(:, 3);
    G_grid(:, 2) = G_grid(:, 3);
    G_medmd_grid = G_grid;

    % Mask below threshold with NaN for transparency
    G_plot = G_grid;
    G_plot(G_plot < threshold) = NaN;

    % pcolor on deformed structured grid (right half + mirrored left half)
    pcolor( Xg_LE_d, Yg_LE_d + yOffset, G_plot); shading flat; hold on
    pcolor(-Xg_LE_d, Yg_LE_d + yOffset, G_plot); shading flat;

    % Contour directly on deformed structured grid
    contour( Xg_LE_d, Yg_LE_d + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 2);
    contour(-Xg_LE_d, Yg_LE_d + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 2);

    % DNS exact overlay (Eulerian)
    contour( Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);
    contour(-Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);

    colormap(flipud(gray)); clim([0 1]); daspect([1 1 1])
    xlim([-1 1]); ylim([2 14]); box on
    xlabel('{\itR}'); title('M-DMD')
    %cb = colorbar;
    %cb.FontName = fontName; cb.FontSize = fontSize;
    %ylabel(cb, 'Gas volume fraction', 'FontName', fontName, 'FontSize', fontSize, 'Rotation', 90);
    set(gca, ax_style{:});

    %======================================================================
    %  FIGURE 2 — POD
    %======================================================================
    figure(fig_pod);

    %______________________________________________________________________
    % Eulerian POD
    %______________________________________________________________________
    subplot(1, 3, 1);
    EG_pod = h5read(EG_file_POD, [grp_EG, file_name]);

    G_grid = reshape(EG_pod, [ny_E, nx_E]);
    G_grid(:, 1) = G_grid(:, 3);
    G_grid(:, 2) = G_grid(:, 3);
    G_epod_grid = G_grid;

    [row, col] = find(threshold < G_grid);
    Xg_masked = Xg_E(sub2ind(size(Xg_E), row, col));
    Yg_masked = Yg_E(sub2ind(size(Yg_E), row, col));
    G_vals    = G_grid(sub2ind(size(G_grid), row, col));

    scatter( Xg_masked, Yg_masked + yOffset, markerSize, G_vals, 'filled'); hold on
    scatter(-Xg_masked, Yg_masked + yOffset, markerSize, G_vals, 'filled');

    contour( Xg_E, Yg_E + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 2);
    contour(-Xg_E, Yg_E + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 2);

    contour( Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);
    contour(-Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);

    colormap(flipud(gray)); clim([0 1]); daspect([1 1 1])
    xlim([-1 1]); ylim([2 14]); box on
    ylabel('{\itZ}'); xlabel('{\itR}'); title('E-POD')
    set(gca, ax_style{:});

    %______________________________________________________________________
    % Lagrangian POD (scattered — no interpolation)
    %______________________________________________________________________
    subplot(1, 3, 2);

    % Positions
    Po_L = h5read(LP_file_POD, [grp_LP, file_name]);
    nL   = length(Po_L) / 2;
    x_L  = Po_L(1:nL);
    y_L  = Po_L(nL+1:end);

    % Gas
    G_lpod = h5read(LG_file_POD, [grp_LG, file_name]);
    lpod_boundary_x = [];
    lpod_boundary_y = [];

    % Scatter on actual positions (no interpolation)
    mask = G_lpod > threshold;
    scatter( x_L(mask), y_L(mask) + yOffset, markerSize, G_lpod(mask), 'filled'); hold on
    scatter(-x_L(mask), y_L(mask) + yOffset, markerSize, G_lpod(mask), 'filled');

    % Approximate 0.5 iso-contour via boundary
    bub = G_lpod > 0.5;
    x_bub = x_L(bub);  y_bub = y_L(bub);
    if numel(x_bub) > 2
        k = boundary(x_bub, y_bub, 0.5);
        lpod_boundary_x = x_bub(k);
        lpod_boundary_y = y_bub(k);
        plot( x_bub(k), y_bub(k) + yOffset, '-y', 'LineWidth', 2);
        plot(-x_bub(k), y_bub(k) + yOffset, '-y', 'LineWidth', 2);
    end

    % DNS exact overlay
    contour( Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);
    contour(-Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);

    colormap(flipud(gray)); clim([0 1]); daspect([1 1 1])
    xlim([-1 1]); ylim([2 14]); box on
    xlabel('{\itR}'); title('L-POD')
    set(gca, ax_style{:});

    %______________________________________________________________________
    % Semi-Lagrangian POD (structured deformed grid — pcolor)
    %______________________________________________________________________
    subplot(1, 3, 3);

    % SL positions (generated via grid_generator_2D instead of LEP file)
    [x_sl, y_sl] = grid_generator_2D(i, P_SL);
    [Xg_LE, Yg_LE] = meshgrid(x_sl, y_sl);

    % SL gas (from POD file)
    G_slpod = h5read(LEG_file_POD, [grp_LEG, file_name]);
    G_grid  = reshape(G_slpod, [ny_LE, nx_LE]);
    G_grid(:, 1) = G_grid(:, 3);
    G_grid(:, 2) = G_grid(:, 3);
    G_mepod_grid = G_grid;

    % Mask below threshold with NaN for transparency
    G_plot = G_grid;
    G_plot(G_plot < threshold) = NaN;

    % pcolor on deformed structured grid (right half + mirrored left half)
    pcolor( Xg_LE, Yg_LE + yOffset, G_plot); shading flat; hold on
    pcolor(-Xg_LE, Yg_LE + yOffset, G_plot); shading flat;

    % Contour directly on deformed structured grid
    contour( Xg_LE, Yg_LE + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 2);
    contour(-Xg_LE, Yg_LE + yOffset, G_grid, [0.5 0.5], '-y', 'LineWidth', 2);

    % DNS exact overlay (Eulerian)
    contour( Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);
    contour(-Xg_E, Yg_E + yOffset, G_E_exact, [0.5 0.5], '--r', 'LineWidth', 2);

    colormap(flipud(gray)); clim([0 1]); daspect([1 1 1])
    xlim([-1 1]); ylim([2 14]); box on
    xlabel('{\itR}'); title('M-POD')
    %cb = colorbar;
    %cb.FontName = fontName; cb.FontSize = fontSize;
    %ylabel(cb, 'Gas volume fraction', 'FontName', fontName, 'FontSize', fontSize, 'Rotation', 90);
    set(gca, ax_style{:});

    snapshot_data{snapshot_counter} = struct( ...
        'snapshot_index', i, ...
        'file_name', file_name, ...
        'y_offset', yOffset, ...
        'dns', struct('x', x_E, 'y', y_E, 'alpha', G_E_exact), ...
        'e_dmd', struct('x', x_E, 'y', y_E, 'alpha', G_edmd_grid), ...
        'l_dmd', struct('x', x_Ld, 'y', y_Ld, 'alpha', G_ldmd, 'boundary_x', ldmd_boundary_x, 'boundary_y', ldmd_boundary_y), ...
        'm_dmd', struct('x', x_sl_d, 'y', y_sl_d, 'alpha', G_medmd_grid), ...
        'e_pod', struct('x', x_E, 'y', y_E, 'alpha', G_epod_grid), ...
        'l_pod', struct('x', x_L, 'y', y_L, 'alpha', G_lpod, 'boundary_x', lpod_boundary_x, 'boundary_y', lpod_boundary_y), ...
        'm_pod', struct('x', x_sl, 'y', y_sl, 'alpha', G_mepod_grid));

end

export_metadata = struct( ...
    'start_num', startNum, ...
    'end_num', endNum, ...
    'nt', nt, ...
    'threshold', threshold, ...
    'vertical_shift', vertical_shift, ...
    'marker_size', markerSize, ...
    't_indices', tIndices, ...
    'rank_EG', rank_EG, ...
    'rank_LEG', rank_LEG, ...
    'rank_LG', rank_LG, ...
    'rank_LP', rank_LP);
export_metadata.var_EG = var_EG;
export_metadata.var_LEG = var_LEG;
export_metadata.var_LG = var_LG;
export_metadata.var_LP = var_LP;

save_interface_export(Re, Bo, snapshot_data, export_metadata);
