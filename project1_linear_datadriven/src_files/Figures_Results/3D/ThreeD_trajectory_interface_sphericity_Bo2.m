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
lw_ax   = 3.5;                 % axis line width
ls_POD  = '-';                 % line style POD
ls_DMD  = '--';                % line style DMD
ls_DNS  = '--';                % line style DNS
fontName    = 'Times New Roman';
fontSize    = 28;
lgdFontSize = 26;
ax_style = {'FontName',fontName,'FontSize',fontSize, ...
            'TickDir','in','TickLength',[0.015 0.015],'LineWidth',lw_ax,'Box','on', ...
            'XMinorTick','on','YMinorTick','on'};
ax_stlye2 = {'FontName',fontName,'FontSize',fontSize, ...
            'TickDir','in','TickLength',[0.04 0.04],'LineWidth',lw_ax,'Box','on', ...
            'XMinorTick','on','YMinorTick','on'};

ax_stlye3 = {'FontName',fontName,'FontSize',fontSize, ...
            'TickDir','in','TickLength',[0.04 0.04],'LineWidth',lw_ax,'Box','on', ...
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

%% Scaling Parameters
v_s_D = sqrt(9.81 * D_bub);   % gravitational velocity scale [m/s]
t_s_D = D_bub / v_s_D;         % gravitational time scale [s]

%% ==================  USER CONFIGURATION  ==================
%  Frames at which to overlay DNS isosurfaces on the trajectory plots
snap_iso = 1:60:nt;          % <-- adjust spacing as desired

%  Three frames for the close-up interface comparison (Figure 2)
snap_interface = [1, 22, 474];   % <-- replace with desired frames

%  Iso-level for the bubble interface
isovalue  = 0.5;
contour_min = 0.45;  % contourf only draws levels >= this value
% =============================================================

%% Reconstruction ranks
%              POD   DMD
rank_LEG  = [  10,    10 ];
rank_LEUV = [  10,    10 ];

%% HDF5 Dataset Group Paths
grp_LEG      = sprintf('/xPOD_rank%d/', rank_LEG(1));
grp_LEUV     = sprintf('/xPOD_rank%d/', rank_LEUV(1));
grp_dmd_LEG  = sprintf('/xDMD_rank%d/', rank_LEG(2));
grp_dmd_LEUV = sprintf('/xDMD_rank%d/', rank_LEUV(2));
grp_EXACT    = '/x_stacked/';

%% Resolve Results Folders
thisFile = mfilename('fullpath');
if isempty(thisFile), scriptDir = pwd;
else,                 scriptDir = fileparts(thisFile);
end

POD_path = fullfile(scriptDir, '..', '..', 'Codes', 'POD', 'Results');
POD_path = char(java.io.File(POD_path).getCanonicalPath());
DMD_path = fullfile(scriptDir, '..', '..', 'Codes', 'DMD', 'Results');
DMD_path = char(java.io.File(DMD_path).getCanonicalPath());

%% Build HDF5 File Paths
%               POD        DMD
var_LEG  = {'mean1',   'mean1'};
var_LEUV = {'mean1',   'mean1'};

base_fmt = @(method, variant) sprintf( ...
    'Results_%s_3Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
    method, startNum, endNum, DR, VR, Re, Bo, variant);

LEG_file_POD = fullfile(POD_path, base_fmt('POD_LEGi2',   var_LEG{1}));
LEV_file_POD = fullfile(POD_path, base_fmt('POD_LEUVWi2', var_LEUV{1}));
LEG_file_DMD = fullfile(DMD_path, base_fmt('DMD_LEGi',    var_LEG{2}));
LEV_file_DMD = fullfile(DMD_path, base_fmt('DMD_LEUVWi',  var_LEUV{2}));

%% Read Grid Dimensions
x_LE  = h5read(LEG_file_POD, '/x');
y_LE  = h5read(LEG_file_POD, '/y');
z_LE  = h5read(LEG_file_POD, '/z');
nx_LE = size(x_LE, 1);
ny_LE = size(y_LE, 1);
nz_LE = size(z_LE, 1);
nPts  = nx_LE * ny_LE * nz_LE;
clear x_LE y_LE z_LE

P_SL = load(fullfile('grid_generator_case_parameters', 'Re168Bo2.mat'));

%% Preallocate
time_nd = (1:nt)' * dt / t_s_D;

% 3-D centroid positions: [x, y, z] for each method
centroid_dns = nan(nt, 3);
centroid_pod = nan(nt, 3);
centroid_dmd = nan(nt, 3);

% Sphericity: [1] ME-POD, [2] DNS, [3] ME-DMD
spher = nan(nt, 3);

slice_z = round(nz_LE / 2);   % mid-plane z-index

%% ================================================================
%  PHASE 1 — Main loop: compute centroids + sphericity
%  ================================================================
fprintf('--- Phase 1: computing centroids & sphericity ---\n');
for i = startNum:nt

    file_name = sprintf('Stacked_%d', i);

    [x_sl, y_sl, z_sl] = grid_generator_3D(i, P_SL);
    [Xg, Yg, Zg] = meshgrid(x_sl, y_sl, z_sl);

    % DNS
    G_dns  = h5read(LEG_file_POD, [grp_EXACT, file_name]);
    G_grid = permute(reshape(G_dns, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    centroid_dns(i,:) = centroid_alpha_3D(G_grid, Xg, Yg, Zg);
    spher(i, 2) = sphericity_alpha_grid_3D(G_grid, x_sl, y_sl, z_sl, isovalue);

    % ME-POD
    G_pod  = h5read(LEG_file_POD, [grp_LEG, file_name]);
    G_grid = permute(reshape(G_pod, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    centroid_pod(i,:) = centroid_alpha_3D(G_grid, Xg, Yg, Zg);
    spher(i, 1) = sphericity_alpha_grid_3D(G_grid, x_sl, y_sl, z_sl, isovalue);

    % ME-DMD
    G_dmd  = h5read(LEG_file_DMD, [grp_dmd_LEG, file_name]);
    G_grid = permute(reshape(G_dmd, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    centroid_dmd(i,:) = centroid_alpha_3D(G_grid, Xg, Yg, Zg);
    spher(i, 3) = sphericity_alpha_grid_3D(G_grid, x_sl, y_sl, z_sl, isovalue);

    if mod(i, 50) == 0, fprintf('  frame %d / %d\n', i, nt); end
end
fprintf('--- Phase 1 complete ---\n');

%% ================================================================
%  FIGURE 1a — ZY planar trajectory  (with DNS bubble contours)
%  ================================================================
col_pod = col_ME;             % ME-POD = green
col_dmd = col_ME;             % ME-DMD = green (dashed)
col_dns = col_DNS;            % DNS reference = red

slice_x = round(nx_LE / 2);   % mid-plane x-index

fig1a = figure('Units','centimeters','Position',[0 0 9 29]);
ax_zy = axes(fig1a);
hold(ax_zy, 'on');  
box(ax_zy, 'on');

% DNS bubble contours at selected frames
for k = snap_iso(:)'
    file_name = sprintf('Stacked_%d', k);
    [x_sl, y_sl, z_sl] = grid_generator_3D(k, P_SL);
    [~, Yg, Zg] = meshgrid(x_sl, y_sl, z_sl);

    G_dns  = h5read(LEG_file_POD, [grp_EXACT, file_name]);
    G_grid = permute(reshape(G_dns, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    G_slice = squeeze(G_grid(:, slice_x, :));

    Z_slice = squeeze(Zg(:, slice_x, :));
    Y_slice = squeeze(Yg(:, slice_x, :));

    contourf(ax_zy, Z_slice, Y_slice, G_slice, linspace(contour_min, 1, 20), 'LineStyle', 'none');
end

% Trajectory lines
h_pod_zy = plot(ax_zy, centroid_pod(:,3), centroid_pod(:,2), ...
    'Color', col_pod, 'LineStyle', ls_POD, 'LineWidth', lw);
h_dmd_zy = plot(ax_zy, centroid_dmd(:,3), centroid_dmd(:,2), ...
    'Color', col_dmd, 'LineStyle', ls_DMD, 'LineWidth', lw);
h_dns_zy = plot(ax_zy, centroid_dns(:,3), centroid_dns(:,2), ...
    'Color', col_dns, 'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
uistack(h_dns_zy, 'top');

colormap(ax_zy, flipud(gray));
clim(ax_zy, [0 1]);
daspect(ax_zy, [1 1 1]);
xlim([2.5 7.5])
ylim([25 50])
%xlabel(ax_zy, '{\itZ} / {\itD}');
%ylabel(ax_zy, '{\itY} / {\itD}');
set(ax_zy, ax_style{:});
%legend([h_pod_zy, h_dns_zy, h_dmd_zy], {'ME-POD', 'DNS', 'ME-DMD'}, 'FontSize', 22, 'FontName', 'Times New Roman', 'Location', 'best');
drawnow

%% ================================================================
%  FIGURE 1b — 3-D trajectory  (with DNS isosurfaces)
%  ================================================================
fig1b = figure('Units','centimeters','Position',[10 0 9 29]);
ax_3d = axes(fig1b);
hold(ax_3d, 'on');  box(ax_3d, 'on');

% DNS isosurfaces at selected frames
for k = snap_iso(:)'
    file_name = sprintf('Stacked_%d', k);
    [x_sl, y_sl, z_sl] = grid_generator_3D(k, P_SL);
    [Xg, Yg, Zg] = meshgrid(x_sl, y_sl, z_sl);

    G_dns  = h5read(LEG_file_POD, [grp_EXACT, file_name]);
    G_grid = permute(reshape(G_dns, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);

    fv = isosurface(Xg, Yg, Zg, G_grid, isovalue);
    if ~isempty(fv.vertices)
        fv.vertices = fv.vertices(:,[1 3 2]);   % swap Y<->Z so rise axis is vertical
        p = patch(ax_3d, fv);
        p.FaceColor = 'k';
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.85;
    end
end

% 3-D trajectory lines
h_pod_3d = plot3(ax_3d, centroid_pod(:,1), centroid_pod(:,3), centroid_pod(:,2), ...
    'Color', col_pod, 'LineStyle', ls_POD, 'LineWidth', lw);
h_dmd_3d = plot3(ax_3d, centroid_dmd(:,1), centroid_dmd(:,3), centroid_dmd(:,2), ...
    'Color', col_dmd, 'LineStyle', ls_DMD, 'LineWidth', lw);
h_dns_3d = plot3(ax_3d, centroid_dns(:,1), centroid_dns(:,3), centroid_dns(:,2), ...
    'Color', col_dns, 'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
uistack(h_dns_3d, 'top');

camlight;
lighting gouraud;
zlim([0 50])
xlim([0 8])
ylim([0 8])
daspect(ax_3d, [1 1 1]);
grid(ax_3d, 'on');
%xlabel(ax_3d, '{\itX} / {\itD}');
%ylabel(ax_3d, '{\itY} / {\itD}');
%zlabel(ax_3d, '{\itZ} / {\itD}');
set(ax_3d, ax_style{:}, 'ZMinorTick', 'on');
%legend([h_pod_3d, h_dns_3d, h_dmd_3d], {'ME-POD', 'DNS', 'ME-DMD'},'FontSize', 22, 'FontName', 'Times New Roman', 'Location', 'best');
view(ax_3d, [-37.5+90 30]);
drawnow

%% ================================================================
%  FIGURE 2 — Interface comparison at three time instances
%  ================================================================
nSnap = numel(snap_interface);
fig2 = figure('Units','centimeters','Position',[20 0 9.66 9.66*nSnap]);

for s = 1:nSnap
    k = snap_interface(s);
    file_name = sprintf('Stacked_%d', k);
    tau_k = time_nd(k);

    [x_sl, y_sl, z_sl] = grid_generator_3D(k, P_SL);
    [Xg, Yg, ~] = meshgrid(x_sl, y_sl, z_sl);
    X_slice = Xg(:,:, slice_z);
    Y_slice = Yg(:,:, slice_z);

    ax = subplot(nSnap, 1, nSnap - s + 1);   % bottom = earliest
    hold(ax, 'on');  box(ax, 'on');

    % --- DNS filled contour ---
    G_dns  = h5read(LEG_file_POD, [grp_EXACT, file_name]);
    G_dns3 = permute(reshape(G_dns, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    G_DNs    = G_dns3(:,:, slice_z);
    contourf(ax, X_slice, Y_slice, G_DNs, linspace(contour_min, 1, 20), 'LineStyle', 'none');

    % --- ME-POD interface contour ---
    G_pod  = h5read(LEG_file_POD, [grp_LEG, file_name]);
    G_pod3 = permute(reshape(G_pod, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    G_s    = G_pod3(:,:, slice_z);
    contour(ax, X_slice, Y_slice, G_s, [isovalue isovalue], ...
        'Color', col_pod, 'LineStyle', '-', 'LineWidth', 4);

    % --- ME-DMD interface contour ---
    G_dmd  = h5read(LEG_file_DMD, [grp_dmd_LEG, file_name]);
    G_dmd3 = permute(reshape(G_dmd, [nx_LE, ny_LE, nz_LE]), [2, 1, 3]);
    G_s    = G_dmd3(:,:, slice_z);
    contour(ax, X_slice, Y_slice, G_s, [isovalue isovalue], ...
        'Color', col_dmd, 'LineStyle', '--', 'LineWidth', 4);

    % --- DNS interface contour ---
    contour(ax, X_slice, Y_slice, G_DNs, [isovalue isovalue], ...
        'Color', col_dns, 'LineStyle', '--', 'LineWidth', 4);

    colormap(ax, flipud(gray));
    clim(ax, [0 1]);
    daspect(ax, [1 1 1]);
    text(ax, 0.05, 0.07, sprintf('\\tau = %.2f', tau_k), ...
        'Units', 'normalized', 'FontName', fontName, 'FontSize', 22,'FontWeight','bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'Color', 'k', 'BackgroundColor', 'w', 'Margin', 2);
    % ylabel(ax, '{\itY} / {\itD}');
    %title(ax, sprintf('\\tau = %.2f', tau_k), 'FontName', 'Times New Roman', 'FontSize', 28);
    set(ax, ax_stlye2{:});

    % if s == 1   % add legend to last-drawn (topmost visual) panel
    % 
    %     xlabel(ax, '{\itX} / {\itD}');
    % 
    % end
end

% Shared colorbar
%cb = colorbar(ax, 'eastoutside');
%ylabel(cb, 'Volume fraction', 'FontSize', 26, 'FontName', 'Times New Roman');
drawnow

% ================================================================
%  FIGURE 3 — Sphericity vs time
%  ================================================================
fig3 = figure('Units','centimeters','Position',[31 0 31 29]);
ax_sp = axes(fig3);
hold(ax_sp, 'on');  box(ax_sp, 'on');

h_sp_pod = plot(ax_sp, time_nd, spher(:,1), 'Color', col_pod, 'LineStyle', ls_POD, 'LineWidth', lw);
h_sp_dns = plot(ax_sp, time_nd, spher(:,2), 'Color', col_dns, 'LineStyle', ls_DNS, 'LineWidth', lw_DNS);
h_sp_dmd = plot(ax_sp, time_nd, spher(:,3), 'Color', col_dmd, 'LineStyle', ls_DMD, 'LineWidth', lw);
uistack(h_sp_dns, 'top');

ax_sp_inset = axes('Parent', fig3, 'Position', [0.28 0.64 0.30 0.24]);
hold(ax_sp_inset, 'on');
box(ax_sp_inset, 'on');

h_sp_dns_inset = plot(ax_sp_inset, time_nd, spher(:,2), 'Color', col_DNS, 'LineStyle', ls_DNS, 'LineWidth', 3.5);
h_sp_pod_inset = plot(ax_sp_inset, time_nd, spher(:,1), 'Color', col_ME, 'LineStyle', ls_POD, 'LineWidth', 3.5);
h_sp_dmd_inset = plot(ax_sp_inset, time_nd, spher(:,3), 'Color', col_ME, 'LineStyle', ls_DMD, 'LineWidth', 3.5);
uistack(h_sp_dns_inset, 'top');

xlim(ax_sp_inset, [0 4]);
ylim(ax_sp_inset, [0.94 1]);
set(ax_sp_inset, 'FontName', fontName, 'FontSize', 28, ...
    'TickDir', 'in', 'LineWidth', 3.0, 'Box', 'on', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', [0.03 0.03]);

xlim(ax_sp, [0 50]);
ylim(ax_sp, [0.94 1.005]);
xlabel(ax_sp, '{\it\tau}');
%ylabel(ax_sp, 'Sphericity');
set(ax_sp, ax_style{:});

annotation(gcf,'rectangle',...
    [0.130244249726177 0.146778042959427 0.0570503833515882 0.714797136038189],...
    'LineWidth',3,...
    'LineStyle','--');


legend([h_sp_pod, h_sp_dns, h_sp_dmd], {'M-POD', 'DNS', 'M-DMD'}, 'NumColumns', 1, 'Orientation', 'horizontal', 'FontSize', 27,'FontName', 'Times New Roman','Location','northeast');
drawnow
% 
fprintf('--- All figures complete ---\n');


print(fig1a, fullfile(scriptDir, 'Figure7a1.eps'), '-depsc2', '-vector');
print(fig1b, fullfile(scriptDir, 'Figure7b1.eps'), '-depsc2', '-vector');
print(fig2,  fullfile(scriptDir, 'Figure7c1.eps'), '-depsc2', '-vector');
print(fig3,  fullfile(scriptDir, 'Figure7d1.eps'), '-depsc2', '-vector');

export_data = struct();
export_data.trajectory = struct('dns', centroid_dns, 'm_pod', centroid_pod, 'm_dmd', centroid_dmd);
export_data.sphericity = struct('dns', spher(:, 2), 'm_pod', spher(:, 1), 'm_dmd', spher(:, 3));
export_data.snapshots = struct('snap_iso', snap_iso, 'snap_interface', snap_interface, 'slice_z', slice_z, 'isovalue', isovalue);

export_metadata = struct( ...
    'start_num', startNum, 'end_num', endNum, 'nt', nt, 'dt', dt, ...
    'D_bub', D_bub, 'v_s_D', v_s_D, 't_s_D', t_s_D, ...
    'rank_LEG', rank_LEG, 'rank_LEUV', rank_LEUV);
export_metadata.var_LEG = var_LEG;
export_metadata.var_LEUV = var_LEUV;
export_metadata.trajectory_definition = '3D centroid trajectory extracted from the alpha-weighted centroid position in the moving-Eulerian frame.';
export_metadata.sphericity_definition = 'Sphericity obtained from the alpha=0.5 isosurface of the reconstructed 3D bubble.';

save_trajectory_sphericity_export_3D(Re, Bo, time_nd, export_data, export_metadata);

fprintf('--- EPS export complete: Figure6a1-Figure6d1 ---\n');
