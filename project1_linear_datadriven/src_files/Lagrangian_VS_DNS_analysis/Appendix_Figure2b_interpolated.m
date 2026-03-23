clc;
clear;
close all

%% Case Setup (Do not change these if you want identical results)
D_bub     = 0.63;       % Bubble diameter (m)
startNum  = 2;          % Starting snapshot number
endNum    = 450;        % Ending snapshot number
nt        = endNum;     % Total number of snapshots
dt        = 0.01307;    % Time step (s)
threshold = 0.5;        % Gas volume fraction threshold for masking

time_index = 120;

%% Dataset folders (EDIT THESE)

Eulerian_EG_folder   = 'dataRe100Bo100_EG';
Eulerian_EV_folder   = 'dataRe100Bo100_EUV';

% Multiple Lagrangian cases (add as many as you want)
lagBase = 'dataRe100Bo100';
lagCaseSuffixes = {
    'Vlinear_Glinear'
    'Vlinear_Gcubic'
    'Vlinear_GLG5th'
    'Vcubic_Glinear'
    'Vcubic_Gcubic'
    'Vcubic_GLG5th'
    'VLG5th_Glinear'
    'VLG5th_Gcubic'
    'VLG5th_GLG5th'
};

%% Scaling Parameters
v_s_D = sqrt(9.81 * D_bub);     % Reference velocity (gravitational scale)
t_s_D = D_bub / v_s_D;          % Reference time scale

%% Resolve Results Folders (use script location, not pwd)
scriptDir = pwd;  % Fallback (e.g., running code as selection)
repoRoot = fullfile(scriptDir, '..');
databaseRoot = fullfile(repoRoot, 'database');

Eulerian_G_dir   = fullfile(databaseRoot, Eulerian_EG_folder, 'data');
Eulerian_V_dir   = fullfile(databaseRoot, Eulerian_EV_folder, 'data');

lagCases = struct('label', {}, 'G_dir', {}, 'V_dir', {}, 'P_dir', {});
availableSuffixes = {};
for k = 1:numel(lagCaseSuffixes)
    suffix = lagCaseSuffixes{k};
    gDir = fullfile(databaseRoot, sprintf('%s_LGHa_%s', lagBase, suffix), 'data');
    vDir = fullfile(databaseRoot, sprintf('%s_LUVHa_%s', lagBase, suffix), 'data');
    pDir = fullfile(databaseRoot, sprintf('%s_LPHa_%s', lagBase, suffix), 'data');

    if isfolder(gDir) && isfolder(vDir) && isfolder(pDir)
        lagCases(end+1) = struct('label', suffix, 'G_dir', gDir, 'V_dir', vDir, 'P_dir', pDir); %#ok<SAGROW>
        availableSuffixes{end+1} = suffix; %#ok<SAGROW>
    end
end

lagCaseSuffixes = availableSuffixes;
nLag = numel(lagCases);
if nLag == 0
    error('No matching Lagrangian datasets were found under %s.', databaseRoot);
end

% Eulerian grid size (needed only for reshaping the stacked vectors)
x_E  = readmatrix(fullfile(Eulerian_G_dir, 'x.csv'));
y_E  = readmatrix(fullfile(Eulerian_G_dir, 'y.csv'));
nx_E = numel(x_E);
ny_E = numel(y_E);

% Eulerian mesh
[Xg_E, Yg_E] = meshgrid(x_E, y_E);

%---------------- line locations & tolerances
col_idx = 28;                    % column index for the vertical line
row_idx = 247;                   % row    index for the horizontal line
x_line  = Xg_E(1, col_idx);
y_line  = Yg_E(row_idx, 1);

file_name = sprintf('Stacked_%d', time_index);

%======================================================================
%  EULERIAN FRAME: read G and V, reshape to grid
%======================================================================
EG = readmatrix(fullfile(Eulerian_G_dir, [file_name '.csv']));
EV = readmatrix(fullfile(Eulerian_V_dir, [file_name '.csv']));

G_grid = reshape(EG, [ny_E, nx_E]);
V_grid = reshape(EV(nx_E * ny_E + 1:2 * nx_E * ny_E), [ny_E, nx_E]);

% Eulerian profiles along the two cut-lines
G_hline_E = G_grid(row_idx, :);   % G along horizontal line (vs x_E)
V_hline_E = V_grid(row_idx, :);   % V along horizontal line (vs x_E)
G_vline_E = G_grid(:, col_idx);   % G along vertical   line (vs y_E)
V_vline_E = V_grid(:, col_idx);   % V along vertical   line (vs y_E)

%======================================================================
%  LAGRANGIAN: interpolate scattered data onto the Eulerian grid,
%  then extract profiles along the same cut-lines
%======================================================================
G_hline_L = nan(nx_E, nLag);   % G along horizontal line for each case
V_hline_L = nan(nx_E, nLag);
G_vline_L = nan(ny_E, nLag);   % G along vertical   line for each case
V_vline_L = nan(ny_E, nLag);

for k = 1:nLag
    LG  = readmatrix(fullfile(lagCases(k).G_dir, [file_name '.csv']));
    LV  = readmatrix(fullfile(lagCases(k).V_dir, [file_name '.csv']));
    LXY = readmatrix(fullfile(lagCases(k).P_dir, [file_name '.csv']));

    LG  = LG(:);  LV  = LV(:);  LXY = LXY(:);
    LV  = LV(length(LG)+1:2*length(LG));
    LX  = LXY(1:length(LG));
    LY  = LXY(length(LG)+1:2*length(LG));

    % Interpolate scattered Lagrangian fields onto Eulerian grid
    F_G   = scatteredInterpolant(LX, LY, LG, 'natural', 'none');
    F_V   = scatteredInterpolant(LX, LY, LV, 'natural', 'none');
    G_int = F_G(Xg_E, Yg_E);
    V_int = F_V(Xg_E, Yg_E);

    G_int(isnan(G_int)) = 0;
    V_int(isnan(V_int)) = 0;

    % Extract profiles from the interpolated grid
    G_hline_L(:, k) = G_int(row_idx, :)';
    V_hline_L(:, k) = V_int(row_idx, :)';
    G_vline_L(:, k) = G_int(:, col_idx);
    V_vline_L(:, k) = V_int(:, col_idx);

    % Save first-case interpolated field for subplot 1
    if k == 1
        G_int1 = G_int;
    end
end

%======================================================================
%  FIGURE  –  three subplots  (same styling as Figure1b)
%======================================================================
figure('Position', [100 100 2200 700]);

% Marker styles for the four Lagrangian cases
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

% =====================================================================
%  SUBPLOT 1: Interpolated field of first Lagrangian case, coloured by G
% =====================================================================
ax1 = subplot(1,3,1);
contourf(ax1, Xg_E, Yg_E, G_int1, 'LineStyle', 'none');
hold(ax1,'on')
% Draw the two cut-lines for reference
plot(ax1, [x_line x_line], ax1.YLim, 'b--', 'LineWidth', 1.5);
plot(ax1, ax1.XLim, [y_line y_line], 'r--', 'LineWidth', 1.5);
daspect(ax1, [1 1 1]);
colormap(ax1, flipud(gray));
cb = colorbar(ax1);
ylabel(cb, '\alpha_g', 'FontName', 'Times New Roman', 'FontSize', 28);
xlabel(ax1, 'R');
ylabel(ax1, 'Z');
ylim(ax1,[2.2 3]); 
xlim(ax1,[0 0.6])
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 28, 'LineWidth', 2);
box(ax1, 'on');

% =====================================================================
%  SUBPLOT 2: Gas volume fraction profiles (Eulerian + Lagrangian)
%             Bottom axis = R (red),  Top axis = Z (blue)
% =====================================================================
ax2 = subplot(1,3,2);
hold(ax2, 'on'); 
box(ax2, 'on');

% --- R-direction (bottom axis, red) ---
plot(ax2, x_E, G_hline_E, 'r-', 'LineWidth', 3, 'DisplayName', 'Eulerian (R)');
for k = 1:nLag
    plot(ax2, x_E, G_hline_L(:,k), 'r', 'Marker', markers{k}, ...
        'LineStyle', 'none', 'MarkerSize', 6, 'MarkerFaceColor', 'r', ...
        'DisplayName', sprintf('%s (R)', lagCaseSuffixes{k}));
end
xlabel(ax2, 'R');
ylabel(ax2, '\alpha_g');
xlim(ax2,[0 0.6])
ylim(ax2,[0 1])
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 28, 'LineWidth', 2);

% --- Z-direction (top axis, blue) ---
ax2t = axes('Position', ax2.Position, 'XAxisLocation', 'top', ...
            'YAxisLocation', 'right', 'Color', 'none');
hold(ax2t, 'on');
plot(ax2t, y_E, G_vline_E, 'b-', 'LineWidth', 3, 'DisplayName', 'Eulerian (Z)');
for k = 1:nLag
    plot(ax2t, y_E, G_vline_L(:,k), 'b', 'Marker', markers{k}, ...
        'LineStyle', 'none', 'MarkerSize', 6, ...
        'DisplayName', sprintf('%s (Z)', lagCaseSuffixes{k}));
end
xlabel(ax2t, 'Z');
set(ax2t, 'FontName', 'Times New Roman', 'FontSize', 28, 'LineWidth', 2);
set(ax2t, 'YTick', []);  % hide right y-axis ticks (shared with ax2)
xlim(ax2t,[2.2 3])
ylim(ax2t,[0 1])

% =====================================================================
%  SUBPLOT 3: Velocity profiles (Eulerian + Lagrangian)
%             Bottom axis = R (red),  Top axis = Z (blue)
% =====================================================================
ax3 = subplot(1,3,3);
hold(ax3, 'on'); box(ax3, 'on');

% --- R-direction (bottom axis, red) ---
plot(ax3, x_E, V_hline_E, 'r-', 'LineWidth', 3, 'DisplayName', 'Eulerian (R)');
for k = 1:nLag
    plot(ax3, x_E, V_hline_L(:,k), 'r', 'Marker', markers{k}, ...
        'LineStyle', 'none', 'MarkerSize', 6, 'MarkerFaceColor', 'r', ...
        'DisplayName', sprintf('%s (R)', lagCaseSuffixes{k}));
end
xlabel(ax3, 'R');
ylabel(ax3, 'V_z');
xlim(ax3,[0 0.6])
set(ax3, 'FontName', 'Times New Roman', 'FontSize', 28, 'LineWidth', 2);

% --- Z-direction (top axis, blue) ---
ax3t = axes('Position', ax3.Position, 'XAxisLocation', 'top', ...
            'YAxisLocation', 'right', 'Color', 'none');
hold(ax3t, 'on');
plot(ax3t, y_E, V_vline_E, 'b-', 'LineWidth', 3, 'DisplayName', 'Eulerian (Z)');
for k = 1:nLag
    plot(ax3t, y_E, V_vline_L(:,k), 'b', 'Marker', markers{k}, ...
        'LineStyle', 'none', 'MarkerSize', 6, ...
        'DisplayName', sprintf('%s (Z)', lagCaseSuffixes{k}));
end
xlabel(ax3t, 'Z');
set(ax3t, 'FontName', 'Times New Roman', 'FontSize', 28, 'LineWidth', 2);
set(ax3t, 'YTick', []);  % hide right y-axis ticks (shared with ax3)
xlim(ax3t,[2.2 3])
