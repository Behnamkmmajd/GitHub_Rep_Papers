clc;
clear;
close all

%% Case Setup (Do not change these if you want identical results)
D_bub     = 0.63;       % Bubble diameter (m)
startNum  = 1;          % Starting snapshot number
endNum    = 450;        % Ending snapshot number
nt        = endNum;     % Total number of snapshots
dt_sim    = 0.05;       % Snapshot spacing in the simulation non-dimensional time
dt_phys   = 0.01307;    % Physical time step between snapshots (s)

%% Dataset folders (EDIT THESE)

Eulerian_EG_folder   = 'dataRe100Bo100_EG';
Eulerian_EV_folder   = 'dataRe100Bo100_EUV';

% Multiple Lagrangian cases (add as many as you want)
lagBase = 'dataRe100Bo100';
lagCaseSuffixes = {
    'Vcubic_Gcubic'
};

%% Gravity-based normalization built from the physical time scale
v_s_D      = sqrt(9.81 * D_bub);     % Reference velocity (gravitational scale)
t_s_D      = D_bub / v_s_D;          % Reference time scale
t_unit_phys = dt_phys / dt_sim;      % Physical time represented by one simulation time unit
v_scale     = D_bub / (t_unit_phys * v_s_D);

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

ns=ny_E*nx_E;
% Eulerian mesh
[Xg_E, Yg_E] = meshgrid(x_E, y_E);

%% Preallocate time + outputs
time_sim  = (1:nt)' * dt_sim;
time_phys = time_sim * t_unit_phys;
time      = time_phys / t_s_D;
circ_level = 0.5;

% Eulerian:
v_E_field    = nan(nt, 1);   % velocity-field average
yc_E         = nan(nt, 1);   % centroid position
v_E_centroid = nan(nt, 1);   % centroid-based rise velocity
circ_E       = nan(nt, 1);   % circularity profile

% Lagrangian: centroid + volume-averaged, one curve per case
v_L_field    = nan(nt, nLag);   % volume-averaged rise velocity
yc_L         = nan(nt, nLag);
v_L_centroid = nan(nt, nLag);
circ_L       = nan(nt, nLag);

%% Time Loop: Calculate Rise Velocity for Each Snapshot
for i = startNum:endNum

    file_name = sprintf('Stacked_%d', i);

    %======================================================================
    %  EULERIAN FRAME: read G and V, reshape to grid
    %======================================================================
    EG = readmatrix(fullfile(Eulerian_G_dir, [file_name '.csv']));
    EV = readmatrix(fullfile(Eulerian_V_dir, [file_name '.csv']));

    G_grid = reshape(EG, [ny_E, nx_E]);
    V_grid = reshape(EV(ns+1:2*ns), [ny_E, nx_E]) * v_scale;

    %----------------------------------------------------------------------
    %  Eulerian rise velocity  –  velocity-field average (rise_velocity_mean)
    %----------------------------------------------------------------------
    v_E_field(i) = rise_velocity_mean_alpha(G_grid, Xg_E, V_grid);
    circ_E(i) = circularity_alpha_grid(G_grid, x_E, y_E, circ_level);

    %----------------------------------------------------------------------
    %  Eulerian rise velocity  –  centroid motion dyc/dt (rise_velocity_centroid)
    %----------------------------------------------------------------------
    if i == startNum
        [yc_E(i), v_E_centroid(i)] = rise_velocity_centroid_alpha( ...
            G_grid, Xg_E, Yg_E, NaN, time(i), time(i));
    else
        [yc_E(i), v_E_centroid(i)] = rise_velocity_centroid_alpha( ...
            G_grid, Xg_E, Yg_E, yc_E(i-1), time(i), time(i-1));
    end

    %======================================================================
    %  LAGRANGIAN FRAMES: scattered G, V, positions passed directly
    %  (no interpolation)
    %======================================================================
    for k = 1:nLag
        LG  = readmatrix(fullfile(lagCases(k).G_dir, [file_name '.csv']));
        LV  = readmatrix(fullfile(lagCases(k).V_dir, [file_name '.csv']));
        LXY = readmatrix(fullfile(lagCases(k).P_dir, [file_name '.csv']));

        LG  = LG(:);
        LV  = LV(:);
        LXY = LXY(:);
        LV  = LV(length(LG)+1:2*length(LG)) * v_scale;
        LX  = LXY(1:length(LG));
        LY  = LXY(length(LG)+1:2*length(LG));

        %------------------------------------------------------------------
        %  Lagrangian rise velocity  –  volume-averaged (rise_velocity_mean)
        %  (ones for XG since scattered points have no radial cell width)
        %------------------------------------------------------------------
        v_L_field(i, k) = rise_velocity_mean_alpha(LG, LX, LV);
        circ_L(i, k) = circularity_alpha_scattered(LG, LX, LY, circ_level);

        %------------------------------------------------------------------
        %  Lagrangian rise velocity  –  centroid motion dyc/dt
        %------------------------------------------------------------------
        if i == startNum
            [yc_L(i, k), v_L_centroid(i, k)] = rise_velocity_centroid_alpha( ...
                LG, LX, LY, NaN, time(i), time(i));
        else
            [yc_L(i, k), v_L_centroid(i, k)] = rise_velocity_centroid_alpha( ...
                LG, LX, LY, yc_L(i-1, k), time(i), time(i-1));
        end
    end

    display(file_name)

end
%%


load("cache_Fig1a.mat")

%Plot Results
t_plot = time(startNum:endNum);
idx = startNum:endNum;

figure('Position', [0 0 1800 900]);

% ---------------- Subplot 1: Rise velocity profile (case 1 only) --------
ax1 = subplot(1,2,1);
hold(ax1,'on'); box(ax1,'on');

plot(ax1, t_plot, v_E_field(idx), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 4);
plot(ax1, t_plot, v_E_centroid(idx), 'Color', 'r', 'LineStyle', '-',  'LineWidth', 4);
plot(ax1, t_plot, v_L_field(idx, 1), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 4);
plot(ax1, t_plot, v_L_centroid(idx, 1), 'Color', 'b', 'LineStyle', '-',  'LineWidth', 4);

xlim(ax1,[t_plot(1) t_plot(end)]); ylim(ax1,[0.0 0.3]);
xlabel(ax1,'Time \tau', 'FontWeight', 'bold');
ylabel(ax1,'Rise Velocity \itV', 'FontWeight', 'bold');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 38, 'FontWeight', 'bold', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 5.5);
legend(ax1, 'DNS–Vol', 'DNS–Cent', ...
       'Lagrang–Vol', 'Lagrang–Cent', ...
       'Orientation','vertical','NumColumns',2,...
       'FontSize',34,'FontName','Times New Roman');

% ---------------- Subplot 2: Circularity profile -----------------------
ax2 = subplot(1,2,2);
hold(ax2,'on'); box(ax2,'on');
plot(ax2, t_plot, circ_E(idx), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 4);
plot(ax2, t_plot, circ_L(idx, 1), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 4);

xlim(ax2,[t_plot(1) t_plot(end)]);
ylim(ax2,[0.6 1.05]);
xlabel(ax2,'Time \tau', 'FontWeight', 'bold');
ylabel(ax2,'Circularity', 'FontWeight', 'bold');
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 38, 'FontWeight', 'bold', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 5.5);
legend(ax2, 'DNS', 'Lagrang', ...
       'Orientation','vertical','NumColumns',1,...
       'FontSize',34,'FontName','Times New Roman');


save('cache_Fig1a.mat', 't_plot','idx','v_E_field','v_E_centroid','v_L_field', ...
    'v_L_centroid','circ_E','circ_L','lagCaseSuffixes','startNum','endNum','nLag','time')