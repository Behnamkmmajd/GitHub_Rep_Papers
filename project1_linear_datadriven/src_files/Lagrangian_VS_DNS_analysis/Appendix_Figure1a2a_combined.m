clc;
clear;
close all;

%% Load cached outputs
S_scattered   = load('cache_Fig1a.mat'); % Lagrangian: scattered
S_interp      = load('cache_Fig2a.mat'); % Lagrangian: interpolated

%% Case selection (use first case: Vlinear_Glinear)
caseIdx = 1;

%% ------------------------------------------------------------------------
%% Figure 1: Rise velocity profiles (6 curves)
fig1 = figure();
set(fig1,'position',[0 0 1300 1300])

ax1 = axes(fig1);
hold(ax1, 'on');
box(ax1, 'on');

% Eulerian (red)
plot(ax1, S_scattered.t_plot, S_scattered.v_E_centroid(S_scattered.idx), ...
    'r-', 'LineWidth', 2.5);
plot(ax1, S_scattered.t_plot, S_scattered.v_E_field(S_scattered.idx), ...
    'r--', 'LineWidth', 2.5);

% Lagrangian interpolated (blue)
plot(ax1, S_interp.t_plot, S_interp.v_L_centroid(S_interp.idx, caseIdx), ...
    'b-', 'LineWidth', 2.5);
plot(ax1, S_interp.t_plot, S_interp.v_L_field(S_interp.idx, caseIdx), ...
    'b--', 'LineWidth', 2.5);

% Lagrangian scattered (green)
plot(ax1, S_scattered.t_plot, S_scattered.v_L_centroid(S_scattered.idx, caseIdx), ...
    'g-', 'LineWidth', 2.5);
plot(ax1, S_scattered.t_plot, S_scattered.v_L_field(S_scattered.idx, caseIdx), ...
    'g--', 'LineWidth', 2.5);

xlim(ax1, [min([S_scattered.t_plot(:); S_interp.t_plot(:)]), max([S_scattered.t_plot(:); S_interp.t_plot(:)])]);
ylim(ax1, [0 0.3]);
xlabel(ax1, 'Time \tau');
ylabel(ax1, 'Non-dimensional Rise Velocity V');

legend(ax1, {
    'Eulerian (Centroid)', 'Eulerian (Averaged)', ...
    'Lagrangian Interp (Centroid)', 'Lagrangian Interp (Averaged)', ...
    'Lagrangian Scatter (Centroid)', 'Lagrangian Scatter (Averaged)'
    }, ...
    'Location', 'northeast', 'NumColumns', 1);

set(ax1,'FontName','Times New Roman','FontSize',34,'FontWeight','bold','XMinorTick','on',...
        'YMinorTick','on','LineWidth',4.5);

%% ------------------------------------------------------------------------
%% Figure 2: Circularity profiles (scattered)
fig2 = figure('Position', [0 0 1200 1200]);
ax2 = axes(fig2);
hold(ax2, 'on');
box(ax2, 'on');

plot(ax2, S_scattered.t_plot, S_scattered.circ_E(S_scattered.idx), ...
    'r-', 'LineWidth', 2.5);
plot(ax2, S_scattered.t_plot, S_scattered.circ_L(S_scattered.idx, caseIdx), ...
    'b-', 'LineWidth', 2.5);

xlim(ax2, [min(S_scattered.t_plot(:)), max(S_scattered.t_plot(:))]);
ylim(ax2, [0.6 1.05]);
xlabel(ax2, 'Time \tau');
ylabel(ax2, 'Circularity');
title(ax2, 'Circularity (Scattered)');
legend(ax2, {'Eulerian', 'Lagrangian'}, 'Location', 'best');
set(ax2,'FontName','Times New Roman','FontSize',34,'FontWeight','bold','XMinorTick','on',...
        'YMinorTick','on','LineWidth',4.5);

%% Figure 3: Circularity profiles (interpolated)
fig3 = figure('Position', [0 0 1200 1200]);
ax3 = axes(fig3);
hold(ax3, 'on');
box(ax3, 'on');

plot(ax3, S_interp.t_plot, S_interp.circ_E(S_interp.idx), ...
    'r-', 'LineWidth', 2.5);
plot(ax3, S_interp.t_plot, S_interp.circ_L(S_interp.idx, caseIdx), ...
    'b-', 'LineWidth', 2.5);

xlim(ax3, [min(S_interp.t_plot(:)), max(S_interp.t_plot(:))]);
ylim(ax3, [0.6 1.05]);
xlabel(ax3, 'Time \tau');
ylabel(ax3, 'Circularity');
title(ax3, 'Circularity (Interpolated)');
legend(ax3, {'Eulerian', 'Lagrangian'}, 'Location', 'best');
set(ax3,'FontName','Times New Roman','FontSize',34,'FontWeight','bold','XMinorTick','on',...
        'YMinorTick','on','LineWidth',4.5);

