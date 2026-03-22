clc
clear
close all

%% Settings
l_s = 1;        % length scale (set 1 for no scaling)
stride = 10;    % plot every 'stride' time step
ms = 5;         % marker size

folder = "Results/TrajOnly_dataRe40Bo40";
file   = "trajectories_0_405_Grid250x750_dataRe40Bo40_Vlinear_Glinear.mat";

%% Load (reads ONE .mat file)
S = load(fullfile(folder, file));

X = permute(S.X0traj, [3, 2, 1]);
C = permute(S.properties_trajec, [3, 2, 1]);

% Limits
xMax = max(S.XX(:)) / l_s;
yMax = max(S.YY(:)) / l_s;

%% Plot
figure(1)
clf
% colormap(flipud(hot))
colorbar

for t = 1:stride:size(X, 3)
    scatter(X(:,1,t)/l_s, X(:,2,t)/l_s, ms, C(:,1,t), 'filled')
    xlim([0 xMax])
    ylim([0 yMax])
    daspect([1 1 1])
    title(sprintf('t index = %d', t))
    drawnow
end