function [x_nd, y_nd, t_nd, box_y] = grid_generator_2D(k, P)
%GRID_GENERATOR_2D  Reconstruct the semi-Lagrangian box grid at frame k.
%
%  [x_nd, y_nd, t_nd]        = grid_generator_2D(k, P)
%  [x_nd, y_nd, t_nd, box_y] = grid_generator_2D(k, P)
%
%  Produces the EXACT same grid that extract_2D_lagrangian.py used for
%  the Lagrangian interpolation, given the pre-computed data stored in
%  case_params.mat.
%
%  INPUTS
%    k  – frame index  (1-based,  1 ≤ k ≤ P.nframes)
%    P  – struct from  load('<case>/case_params.mat')
%
%  OUTPUTS
%    x_nd  – (Nx, 1)  box x-coordinates / D   (constant across frames)
%    y_nd  – (Ny, 1)  absolute box y-coords / D at frame k
%    t_nd  – scalar   non-dimensional time at frame k
%    box_y – scalar   box centre y / D at frame k
%
%  The full 2-D mesh for frame k (identical to the extraction grid):
%    [X, Y] = meshgrid(x_nd, y_nd);     % size (Ny, Nx)
%
%  EXAMPLE
%    P = load('Bo100Re100Sc1DR10VR10/case_params.mat');
%    [x, y, t] = grid_generator_2D(100, P);
%    [X, Y] = meshgrid(x, y);

    % ── Box x-coordinates (same for all frames) ──
    x_nd = P.x_nd(:);

    % ── Relative y-coordinates (same for all frames) ──
    y_rel = P.y_rel_nd(:);

    % ── Integrate rise velocity → box centre at frame k ──
    V  = P.V_rise_nd(:);
    dt = P.dt_nd;
    box_y = P.box_y0_nd;
    for i = 2:k
        box_y = box_y + V(i) * dt;
    end

    % ── Absolute y at frame k ──
    y_nd = y_rel + box_y;

    % ── Non-dimensional time ──
    t_nd = (k - 1) * dt;
end
