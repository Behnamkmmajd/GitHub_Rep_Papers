function [x_nd, y_nd, z_nd, t_nd, box_centre] = grid_generator_3D(k, P)
%GRID_GENERATOR_3D  Reconstruct the 3-D semi-Lagrangian box grid at frame k.
%
%  [x_nd, y_nd, z_nd, t_nd]             = grid_generator_3D(k, P)
%  [x_nd, y_nd, z_nd, t_nd, box_centre] = grid_generator_3D(k, P)
%
%  Produces the EXACT same grid that extract_3D_lagrangian.py used for
%  the Lagrangian interpolation, given the pre-computed data stored in
%  case_params.mat.
%
%  The box tracks the bubble centroid in all three directions (x, y, z).
%  At each frame the box centre is obtained by integrating the bubble
%  velocity (V_x, V_rise, V_z) from the initial centroid position.
%
%  INPUTS
%    k  – frame index  (1-based,  1 ≤ k ≤ P.nframes)
%    P  – struct from  load('<case>/case_params.mat')
%
%  OUTPUTS
%    x_nd       – (Nx, 1)  absolute box x-coordinates / D at frame k
%    y_nd       – (Ny, 1)  absolute box y-coordinates / D at frame k
%    z_nd       – (Nz, 1)  absolute box z-coordinates / D at frame k
%    t_nd       – scalar    non-dimensional time at frame k
%    box_centre – [1×3]     box centre [x, y, z] / D at frame k
%
%  All three coordinate vectors shift with the bubble every frame.
%
%  EXAMPLE
%    P = load('3D_case_Re168_Bo12/case_params.mat');
%    [x, y, z, t, bc] = grid_generator_3D(100, P);
%    fprintf('Box centre: (%.4f, %.4f, %.4f) D\n', bc);

    % ── Relative coordinates (same for all frames) ──
    x_rel = P.x_rel_nd(:);
    y_rel = P.y_rel_nd(:);
    z_rel = P.z_rel_nd(:);

    % ── Integrate velocities → box centre at frame k ──
    Vx = P.V_x_nd(:);
    Vy = P.V_rise_nd(:);
    Vz = P.V_z_nd(:);
    dt = P.dt_nd;

    box_x = P.box_x0_nd;
    box_y = P.box_y0_nd;
    box_z = P.box_z0_nd;
    for i = 2:k
        box_x = box_x + Vx(i) * dt;
        box_y = box_y + Vy(i) * dt;
        box_z = box_z + Vz(i) * dt;
    end

    % ── Absolute coordinates at frame k ──
    x_nd = x_rel + box_x;
    y_nd = y_rel + box_y;
    z_nd = z_rel + box_z;

    % ── Non-dimensional time ──
    t_nd = (k - 1) * dt;

    % ── Box centre ──
    box_centre = [box_x, box_y, box_z];
end
