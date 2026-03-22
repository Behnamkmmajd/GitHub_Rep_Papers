function vy_mean_k = rise_velocity_mean_alpha_3D(alpha, vy)
%RISE_VELOCITY_MEAN_ALPHA_3D  Volume-averaged rise velocity for full 3-D fields.
%
%   Unlike the 2-D axisymmetric version (which uses radial weighting),
%   each 3-D cell already represents its true volume, so the mean
%   velocity is simply  v_mean = sum(alpha * vy) / sum(alpha).
%
%   Inputs:
%     alpha - 3-D gas-fraction field
%     vy    - 3-D vertical velocity field
%
%   Output:
%     vy_mean_k - mean rise velocity

    vy_mean_k = sum(alpha(:) .* vy(:)) / sum(alpha(:));
end
