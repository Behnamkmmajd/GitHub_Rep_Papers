function vz_mean_k = rise_velocity_mean_alpha(alpha, XG, vz)
%RISE_VELOCITY_MEAN_ALPHA  Volume-averaged rise velocity using alpha-weighting.
%   (2-D axisymmetric, no threshold mask — alpha weights all cells)
%
%   Inputs:
%     alpha - 2-D gas-fraction field
%     XG    - 2-D mesh of r (radial) coordinates
%     vz    - 2-D vertical velocity field
%
%   Output:
%     vz_mean_k - mean rise velocity: sum(alpha*r*vz) / sum(alpha*r)

    vz_mean_k = sum(alpha(:) .* XG(:) .* vz(:)) / sum(alpha(:) .* XG(:));
end
