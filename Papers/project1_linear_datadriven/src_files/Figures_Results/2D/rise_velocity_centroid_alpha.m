function [centroid_z_k, vz_centroid_k] = rise_velocity_centroid_alpha(alpha, XG, YG, centroid_z_prev, t_k, t_km1)
%RISE_VELOCITY_CENTROID_ALPHA  Centroid rise velocity using alpha-weighting.
%   (2-D axisymmetric, no threshold mask — alpha weights all cells)
%
%   Inputs:
%     alpha           - 2-D gas-fraction field at time k
%     XG              - 2-D mesh of r (radial) coordinates
%     YG              - 2-D mesh of y (vertical) coordinates
%     centroid_z_prev - centroid position at previous step (NaN for first step)
%     t_k             - time at current step
%     t_km1           - time at previous step
%
%   Outputs:
%     centroid_z_k    - centroid position at time k
%     vz_centroid_k   - rise velocity dy/dt (NaN if not computable)

    % Centroid: y_c = sum(alpha * r * y) / sum(alpha * r)
    centroid_z_k = sum(alpha(:) .* XG(:) .* YG(:)) / sum(alpha(:) .* XG(:));

    % Rise velocity: dc/dt
    vz_centroid_k = NaN;
    dt_c = t_k - t_km1;
    if ~isnan(centroid_z_k) && ~isnan(centroid_z_prev) && dt_c ~= 0
        vz_centroid_k = (centroid_z_k - centroid_z_prev) / dt_c;
    end
end
