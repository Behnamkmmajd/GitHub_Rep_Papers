function [centroid_y_k, vy_centroid_k] = rise_velocity_centroid_alpha_3D(alpha, YG, centroid_y_prev, t_k, t_km1)
%RISE_VELOCITY_CENTROID_ALPHA_3D  Centroid rise velocity for full 3-D fields.
%
%   Unlike the 2-D axisymmetric version (which uses radial weighting),
%   each 3-D cell already represents its true volume, so the centroid
%   is simply  y_c = sum(alpha * y) / sum(alpha).
%
%   Inputs:
%     alpha           - 3-D gas-fraction field
%     YG              - 3-D mesh of y (vertical) coordinates
%     centroid_y_prev - centroid position at previous step (NaN for first step)
%     t_k             - time at current step
%     t_km1           - time at previous step
%
%   Outputs:
%     centroid_y_k    - centroid position at time k
%     vy_centroid_k   - rise velocity dy/dt (NaN if not computable)

    centroid_y_k = sum(alpha(:) .* YG(:)) / sum(alpha(:));

    vy_centroid_k = NaN;
    dt_c = t_k - t_km1;
    if ~isnan(centroid_y_k) && ~isnan(centroid_y_prev) && dt_c ~= 0
        vy_centroid_k = (centroid_y_k - centroid_y_prev) / dt_c;
    end
end
