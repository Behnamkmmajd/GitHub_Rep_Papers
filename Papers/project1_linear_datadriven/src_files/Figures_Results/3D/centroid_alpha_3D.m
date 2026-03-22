function pos = centroid_alpha_3D(alpha, XG, YG, ZG)
%CENTROID_ALPHA_3D  3-D centroid position from alpha-weighted field.
%
%   pos = centroid_alpha_3D(alpha, XG, YG, ZG)
%
%   Returns pos = [x_c, y_c, z_c] where each component is
%   sum(alpha * coord) / sum(alpha).  No radial weighting (full 3-D).

    w = alpha(:);
    S = sum(w);
    pos = [sum(w .* XG(:)) / S, ...
           sum(w .* YG(:)) / S, ...
           sum(w .* ZG(:)) / S];
end
