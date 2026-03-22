function sph = sphericity_alpha_grid_3D(alpha, x, y, z, level)
%SPHERICITY_ALPHA_GRID_3D  Sphericity from a 3-D gridded alpha field.
%
%   sph = sphericity_alpha_grid_3D(alpha, x, y, z, level)
%
%   Extracts the isosurface at the given level, computes surface area
%   and enclosed volume from the triangulation, then returns:
%
%       sph = pi^(1/3) * (6*V)^(2/3) / A_s
%
%   sph = 1 for a perfect sphere;  sph < 1 for deformed shapes.

    if nargin < 5, level = 0.5; end

    sph = NaN;
    if isempty(alpha), return; end

    [X, Y, Z] = meshgrid(x, y, z);

    % Extract isosurface triangulation
    fv = isosurface(X, Y, Z, alpha, level);
    if isempty(fv.vertices) || size(fv.faces, 1) < 4, return; end

    verts = fv.vertices;   % [N x 3]
    faces = fv.faces;      % [M x 3]

    % Triangle vertices
    a = verts(faces(:,1), :);
    b = verts(faces(:,2), :);
    c = verts(faces(:,3), :);

    % Surface area  =  sum of triangle areas
    cr = cross(b - a, c - a, 2);
    A_s = sum(sqrt(sum(cr.^2, 2))) / 2;

    % Enclosed volume via divergence theorem (signed tetrahedra with origin)
    V = abs(sum(dot(a, cr, 2))) / 6;

    if A_s == 0 || V == 0, return; end

    sph = pi^(1/3) * (6 * V)^(2/3) / A_s;
end
