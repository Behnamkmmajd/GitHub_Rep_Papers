function circ = circularity_alpha_grid(alpha, x, y, level)
%CIRCULARITY_ALPHA_GRID  Circularity from gridded alpha on an axisymmetric half-plane.
%
%   circ = circularity_alpha_grid(alpha, x, y, level)
%
%   Before contouring, the first grid column (closest to the symmetry axis)
%   is replaced with the second column.  This prevents a spurious closing
%   segment along r ≈ 0 from inflating the perimeter.
%
%   A_full = 2 * A_half,  P_full = 2 * P_half,  circ = 4*pi*A / P^2.

    if nargin < 4, level = 0.5; end

    circ = NaN;
    if isempty(alpha) || isempty(x) || isempty(y), return; end

    xv = x(:)';
    yv = y(:)';
    A  = alpha;

    % Replace first column with second column (remove axis artifact)
    A(:, 1) = A(:, 2);

    % Extract iso-contour
    C = contourc(xv, yv, A, [level level]);
    if isempty(C), return; end

    % Pick the longest contour segment
    best = [];  k = 1;
    while k < size(C, 2)
        n = C(2, k);
        if n < 3 || k + n > size(C, 2), break; end
        seg = C(:, k+1:k+n);
        if size(seg, 2) > size(best, 2), best = seg; end
        k = k + n + 1;
    end
    if isempty(best) || size(best, 2) < 3, return; end

    xb = best(1,:);
    yb = best(2,:);

    % Half-domain area and perimeter
    A_half = polyarea(xb, yb);
    P_half = sum(hypot(diff(xb), diff(yb)));

    if P_half == 0, return; end

    circ = (4 * pi * 2 * A_half) / (2 * P_half)^2;
end
