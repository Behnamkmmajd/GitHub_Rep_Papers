function circ = circularity_alpha_scattered(alpha, x, y, level)
%CIRCULARITY_ALPHA_SCATTERED  Circularity from scattered alpha (axisymmetric).
%
%   circ = circularity_alpha_scattered(alpha, x, y, level)
%
%   Mirror points across r = 0, threshold at level, use boundary()
%   (concave hull) to capture indentations, circ = 4*pi*A/P^2.

    if nargin < 4, level = 0.5; end

    circ = NaN;
    if isempty(alpha), return; end

    mask = alpha(:) >= level;
    if nnz(mask) < 3, return; end

    xs = x(mask);  ys = y(mask);

    % Mirror about symmetry axis
    xf = [xs(:); -xs(:)];
    yf = [ys(:);  ys(:)];

    % Concave boundary (shrink factor 0.8: captures concavities)
    try,  idx = boundary(xf, yf, 0.8);  catch, return;  end

    xb = xf(idx);  yb = yf(idx);

    A = polyarea(xb, yb);
    P = sum(hypot(diff(xb), diff(yb)));
    if P == 0, return; end

    circ = (4 * pi * A) / P^2;
end
