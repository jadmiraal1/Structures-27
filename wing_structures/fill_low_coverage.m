function [L, filled] = fill_low_coverage(L, coverage, cmin)
% FILL_LOW_COVERAGE  Replace strips whose classification coverage is below
% CMIN with values interpolated linearly (in s) from the well-covered
% strips. Force and moment per unit span (f, m) are interpolated, then
% F = f*ds and M = m*ds are rebuilt. Strips beyond the last good strip take
% the nearest good value.
%   filled : logical, true for strips that were replaced
filled = coverage(:) < cmin;
good = ~filled;
if nnz(good) < 2
    error('fill_low_coverage:good', 'Fewer than two well-covered strips; cannot interpolate.');
end
s = L.s_mid;
for c = 1:3
    L.f(filled, c) = interp_fill(s(good), L.f(good, c), s(filled));
    L.m(filled, c) = interp_fill(s(good), L.m(good, c), s(filled));
end
L.F = L.f .* L.ds;
L.M = L.m .* L.ds;
L.filled = filled;
end

function yq = interp_fill(x, y, xq)
yq = interp1(x, y, xq, 'linear');
nanq = isnan(yq);
if any(nanq)
    yq(nanq) = interp1(x, y, xq(nanq), 'nearest', 'extrap');
end
end
