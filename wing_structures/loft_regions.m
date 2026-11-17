function L = loft_regions(s, y_top, y_bot, slope_min)
% LOFT_REGIONS  Spanwise regions where the wing is lofted, i.e. where both
% skin lines slope in front view (the whole section rises or drops).
%   s            : spanwise positions [m], uniform spacing
%   y_top, y_bot : height of the top and bottom skin lines in front view [m]
%   slope_min    : |slope| of both lines above which the wing counts as lofted [rad]
%
% Returns L with one row [start end] per region [m]. Each region is widened
% by one grid step at both ends, so start and end lie on constant section.
% Regions that touch after widening are merged.
s = s(:); ds = s(2) - s(1);
% lofted = BOTH skin lines slope: the whole section moves up or down. A slope in
% one line only (the section changing thickness or shape) does not count.
sl = min(abs(atan(gradient(y_top(:), ds))), abs(atan(gradient(y_bot(:), ds))));
on = double(sl > slope_min);
d = diff([0; on; 0]);
i0 = find(d == 1); i1 = find(d == -1) - 1;
L = zeros(numel(i0), 2);
for k = 1:numel(i0)
    L(k,:) = [s(max(i0(k) - 1, 1)), s(min(i1(k) + 1, numel(s)))];
end
k = 1;
while k < size(L, 1)
    if L(k+1,1) <= L(k,2)
        L(k,2) = L(k+1,2); L(k+1,:) = [];
    else
        k = k + 1;
    end
end
end
