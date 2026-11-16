function [cx, cy, A] = polygon_centroid(x, y)
% POLYGON_CENTROID  Area centroid of a closed polygon (vertices in order,
% last not repeated). A is the unsigned area.
x = x(:); y = y(:);
x2 = x([2:end 1]); y2 = y([2:end 1]);
cr = x.*y2 - x2.*y;
As = sum(cr) / 2;
cx = sum((x + x2).*cr) / (6*As);
cy = sum((y + y2).*cr) / (6*As);
A = abs(As);
end
