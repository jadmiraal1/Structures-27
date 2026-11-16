function [d, Q] = point_triangle_distance(P, a, b, c)
% POINT_TRIANGLE_DISTANCE  Exact distance from points P (N x 3) to
% triangles with vertices a, b, c. Each of a, b, c is either 1 x 3 (one
% triangle for all points) or N x 3 (one triangle per point, pairwise).
% Returns distance d (N x 1) and closest points Q (N x 3).
%
% Voronoi-region method from Ericson, "Real-Time Collision Detection",
% sec. 5.1.5, vectorised over points.

N  = size(P, 1);
if size(a,1) == 1, a = repmat(a, N, 1); b = repmat(b, N, 1); c = repmat(c, N, 1); end
ab = b - a; ac = c - a; bc = c - b;
Q  = zeros(N, 3);
todo = true(N, 1);

ap = P - a;
d1 = sum(ap.*ab, 2); d2 = sum(ap.*ac, 2);
m = todo & d1 <= 0 & d2 <= 0;                 % vertex region a
Q(m,:) = a(m,:); todo(m) = false;

bp = P - b;
d3 = sum(bp.*ab, 2); d4 = sum(bp.*ac, 2);
m = todo & d3 >= 0 & d4 <= d3;                % vertex region b
Q(m,:) = b(m,:); todo(m) = false;

vc = d1.*d4 - d3.*d2;
m = todo & vc <= 0 & d1 >= 0 & d3 <= 0;       % edge ab
v = d1(m) ./ (d1(m) - d3(m));
Q(m,:) = a(m,:) + v(:) .* ab(m,:); todo(m) = false;

cp = P - c;
d5 = sum(cp.*ab, 2); d6 = sum(cp.*ac, 2);
m = todo & d6 >= 0 & d5 <= d6;                % vertex region c
Q(m,:) = c(m,:); todo(m) = false;

vb = d5.*d2 - d1.*d6;
m = todo & vb <= 0 & d2 >= 0 & d6 <= 0;       % edge ac
w = d2(m) ./ (d2(m) - d6(m));
Q(m,:) = a(m,:) + w(:) .* ac(m,:); todo(m) = false;

va = d3.*d6 - d5.*d4;
m = todo & va <= 0 & (d4 - d3) >= 0 & (d5 - d6) >= 0;   % edge bc
w = (d4(m) - d3(m)) ./ ((d4(m) - d3(m)) + (d5(m) - d6(m)));
Q(m,:) = b(m,:) + w(:) .* bc(m,:); todo(m) = false;

m = todo;                                      % face interior
den = va(m) + vb(m) + vc(m);
v = vb(m) ./ den; w = vc(m) ./ den;
Q(m,:) = a(m,:) + v(:) .* ab(m,:) + w(:) .* ac(m,:);

d = sqrt(sum((P - Q).^2, 2));
end
