function L = integrate_strips(xyz, Ff, Amag, s, edges, ref_xy, z_of_s)
% INTEGRATE_STRIPS  Sum face forces into spanwise strips.
%   xyz    : N x 3 face centroids (CFD frame)
%   Ff     : N x 3 face forces [N]
%   Amag   : N x 1 face areas [m^2]
%   s      : N x 1 spanwise coordinate of each face [m from centreline]
%   edges  : strip edges in s (uniform spacing; the last strip may be
%            shorter, faces just beyond the last edge go into it)
%   ref_xy : [x y] of the moment reference line (constant along span);
%            moments for strip k are about [ref_xy(1) ref_xy(2) z_of_s(s_mid(k))]
%   z_of_s : function handle mapping s to CFD z
%
% Returns per strip: s_mid, ds, F (nb x 3), M (nb x 3), area, and per-unit-
% span versions f = F/ds, m = M/ds. Move moments to another reference
% point r with  M_r = M + (r0 - r) x F.

edges = edges(:);
nb = numel(edges) - 1;
ds = diff(edges);
bin = floor((s - edges(1)) / ds(1)) + 1;
bin(bin == nb + 1 & s < edges(1) + nb*ds(1)) = nb;   % shortened last strip
ok  = bin >= 1 & bin <= nb;

L.s_mid = (edges(1:end-1) + edges(2:end)) / 2;
L.ds = ds;
L.F = zeros(nb, 3); L.M = zeros(nb, 3); L.area = zeros(nb, 1); L.nfaces = zeros(nb, 1);
for k = 1:nb
    m = ok & bin == k;
    if ~any(m), continue; end
    r0 = [ref_xy(1) ref_xy(2) z_of_s(L.s_mid(k))];
    r  = xyz(m,:) - r0;
    L.F(k,:) = sum(Ff(m,:), 1);
    L.M(k,:) = sum(cross(r, Ff(m,:), 2), 1);
    L.area(k) = sum(Amag(m));
    L.nfaces(k) = nnz(m);
end
L.f = L.F ./ ds;
L.m = L.M ./ ds;
L.ref_xy = ref_xy;
L.outside = nnz(~ok);
end
