function S = thinwall_closed_section(xz, t, E, G, rho)
% THINWALL_CLOSED_SECTION  Modulus-weighted properties of a single-cell
% closed thin-walled section, built from straight wall segments.
%
%   xz  : N x 2 wall MIDLINE points [x z] in metres, in order around the
%         cell (closing segment from last to first point is added).
%         x = chordwise, z = thickness direction (downforce acts along z).
%   t   : wall thickness, scalar or N x 1 (one per segment i -> i+1)
%   E   : axial (spanwise) modulus per segment, scalar or N x 1
%   G   : in-plane shear modulus per segment, scalar or N x 1
%   rho : density per segment, scalar or N x 1
%
%   Bending stiffnesses are about axes through the modulus-weighted
%   centroid:
%     EIxx  resists bending from z loads (downforce)
%     EIzz  resists bending from x loads (drag)
%     EIxz  product term (non-zero for cambered sections)
%     EIxx_eff = EIxx - EIxz^2/EIzz : effective stiffness for a pure z load
%                when the section is free to deflect in x as well.
%   GJ from Bredt-Batho: GJ = 4 A^2 / sum(L/(G t)).
%   Straight segments are integrated exactly (including each segment's own
%   second moment), so a polygon contour is exact; a curved airfoil is
%   approximated by its point density.

N = size(xz, 1);
t   = expand(t, N);  E = expand(E, N);
G   = expand(G, N);  rho = expand(rho, N);

p1 = xz;
p2 = xz([2:N 1], :);
d  = p2 - p1;
dx = d(:,1); dz = d(:,2);
L  = sqrt(dx.^2 + dz.^2);
xm = (p1(:,1) + p2(:,1)) / 2;
zm = (p1(:,2) + p2(:,2)) / 2;

w   = E .* t .* L;              % modulus-weighted segment "area"
EA  = sum(w);
xc  = sum(w .* xm) / EA;
zc  = sum(w .* zm) / EA;

S.EA    = EA;
S.xc    = xc;
S.zc    = zc;
S.EIxx  = sum(w .* ((zm - zc).^2 + dz.^2/12));
S.EIzz  = sum(w .* ((xm - xc).^2 + dx.^2/12));
S.EIxz  = sum(w .* ((xm - xc).*(zm - zc) + dx.*dz/12));
S.EIxx_eff = S.EIxx - S.EIxz^2 / S.EIzz;

% Enclosed area of the midline (shoelace)
Aenc  = 0.5 * abs(sum(p1(:,1).*p2(:,2) - p2(:,1).*p1(:,2)));
S.Aenc = Aenc;
S.GJ   = 4 * Aenc^2 / sum(L ./ (G .* t));

S.mass_per_len = sum(rho .* t .* L);   % kg/m
S.perimeter    = sum(L);
end

function v = expand(v, N)
if isscalar(v), v = v * ones(N, 1); else, v = v(:); end
end
