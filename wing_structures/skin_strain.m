function [e_hi, e_lo] = skin_strain(xy, xc, yc, EIxx, EIxz, EIzz, M1, M2)
% SKIN_STRAIN  Largest and smallest axial strain in the skin of one section
% under bending about both section axes (unsymmetric bending).
%   xy               : section outline [x forward, y up], M x 2
%   xc, yc           : modulus-weighted centroid of the skin (scalars)
%   EIxx, EIxz, EIzz : bending stiffnesses about axes through the centroid
%                      (EIxx resists vertical loads, EIzz chordwise loads,
%                      EIxz is the product term). Add a spar's EI to EIxx
%                      before calling.
%   M1 : bending moment from vertical loads, + = hogging (halfspan_beam M_hog)
%   M2 : bending moment from chordwise loads (halfspan_beam M_lag); 0 to ignore
% The stiffnesses and moments may be vectors, one entry per load case.
%
%   strain(x, y) = a (y - yc) + b (x - xc),
%   [EIxx EIxz; EIxz EIzz] [a; b] = [M1; M2]
%
% With M2 = 0 this is a vertical load on a section that is free to deflect
% sideways: a = M1 / (EIxx - EIxz^2/EIzz), b = -a EIxz / EIzz, so the
% neutral axis is tilted by EIxz/EIzz instead of being horizontal.
%   e_hi, e_lo : max and min strain over the outline, one row per load case
EIxx = EIxx(:); EIxz = EIxz(:); EIzz = EIzz(:); M1 = M1(:); M2 = M2(:);
D = EIxx .* EIzz - EIxz.^2;
a = (M1 .* EIzz - M2 .* EIxz) ./ D;
b = (M2 .* EIxx - M1 .* EIxz) ./ D;
e = (xy(:,2) - yc) * a' + (xy(:,1) - xc) * b';       % points x load cases
e_hi = max(e, [], 1)'; e_lo = min(e, [], 1)';
end
