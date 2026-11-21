function P = station_props(outlines, s_st, lam, pl, x_plate, gap)
% STATION_PROPS  Skin + plate-spar properties at each spanwise station.
%   outlines : cell array of section loops [x forward, y up] (outer mould line)
%   s_st     : station positions [m]
%   lam      : skin laminate (laminate_clt)
%   pl       : plate struct: t, E, G, rho
%   x_plate  : x of the (flat, vertical) spar plate [m]
%   gap      : clearance between plate and skin inner surface, each side [m]
%
% Skin and plate are joined only through ribs, so each bends about its own
% neutral axis (EI = EI_skin + EI_plate) and a vertical load splits by EI.
% The plate has negligible torsional stiffness, so the combined elastic
% axis is the EI-weighted mean of the skin shear centre and the plate:
%   x_ea = (EI_s x_sc + EI_p x_p) / (EI_s + EI_p)
ns = numel(s_st);
z = zeros(ns, 1);
P = struct('s', s_st(:), 'EI_s', z, 'GJ_s', z, 'x_sc', z, 'y_c', z, 'y_top', z, ...
    'y_bot', z, 'A_enc', z, 'm_s', z, 'h_p', z, 'y_p', z, 'EI_p', z, 'GJ_p', z, ...
    'm_p', z, 'EI', z, 'GJ', z, 'x_ea', z, 'y_ea', z, 'x_c', z, 'EIxx_s', z, 'EIxz_s', z, 'EIzz_s', z);
for k = 1:ns
    xy = outlines{k};
    Sk = thinwall_closed_section(xy, lam.h, lam.Ex, lam.Gxy, lam.rho);
    P.EI_s(k) = Sk.EIxx_eff; P.GJ_s(k) = Sk.GJ; P.y_c(k) = Sk.zc;
    % full skin bending stiffnesses about its centroid, for unsymmetric-bending strains (skin_strain)
    P.x_c(k) = Sk.xc; P.EIxx_s(k) = Sk.EIxx; P.EIxz_s(k) = Sk.EIxz; P.EIzz_s(k) = Sk.EIzz;
    P.A_enc(k) = Sk.Aenc; P.m_s(k) = Sk.mass_per_len;
    P.y_top(k) = max(xy(:,2)); P.y_bot(k) = min(xy(:,2));
    P.x_sc(k) = shear_center_closed(xy, lam.h);

    D = section_depth(xy, 2001);
    yh = interp1(D.x, D.y_hi, x_plate); yl = interp1(D.x, D.y_lo, x_plate);
    h = (yh - yl) - 2*lam.h - 2*gap;
    P.h_p(k) = h; P.y_p(k) = (yh + yl)/2;
    P.EI_p(k) = pl.E * pl.t * h^3 / 12;
    P.GJ_p(k) = pl.G * h * pl.t^3 / 3 * (1 - 0.63*pl.t/h);   % thin rectangle
    P.m_p(k)  = pl.rho * pl.t * h;
end
P.EI = P.EI_s + P.EI_p;
P.GJ = P.GJ_s + P.GJ_p;
P.x_ea = (P.EI_s .* P.x_sc + P.EI_p * x_plate) ./ P.EI;
P.y_ea = (P.EI_s .* P.y_c + P.EI_p .* P.y_p) ./ P.EI;
P.x_plate = x_plate;
end
