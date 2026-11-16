function [T, Tskin, info] = load_main_element(cfg)
% LOAD_MAIN_ELEMENT  Read the main element STL into CFD coordinates and
% separate the skin from the spanwise end caps.
%   T     : all triangles (nT x 9) in CFD frame [m]
%   Tskin : triangles excluding flat end caps (normal along span axis at
%           the spanwise extremes). The caps coincide with the endplate
%           inner faces and would otherwise claim endplate CFD faces.
%   info  : spanwise extent s_min, s_max of the STL [m from centreline]
T = read_stl(cfg.stl, cfg.stl_scale);
T = transform_tris(T, cfg.stl_R, cfg.stl_t);

nrm = cross(T(:,4:6) - T(:,1:3), T(:,7:9) - T(:,1:3), 2);
nrm = nrm ./ sqrt(sum(nrm.^2, 2));
zc  = (T(:,3) + T(:,6) + T(:,9)) / 3;
s_c = cfg.side * (zc - cfg.z_sym);
sv  = cfg.side * ([T(:,3); T(:,6); T(:,9)] - cfg.z_sym);
info.s_min = min(sv); info.s_max = max(sv);
cap = abs(nrm(:,3)) > 0.99 & (abs(s_c - info.s_min) < 1e-3 | abs(s_c - info.s_max) < 1e-3);
Tskin = T(~cap, :);
info.n_cap = nnz(cap);
end
