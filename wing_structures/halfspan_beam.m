function B = halfspan_beam(st, ld, opt)
% HALFSPAN_BEAM  Bending and torsion of one half of a wing that is
% symmetric about the centre plane and supported at one mount per side.
%
% Model: beam along the projected span s (0 = centre plane, s_tip = tip).
%   - Centre plane (symmetric loading): slope = 0, torque = 0, no vertical
%     support.
%   - Mount at s_m: vertical/drag reaction and pitch (torque) reaction.
%     opt.bc = 'pinned' (no bending moment reaction) or 'clamped' (also
%     slope = 0 at the mount; inner-span moment found from that condition).
%   The overhang (s > s_m) is statically determinate either way, so the
%   bending moment and torque there do not depend on stiffness.
%   Out-of-plane curvature of the spar (the rise in front view) lies in the
%   plane of the vertical load and needs no axial force on a free
%   overhang, so a straight projected beam gives the same moments.
%
% st  : stations, fields s, EI, GJ, x_ea, y_ea (column vectors)
% ld  : point loads, fields r (n x 3, CFD frame), F (n x 3), C (n x 3 couple
%       about r), s (n x 1); applied as given (scale before calling)
% opt : s_m, s_tip, z_of_s (handle), side (+1 right / -1 left), bc,
%       r_mount_xy ([x y] of the mount on the spar), ns (grid points)
%
% Returns grid s and, per point: M_hog (bending moment, + = hogging, tip
% loads down -> top fibres in tension), T (torque about the elastic axis,
% + = nose up), w (deflection, + = down), phi (twist, + = nose up, rad),
% kappa = M_hog/EI, and the mount reactions.

if ~isfield(opt, 'ns'), opt.ns = 1200; end
if ~isfield(opt, 'bc'), opt.bc = 'pinned'; end
s = unique([linspace(0, opt.s_tip, opt.ns)'; opt.s_m]);
EI  = interp1(st.s, st.EI,  s, 'linear', 'extrap');
GJ  = interp1(st.s, st.GJ,  s, 'linear', 'extrap');
xea = interp1(st.s, st.x_ea, s, 'linear', 'extrap');
yea = interp1(st.s, st.y_ea, s, 'linear', 'extrap');
% hold end values instead of extrapolating far
EI(s < st.s(1))  = st.EI(1);  EI(s > st.s(end))  = st.EI(end);
GJ(s < st.s(1))  = st.GJ(1);  GJ(s > st.s(end))  = st.GJ(end);
xea(s < st.s(1)) = st.x_ea(1); xea(s > st.s(end)) = st.x_ea(end);
yea(s < st.s(1)) = st.y_ea(1); yea(s > st.s(end)) = st.y_ea(end);

n = numel(s);
r_m = [opt.r_mount_xy(1) opt.r_mount_xy(2) opt.z_of_s(opt.s_m)];
R   = -sum(ld.F, 1);                         % mount reaction force
Mint = zeros(n, 3);
for i = 1:n
    rc = [xea(i) yea(i) opt.z_of_s(s(i))];
    o = ld.s > s(i);
    Mi = sum(cross(ld.r(o,:) - rc, ld.F(o,:), 2), 1) + sum(ld.C(o,:), 1);
    if s(i) < opt.s_m
        Mi = Mi + cross(r_m - rc, R);
    end
    Mint(i,:) = Mi;
end
% pitch reaction at the mount: makes the torque at the centre plane zero
CRz = -Mint(1,3);
inner = s < opt.s_m;
Mint(inner,3) = Mint(inner,3) + CRz;

M_hog = opt.side * Mint(:,1);
T     = Mint(:,3);
Mc = 0;
if strcmpi(opt.bc, 'clamped')
    % slope zero at centre and at mount -> integral of M/EI over inner span = 0
    si = s(inner | s == opt.s_m);
    Mc = -trapz(si, M_hog(inner | s == opt.s_m) ./ EI(inner | s == opt.s_m)) / ...
          trapz(si, 1 ./ EI(inner | s == opt.s_m));
    M_hog(inner) = M_hog(inner) + Mc;
end

kappa = M_hog ./ EI;
slope = cumtrapz(s, kappa);                  % w'(0) = 0 (symmetry)
im = find(s == opt.s_m, 1);
if strcmpi(opt.bc, 'clamped')
    % overhang starts with zero slope at the mount
    slope(s >= opt.s_m) = slope(s >= opt.s_m) - slope(im);
end
w = cumtrapz(s, slope);
w = w - w(im);                               % w(s_m) = 0

phi = cumtrapz(s, T ./ GJ);
phi = phi - phi(im);                         % phi(s_m) = 0

B = struct('s', s, 'M_hog', M_hog, 'T', T, 'w', w, 'phi', phi, 'kappa', kappa, ...
    'EI', EI, 'GJ', GJ, 'x_ea', xea, 'y_ea', yea, 'R', R, 'CRz', CRz, 'Mc', Mc, ...
    'im', im, 'bc', opt.bc);
end
