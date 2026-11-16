% RUN_TESTS_STAGE2  Shear centre and half-span beam against exact results.
% Run from the wing_structures folder:  addpath('tests'); run_tests_stage2
nfail = 0;

%% 1. Shear centre of a thin circle = its centre; resultant = unit shear
th = linspace(0, 2*pi, 2001)'; th(end) = [];
xy = [0.3 + 0.05*cos(th), 0.1 + 0.05*sin(th)];
[xs, inf1] = shear_center_closed(xy, 1e-3);
nfail = nfail + check('Shear centre: circle at centre', xs, 0.3, 1e-6);
nfail = nfail + check('Shear centre: |resultant| = 1', abs(inf1.Fy), 1, 1e-4);
nfail = nfail + check('Shear centre: no horizontal resultant', abs(inf1.Fx) + 1, 1, 1e-4);

%% 2. Translation invariance, and symmetric box at its centre
xyb = rect_tube_midline(0.1, 0.04, 1e-3) + [0.5 0.2];
nfail = nfail + check('Shear centre: symmetric box', shear_center_closed(xyb, 1e-3), 0.5, 1e-9);
nfail = nfail + check('Shear centre: translation', shear_center_closed(xyb + [0.1 0], 1e-3), 0.6, 1e-9);

%% 3. Thicker right web moves the shear centre right
% segments of rect_tube_midline: bottom, right, top, left
xs_thick = shear_center_closed(xyb, [1e-3; 3e-3; 1e-3; 1e-3]);
nfail = nfail + check('Shear centre: moves toward thick web', double(xs_thick > 0.5), 1, 0);

%% 4. Beam: uniform load on the overhang only, constant EI (pinned mount)
L = 0.65; sm = 0.2; a = L - sm; EI = 3000; GJ = 2000; q = 100; t0 = 10;
ns = 2000; edges = linspace(sm, L, ns+1)'; sk = (edges(1:end-1) + edges(2:end))/2; ds = diff(edges);
st = struct('s', [0; L], 'EI', [EI; EI], 'GJ', [GJ; GJ], 'x_ea', [0.8; 0.8], 'y_ea', [0.1; 0.1]);
ld.s = sk; ld.r = [0.8*ones(ns,1) 0.1*ones(ns,1) -sk];
ld.F = [zeros(ns,1) -q*ds zeros(ns,1)];
ld.C = [zeros(ns,2) t0*ds];
opt = struct('s_m', sm, 's_tip', L, 'z_of_s', @(s) -s, 'side', -1, 'bc', 'pinned', ...
             'r_mount_xy', [0.8 0.1], 'ns', 4001);
Bp = halfspan_beam(st, ld, opt);
nfail = nfail + check('Beam: mount moment q a^2/2', Bp.M_hog(Bp.im), q*a^2/2, 1e-4);
nfail = nfail + check('Beam: inner span moment constant', Bp.M_hog(1), q*a^2/2, 1e-4);
nfail = nfail + check('Beam (pinned): tip deflection', Bp.w(end), q*a^4/(8*EI) + q*a^3*sm/(2*EI), 1e-4);
nfail = nfail + check('Beam: tip twist t0 a^2/(2GJ)', Bp.phi(end), t0*a^2/(2*GJ), 1e-4);
nfail = nfail + check('Beam: reaction = total load (up)', Bp.R(2), q*a, 1e-9);
opt.bc = 'clamped';
Bc = halfspan_beam(st, ld, opt);
nfail = nfail + check('Beam (clamped): tip deflection', Bc.w(end), q*a^4/(8*EI), 1e-4);

%% 5. Torque over the whole half span: centre twist and zero centre torque
ns2 = 3000; e2 = linspace(0, L, ns2+1)'; s2 = (e2(1:end-1) + e2(2:end))/2; d2 = diff(e2);
ld2.s = s2; ld2.r = [0.8*ones(ns2,1) 0.1*ones(ns2,1) -s2];
ld2.F = zeros(ns2, 3); ld2.C = [zeros(ns2,2) t0*d2];
opt.bc = 'pinned';
B2 = halfspan_beam(st, ld2, opt);
nfail = nfail + check('Beam: centre torque = 0', abs(B2.T(1)) + 1, 1, 1e-6);
nfail = nfail + check('Beam: centre twist t0 sm^2/(2GJ)', B2.phi(1), t0*sm^2/(2*GJ), 5e-3);
nfail = nfail + check('Beam: tip twist (whole span)', B2.phi(end), t0*a^2/(2*GJ), 1e-3);

fprintf('\n%d test(s) failed.\n', nfail);
