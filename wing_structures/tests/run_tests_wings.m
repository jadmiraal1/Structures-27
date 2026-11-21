% RUN_TESTS_WINGS  Tests for the pieces added to run more than one wing:
% wing_config, skin_strain / beam_skin_strains, halfspan_beam M_lag.
% Run from wing_structures:  addpath('tests'); run_tests_wings
addpath(fileparts(mfilename('fullpath')));
fprintf('--- wing_config ---\n');
wing_config('front_wing'); c = wing_config();
check('default/front config name', strcmp(c.name, 'front_wing'), true, 0);
if exist('rear_wing_config', 'file')
    wing_config('rear_wing'); c = wing_config();
    check('switch to rear wing', strcmp(c.name, 'rear_wing'), true, 0);
    wing_config('front_wing');
end

fprintf('--- skin_strain: thin rectangular box, tilted ---\n');
% box b x h, wall t, modulus E, rotated by th about its centroid; a vertical
% moment M bends it about both principal axes: strain = M cos(th) v / EIu
% - M sin(th) u / EIv in principal coordinates (u along b, v along h).
b = 0.30; h = 0.04; t = 1e-3; E = 50e9; th = deg2rad(10); M = 100;
n = 400; u = [linspace(-b/2, b/2, n)'; b/2*ones(n,1); linspace(b/2, -b/2, n)'; -b/2*ones(n,1)];
v = [-h/2*ones(n,1); linspace(-h/2, h/2, n)'; h/2*ones(n,1); linspace(h/2, -h/2, n)'];
xy = [u*cos(th) - v*sin(th), u*sin(th) + v*cos(th)];
S = thinwall_closed_section(xy, t, E, 20e9, 1500);
EIu = E*t*(2*b*(h/2)^2 + 2*h^3/12); EIv = E*t*(2*h*(b/2)^2 + 2*b^3/12);
e_exact = M*cos(th)*v/EIu - M*sin(th)*u/EIv;
[e_hi, e_lo] = skin_strain(xy, S.xc, S.zc, S.EIxx, S.EIxz, S.EIzz, M, 0);
% (a vertical moment with free sideways deflection; exact solution resolves M on the principal axes)
check('max strain, vertical moment, tilted box', e_hi, max(e_exact), 2e-3);
check('min strain, vertical moment, tilted box', e_lo, min(e_exact), 2e-3);
check('EIxx_eff matches principal-axis result', S.EIxx_eff, 1/(cos(th)^2/EIu + sin(th)^2/EIv), 2e-3);
% chordwise moment only, untilted box: strain = M2 * x / EIv
xy0 = [u v]; S0 = thinwall_closed_section(xy0, t, E, 20e9, 1500);
[e_hi, e_lo] = skin_strain(xy0, S0.xc, S0.zc, S0.EIxx, S0.EIxz, S0.EIzz, 0, M);
check('chordwise moment only: max strain', e_hi, M*(b/2)/EIv, 2e-3);
check('chordwise moment only: min strain', e_lo, -M*(b/2)/EIv, 2e-3);
% both moments, tilted: superpose in principal axes. M1 = int(sigma*y), M2 = int(sigma*x)
M2 = 40; Mu = M*cos(th) - M2*sin(th); Mv = M*sin(th) + M2*cos(th);   % about u-axis (int sigma v), v-axis (int sigma u)
e_exact = Mu*v/EIu + Mv*u/EIv;
[e_hi, e_lo] = skin_strain(xy, S.xc, S.zc, S.EIxx, S.EIxz, S.EIzz, M, M2);
check('both moments, tilted box: max strain', e_hi, max(e_exact), 2e-3);
check('both moments, tilted box: min strain', e_lo, min(e_exact), 2e-3);

fprintf('--- halfspan_beam M_lag and beam_skin_strains ---\n');
% straight untilted box wing, left half (side = -1), tip force [-D -L 0]
ns = 6; st.s = linspace(0, 0.5, ns)'; z1 = ones(ns,1);
P = struct('s', st.s, 'EI', S0.EIxx*z1, 'GJ', S0.GJ*z1, 'x_ea', 0*z1, 'y_ea', 0*z1, 'y_c', 0*z1, 'x_c', 0*z1, ...
    'y_top', h/2*z1, 'y_bot', -h/2*z1, 'EI_p', 0*z1, 'EIxx_s', S0.EIxx*z1, 'EIxz_s', 0*z1, 'EIzz_s', S0.EIzz*z1);
Lf = 200; Df = 80; s_tip = 0.5; s_m = 0.2;
ld = struct('s', s_tip, 'r', [0 0 -s_tip], 'F', [-Df -Lf 0], 'C', [0 0 0]);
opt = struct('s_m', s_m, 's_tip', s_tip, 'z_of_s', @(s) -s, 'side', -1, 'bc', 'pinned', 'r_mount_xy', [0 0], 'ns', 500);
B = halfspan_beam(P, ld, opt);
i = find(B.s >= 0.35, 1);
check('M_hog on the overhang', B.M_hog(i), Lf*(s_tip - B.s(i)), 1e-9);
check('M_lag on the overhang (tip pushed aft, + = LE in tension)', B.M_lag(i), Df*(s_tip - B.s(i)), 1e-9);
outl = repmat({xy0}, ns, 1);
[et, eb] = beam_skin_strains(P, B, outl, 'horizontal', false);
check('horizontal model: top = kappa*h/2', et(i), Lf*(s_tip - B.s(i))/S0.EIxx*h/2, 1e-9);
check('horizontal model: bottom = -kappa*h/2', eb(i), -Lf*(s_tip - B.s(i))/S0.EIxx*h/2, 1e-9);
[et2, eb2] = beam_skin_strains(P, B, outl, 'unsymmetric', false);
check('unsymmetric = horizontal (symmetric section), top', et2(i), et(i), 1e-6);
check('unsymmetric = horizontal (symmetric section), bottom', eb2(i), eb(i), 1e-6);
[et3, ~] = beam_skin_strains(P, B, outl, 'unsymmetric', true);
check('drag adds M_lag*b/2/EIzz at the corner', et3(i), et(i) + Df*(s_tip - B.s(i))*(b/2)/S0.EIzz, 1e-6);
fprintf('run_tests_wings done.\n');
