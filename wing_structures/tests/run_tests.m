% RUN_TESTS  Check each building block against a closed-form result.
% Run from the wing_structures folder:  addpath('tests'); run_tests
clear; nfail = 0;

%% 1. CLT: single 0-deg ply recovers ply properties
ply = struct('E1',135e9,'E2',10e9,'G12',5e9,'nu12',0.3,'t',0.125e-3,'rho',1600);
lam = laminate_clt(ply, 0);
nfail = nfail + check('CLT 0-deg Ex = E1',  lam.Ex,  ply.E1,  1e-10);
nfail = nfail + check('CLT 0-deg Gxy = G12',lam.Gxy, ply.G12, 1e-10);

%% 2. CLT: woven ply (E1=E2) at [0 90 90 0] -> Ex = E1 exactly
wov = struct('E1',55e9,'E2',55e9,'G12',4e9,'nu12',0.05,'t',0.25e-3,'rho',1450);
lam = laminate_clt(wov, [0 90 90 0]);
nfail = nfail + check('CLT woven 0/90 Ex = E1', lam.Ex, wov.E1, 1e-10);

%% 3. CLT: woven at [45 -45 -45 45] -> Gxy = E/(2(1+nu12))
lam = laminate_clt(wov, [45 -45 -45 45]);
nfail = nfail + check('CLT woven +/-45 Gxy', lam.Gxy, wov.E1/(2*(1+wov.nu12)), 1e-10);

%% 4. Thin-wall circle: EI = E*pi*r^3*t, GJ = 2*pi*r^3*t*G
r = 0.05; t = 1e-3; E = 70e9; G = 26e9; n = 4000;
th = linspace(0, 2*pi, n+1)'; th(end) = [];
S = thinwall_closed_section([r*cos(th) r*sin(th)], t, E, G, 2700);
nfail = nfail + check('Circle EIxx', S.EIxx, E*pi*r^3*t, 1e-5);
nfail = nfail + check('Circle GJ',   S.GJ,   2*pi*r^3*t*G, 1e-5);
nfail = nfail + check('Circle EIxz ~ 0', S.EIxz/S.EIxx + 1, 1, 1e-8);

%% 5. Rectangular tube vs thin-wall closed form (midline dims)
b = 0.0254; h = 0.0508; t = 1.65e-3;
S = thinwall_closed_section(rect_tube_midline(b, h, t), t, E, G, 2700);
bm = b - t; hm = h - t;
I_tw = t*hm^3/6 + bm*t*hm^2/2;
nfail = nfail + check('Rect tube EIxx (thin-wall)', S.EIxx, E*I_tw, 1e-12);
GJ_tw = 4*(bm*hm)^2 * t * G / (2*(bm + hm));
nfail = nfail + check('Rect tube GJ (Bredt)', S.GJ, GJ_tw, 1e-12);
I_exact = (b*h^3 - (b-2*t)*(h-2*t)^3)/12;
fprintf('   info: thin-wall vs exact solid-formula I differ by %.2f%%\n', ...
        100*(I_tw - I_exact)/I_exact);

%% 6. Modulus weighting: top flange 3x stiffer shifts centroid up
% Rectangle midline bm x hm, all walls t; top flange E_top = 3E, others E.
xz = rect_tube_midline(b, h, t);       % segments: bottom, right, top, left
Ev = [E; E; 3*E; E];
S = thinwall_closed_section(xz, t, Ev, G, 2700);
% Analytic: z-moments of E*A (webs centred at 0)
zc_expect = (3*E*bm*t*(hm/2) + E*bm*t*(-hm/2)) / (E*t*(2*hm + bm) + 3*E*bm*t);
nfail = nfail + check('Modulus-weighted centroid', S.zc, zc_expect, 1e-12);

%% 7. Capped sandwich web GJ == thin-wall box with matching walls
face = laminate_clt(wov, [45 -45 -45 45]);   % one face
capl = laminate_clt(wov, [0 0]);
c = 6e-3; htot = 0.03;
cap = struct('lam', capl, 'w', c + 2*face.h);
Sw = sandwich_web_section(htot, face, c, 48, cap, 0);
d  = c + face.h; hw = htot - 2*capl.h;
% Equivalent box: walls 1,3 = caps (width d), walls 2,4 = faces (height hw+tc)
xzb = [-d/2 -(hw+capl.h)/2; d/2 -(hw+capl.h)/2; d/2 (hw+capl.h)/2; -d/2 (hw+capl.h)/2];
tb  = [capl.h; face.h; capl.h; face.h];
Gb  = [capl.Gxy; face.Gxy; capl.Gxy; face.Gxy];
Sb  = thinwall_closed_section(xzb, tb, 1, Gb, 1);
nfail = nfail + check('Sandwich capped GJ vs box', Sw.GJ, Sb.GJ, 1e-12);

%% 8. Uncapped sandwich plate GJ == Bredt box with rigid edges (limit check)
Sw0 = sandwich_web_section(htot, face, c, 48, [], 0);
Sb0 = thinwall_closed_section([-d/2 -htot/2; d/2 -htot/2; d/2 htot/2; -d/2 htot/2], ...
        [1e3; face.h; 1e3; face.h], 1, face.Gxy, 1);
nfail = nfail + check('Sandwich plate GJ vs stiff-edge box', Sw0.GJ, Sb0.GJ, 1e-3);

fprintf('\n%d test(s) failed.\n', nfail);
