% MAIN_SIZING  Stage 2: size the skin layup, aluminium plate spar and spar
% position for the wing main element.
%
% For every combination of skin layup x plate gauge x spar position:
%   run A  all LOW stiffness      -> tip deflection and twist at LIMIT load
%   run B  plate HIGH, skin LOW   -> plate stress at ULTIMATE (plate takes most)
%   run C  skin HIGH, plate LOW   -> skin strain, rib spacing for skin
%                                    buckling, mount rib load at ULTIMATE
% (Bounding runs because skin and plate share load in proportion to EI.)
%
% Requires results/<name>_geometry.mat (main_geometry) and
% results/<name>_loads.mat (main_loads). Run from wing_structures.
clear; close all;
cfg = wing_config();                  % active wing: see wing_config.m
M   = materials_input();
Gm  = load(fullfile('results', [cfg.name '_geometry.mat']));
Ld  = load(fullfile('results', [cfg.name '_loads.mat'])); Ld = Ld.loads;
z_of_s = @(s) cfg.z_sym + s / cfg.side;
fprintf('Wing: %s\n', cfg.name);
% Skin strain model (see beam_skin_strains); defaults keep the original method
strain_model = 'horizontal'; drag_bending = false;
if isfield(cfg, 'strain_model'), strain_model = cfg.strain_model; end
if isfield(cfg, 'drag_bending'), drag_bending = cfg.drag_bending; end

%% Loads at CFD speed: strips + flap/endplate lump at the tip
L = Ld.L;
nk = numel(L.s_mid);
ld0.s = [L.s_mid; cfg.side * (Ld.r_ep(3) - cfg.z_sym)];   % strips + tip lump
ld0.r = [repmat(L.ref_xy, nk, 1) z_of_s(L.s_mid); Ld.r_ep];
ld0.F = [L.F; Ld.F_ep];
ld0.C = [L.M; Ld.M_ep];
s_tip = ld0.s(end);
scale = @(ld, k) struct('s', ld.s, 'r', ld.r, 'F', k*ld.F, 'C', k*ld.C);
ldL = scale(ld0, Ld.k_limit);
ldU = scale(ld0, Ld.k_ult);

%% Materials (bounds)
tw = {bound_props(M.twill, 1), bound_props(M.twill, 2)};
al = {bound_props(M.alum, 1),  bound_props(M.alum, 2)};
eps_t = M.twill.eps_t(1); eps_c = M.twill.eps_c(1); Fty = M.alum.Fty(1);

x_LE = max(cellfun(@(o) max(o(:,1)), Gm.outlines));
chord_root = Gm.G.chord(1);
gauges = cfg.plate_gauges_in;

opt = struct('s_m', cfg.mount_s, 's_tip', s_tip, 'z_of_s', z_of_s, 'side', cfg.side, ...
             'bc', cfg.mount_bc, 'ns', 700);
R = [];
fprintf('Sweeping %d skin x %d spar positions x %d gauges ...\n', numel(cfg.skin_options), ...
    numel(cfg.spar_aft_LE), numel(gauges));
for j = 1:numel(cfg.skin_options)
    lay = cfg.skin_options{j};
    lamL = laminate_clt(tw{1}, lay); lamH = laminate_clt(tw{2}, lay);
    for xi = 1:numel(cfg.spar_aft_LE)
        xp = x_LE - cfg.spar_aft_LE(xi);
        for gi = 1:numel(gauges)
            t = gauges(gi) * 0.0254;
            mk = @(a) struct('t', t, 'E', a.E, 'G', a.G, 'rho', a.rho);
            PA = station_props(Gm.outlines, Gm.G.s, lamL, mk(al{1}), xp, cfg.bond_gap);
            PB = station_props(Gm.outlines, Gm.G.s, lamL, mk(al{2}), xp, cfg.bond_gap);
            PC = station_props(Gm.outlines, Gm.G.s, lamH, mk(al{1}), xp, cfg.bond_gap);
            if any(PA.h_p <= 0), continue; end
            om = opt;
            om.r_mount_xy = [xp interp1(PA.s, PA.y_p, cfg.mount_s, 'linear', 'extrap')];

            BA = halfspan_beam(PA, ldL, om);
            BB = halfspan_beam(PB, ldU, om);
            BC = halfspan_beam(PC, ldU, om);
            hold_i = @(P, f, s) interp1(P.s, P.(f), min(max(s, P.s(1)), P.s(end)));

            % plate stress (run B)
            sig_p = al{2}.E * abs(BB.kappa) .* hold_i(PB, 'h_p', BB.s) / 2;
            % skin strains (run C): top fibre = +kappa*(y_top-y_c) for hogging
            [e_top, e_bot] = beam_skin_strains(PC, BC, Gm.outlines, strain_model, drag_bending);
            e_t = max([e_top; e_bot; 0]);  e_c = max([-e_top; -e_bot; 0]);
            ecomp = max(-e_top, -e_bot); ecomp(ecomp < 1e-12) = 1e-12;
            a_rib = pi * sqrt(lamH.D(1,1) ./ (lamH.h * lamH.Ex * ecomp));   % wide-column, flat (conservative)
            tau = abs(BC.T) ./ (2 * hold_i(PC, 'A_enc', BC.s) * lamH.h);
            % mount rib: skin's share of the mount reaction must pass through ribs into the plate
            im = find(PC.s >= cfg.mount_s, 1); if isempty(im), im = numel(PC.s); end
            share_s = interp1(PC.s, PC.EI_s ./ PC.EI, cfg.mount_s, 'linear', 'extrap');
            % mass, half wing (stations held to the ends)
            sg = linspace(0, s_tip, 400)';
            m_s = trapz(sg, hold_i(PA, 'm_s', sg)); m_p = trapz(sg, hold_i(PA, 'm_p', sg));

            r.layup = mat2str(lay); r.nply = numel(lay); r.t_skin_mm = 1e3*lamL.h;
            r.spar_aft_mm = 1e3*cfg.spar_aft_LE(xi); r.spar_xc = cfg.spar_aft_LE(xi)/chord_root;
            r.t_in = gauges(gi); r.t_mm = 1e3*t;
            r.h_p_min_mm = 1e3*min(PA.h_p);
            r.m_skin = m_s; r.m_plate = m_p; r.m_tot = m_s + m_p;
            r.w_tip_mm = 1e3*BA.w(end); r.w_max_mm = 1e3*max(abs(BA.w));
            r.twist_tip_deg = rad2deg(BA.phi(end));
            r.M_mount_ult = BB.M_hog(BB.im);
            r.sig_p_MPa = max(sig_p)/1e6; r.MS_plate = Fty/max(sig_p) - 1;
            r.e_t = e_t; r.e_c = e_c; r.MS_skin = min(eps_t/max(e_t,eps), eps_c/max(e_c,eps)) - 1;
            r.rib_max_mm = 1e3*min(a_rib);
            r.tau_skin_MPa = max(tau)/1e6;
            r.V_mount_rib = share_s * abs(BC.R(2)); r.T_mount = abs(BC.CRz);
            r.skin_share = share_s;
            r.pass = r.MS_plate >= 0 && r.MS_skin >= 0;
            R = [R; r]; %#ok<AGROW>
        end
    end
end

%% Write full table
fid = fopen(fullfile('results', [cfg.name '_sizing.csv']), 'w');
f = fieldnames(R);
fprintf(fid, '%s\n', strjoin(f', ','));
for i = 1:numel(R)
    v = cellfun(@(q) R(i).(q), f, 'UniformOutput', false);
    v = cellfun(@(q) tern_str(q), v, 'UniformOutput', false);
    fprintf(fid, '%s\n', strjoin(v', ','));
end
fclose(fid);

%% Summary
fprintf('\nLoads (half wing): limit factor %.3f, ultimate factor %.3f (n = %.2f, SF = %.2f; PLACEHOLDERS)\n', ...
    Ld.k_limit, Ld.k_ult, cfg.n_load, cfg.SF);
fprintf('Mount at s = %.2f mm (%s). Bending moment at mount, ultimate: %.1f N m (statically determinate)\n', ...
    1e3*cfg.mount_s, cfg.mount_bc, R(1).M_mount_ult);
fprintf('Material properties are PLACEHOLDERS (materials_input.m) - treat numbers as structure, not answers.\n');
for j = 1:numel(cfg.skin_options)
    lay = mat2str(cfg.skin_options{j});
    fprintf('\nSkin %s (%.2f mm):\n', lay, R(find(strcmp({R.layup}, lay), 1)).t_skin_mm);
    fprintf('  %8s %6s %8s %8s %9s %9s %10s %8s %9s %9s %10s\n', 'spar aft', 'x/c', 'plate', 'h_min', ...
        'mass[g]', 'tip w[mm]', 'twist[deg]', 'MS plt', 'MS skin', 'rib max', 'skin share');
    for xi = 1:numel(cfg.spar_aft_LE)
        sel = find(strcmp({R.layup}, lay) & abs([R.spar_aft_mm] - 1e3*cfg.spar_aft_LE(xi)) < 1e-6 & [R.pass]);
        if isempty(sel)
            fprintf('  %6.0fmm  -- no gauge passes --\n', 1e3*cfg.spar_aft_LE(xi)); continue;
        end
        [~, b] = min([R(sel).m_tot]); q = R(sel(b));
        fprintf('  %6.0fmm %6.3f %6.3fin %6.1fmm %9.0f %9.2f %10.3f %8.2f %9.2f %7.0fmm %10.2f\n', ...
            q.spar_aft_mm, q.spar_xc, q.t_in, q.h_p_min_mm, 1e3*q.m_tot, q.w_tip_mm, q.twist_tip_deg, ...
            q.MS_plate, q.MS_skin, q.rib_max_mm, q.skin_share);
    end
end
fprintf('\n(Lightest passing plate gauge per row. mass = skin + plate, half wing, no ribs/bond/mounts.\n');
fprintf(' tip w / twist at LIMIT load; MS = margin of safety at ULTIMATE; rib max = largest rib\n');
fprintf(' spacing before the compression skin buckles as a flat wide column (curvature ignored).\n');
fprintf(' skin share = fraction of bending carried by the skin at the mount.)\n');

%% Plots: tip deflection and twist vs gauge, for one spar position
xi0 = find(abs(cfg.spar_aft_LE - 0.080) < 1e-6, 1); if isempty(xi0), xi0 = ceil(numel(cfg.spar_aft_LE)/2); end
figure;
for j = 1:numel(cfg.skin_options)
    lay = mat2str(cfg.skin_options{j});
    sel = strcmp({R.layup}, lay) & abs([R.spar_aft_mm] - 1e3*cfg.spar_aft_LE(xi0)) < 1e-6;
    subplot(2,1,1); hold on; plot([R(sel).t_mm], [R(sel).w_tip_mm], 'o-', 'DisplayName', ['skin ' lay]);
    subplot(2,1,2); hold on; plot([R(sel).t_mm], [R(sel).twist_tip_deg], 'o-', 'DisplayName', ['skin ' lay]);
end
subplot(2,1,1); grid on; ylabel('Tip deflection [mm]'); legend('Location', 'northeast');
title(sprintf('Limit load, spar %.0f mm aft of LE', 1e3*cfg.spar_aft_LE(xi0)));
subplot(2,1,2); grid on; ylabel('Tip twist [deg]'); xlabel('Plate thickness [mm]');
print('-dpng', fullfile('results', [cfg.name '_sizing_stiffness.png']));

% Diagrams for a baseline case
lam0 = laminate_clt(tw{1}, cfg.skin_options{min(2, end)});
P0 = station_props(Gm.outlines, Gm.G.s, lam0, struct('t', 0.125*0.0254, 'E', al{1}.E, 'G', al{1}.G, ...
    'rho', al{1}.rho), x_LE - cfg.spar_aft_LE(xi0), cfg.bond_gap);
o0 = opt; o0.r_mount_xy = [x_LE - cfg.spar_aft_LE(xi0) interp1(P0.s, P0.y_p, cfg.mount_s, 'linear', 'extrap')];
B0 = halfspan_beam(P0, ldU, o0);
figure;
subplot(2,1,1); plot(1e3*B0.s, B0.M_hog); grid on; ylabel('Bending moment [N m]');
title(sprintf('Ultimate load: skin %s, 0.125 in plate, spar %.0f mm aft of LE', ...
    mat2str(cfg.skin_options{min(2,end)}), 1e3*cfg.spar_aft_LE(xi0)));
subplot(2,1,2); plot(1e3*B0.s, B0.T); grid on; ylabel('Torque about elastic axis [N m]');
xlabel('Span from centre plane [mm]');
print('-dpng', fullfile('results', [cfg.name '_diagrams.png']));
% Rib spacing needed along the span and mount loads, for each skin, at the
% baseline spar position with the thinnest plate (ultimate, run C bounds)
fprintf('\nAllowable rib spacing along the span [mm] (flat wide-column, ultimate), spar %.0f mm aft of LE:\n', 1e3*cfg.spar_aft_LE(xi0));
sq = [0.05 0.15 0.2069 0.25 0.35 0.45 0.55 0.62];
fprintf('  %-14s', 'skin'); fprintf('%7.0f', 1e3*sq); fprintf('   | mount rib V [N]  mount T [N m]\n');
for j = 1:numel(cfg.skin_options)
    lay = cfg.skin_options{j};
    lamH = laminate_clt(tw{2}, lay);
    PCj = station_props(Gm.outlines, Gm.G.s, lamH, struct('t', gauges(1)*0.0254, 'E', al{1}.E, 'G', al{1}.G, ...
        'rho', al{1}.rho), x_LE - cfg.spar_aft_LE(xi0), cfg.bond_gap);
    oj = opt; oj.r_mount_xy = [x_LE - cfg.spar_aft_LE(xi0) interp1(PCj.s, PCj.y_p, cfg.mount_s, 'linear', 'extrap')];
    BCj = halfspan_beam(PCj, ldU, oj);
    hi = @(f, s) interp1(PCj.s, PCj.(f), min(max(s, PCj.s(1)), PCj.s(end)));
    ec = max(-BCj.kappa .* (hi('y_top', BCj.s) - hi('y_c', BCj.s)), BCj.kappa .* (hi('y_c', BCj.s) - hi('y_bot', BCj.s)));
    ec = max(ec, 1e-12);
    a = pi * sqrt(lamH.D(1,1) ./ (lamH.h * lamH.Ex * ec));
    sel = strcmp({R.layup}, mat2str(lay)) & abs([R.spar_aft_mm] - 1e3*cfg.spar_aft_LE(xi0)) < 1e-6 & [R.t_in] == gauges(1);
    fprintf('  %-14s', mat2str(lay)); fprintf('%7.0f', 1e3*interp1(BCj.s, a, sq)); 
    fprintf('   | %10.0f %14.1f\n', R(sel).V_mount_rib, R(sel).T_mount);
end
fprintf('(mount rib V = skin share of the mount reaction, carried by ribs into the plate;\n mount T = pitch torque the mount reacts, ultimate)\n');
save(fullfile('results', [cfg.name '_sizing.mat']), 'R', '-v7');
