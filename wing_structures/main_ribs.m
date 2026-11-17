% MAIN_RIBS  Rib layout for the chosen design, and how extra plies change it.
%
% For each skin variant in cfg.rib_variants (same spar plate):
%   1. Beam at LIMIT (low stiffness)      -> tip deflection, twist
%      Beam at ULTIMATE (skin high bound) -> skin force N(s), skin strain
%      Beam at ULTIMATE (plate high)      -> plate stress
%   2. Allowable rib spacing from skin buckling between ribs:
%         a_max(s) = (pi / K) * sqrt(D11 / N_comp)   (flat wide column;
%         K = cfg.rib_K, 1 for pinned ribs; N at cfg.rib_check_load.
%         Ignores skin curvature, so conservative)
%   3. Loft corners: where the skin turns in front view, the skin force
%      pushes sideways with N*dtheta per unit chord. Corners above
%      cfg.kink_rib_min get their own rib.
%   4. Ribs: mount + corner ribs fixed, the rest spaced so that no bay is
%      longer than a_max anywhere inside it.
%
% Needs results from main_geometry and main_loads. Run from wing_structures.
clear; close all;
cfg = front_wing_config();
M   = materials_input();
Gm  = load(fullfile('results', [cfg.name '_geometry.mat']));
Ld  = load(fullfile('results', [cfg.name '_loads.mat'])); Ld = Ld.loads;
z_of_s = @(s) cfg.z_sym + s / cfg.side;
D = cfg.rib_design;
% Lofted regions (wing rises or drops in front view): no ribs inside them when cfg.rib_no_loft
no_loft = isfield(cfg, 'rib_no_loft') && cfg.rib_no_loft;
lofts = zeros(0, 2);
if no_loft
    lofts = loft_regions(Gm.front.s, Gm.front.y_top, Gm.front.y_bot, cfg.loft_slope_min);
end

L = Ld.L; nk = numel(L.s_mid);
ld0.s = [L.s_mid; cfg.side * (Ld.r_ep(3) - cfg.z_sym)];
ld0.r = [repmat(L.ref_xy, nk, 1) z_of_s(L.s_mid); Ld.r_ep];
ld0.F = [L.F; Ld.F_ep]; ld0.C = [L.M; Ld.M_ep];
s_tip = ld0.s(end);
scale = @(k) struct('s', ld0.s, 'r', ld0.r, 'F', k*ld0.F, 'C', k*ld0.C);
ldL = scale(Ld.k_limit); ldU = scale(Ld.k_ult);

tw = {bound_props(M.twill, 1), bound_props(M.twill, 2)};
al = {bound_props(M.alum, 1),  bound_props(M.alum, 2)};
x_LE = max(cellfun(@(o) max(o(:,1)), Gm.outlines));
xp = x_LE - D.spar_aft;
plate = @(a) struct('t', D.plate_in*0.0254, 'E', a.E, 'G', a.G, 'rho', a.rho);
opt = struct('s_m', cfg.mount_s, 's_tip', s_tip, 'z_of_s', z_of_s, 'side', cfg.side, ...
             'bc', cfg.mount_bc, 'ns', 1300);
hold_i = @(P, f, s) interp1(P.s, P.(f), min(max(s, P.s(1)), P.s(end)));

nv = numel(cfg.rib_variants);
V = struct('name', {}, 'ribs', {}, 'bays', {}, 'a', {}, 's', {}, 'corners', {});
for v = 1:nv
    var = cfg.rib_variants{v};
    zones = var.zones;
    PA = zoned_props(Gm.outlines, Gm.G.s, tw{1}, D.base_layup, zones, plate(al{1}), xp, cfg.bond_gap);
    PB = zoned_props(Gm.outlines, Gm.G.s, tw{1}, D.base_layup, zones, plate(al{2}), xp, cfg.bond_gap);
    PC = zoned_props(Gm.outlines, Gm.G.s, tw{2}, D.base_layup, zones, plate(al{1}), xp, cfg.bond_gap);
    o = opt; o.r_mount_xy = [xp interp1(PA.s, PA.y_p, cfg.mount_s, 'linear', 'extrap')];
    BA = halfspan_beam(PA, ldL, o);
    BB = halfspan_beam(PB, ldU, o);
    BC = halfspan_beam(PC, ldU, o);
    sg = BC.s;

    % skin laminate (high bound) at every grid point, by zone
    D11 = zeros(size(sg)); Ext = D11; ar = D11;
    lamH = @(lay) laminate_clt(tw{2}, lay); lamL = @(lay) laminate_clt(tw{1}, lay);
    for i = 1:numel(sg)
        lay = D.base_layup;
        for z = 1:numel(zones)
            if sg(i) >= zones(z).s(1) && sg(i) <= zones(z).s(2), lay = zones(z).layup; end
        end
        lh = lamH(lay); ll = lamL(lay);
        D11(i) = lh.D(1,1); Ext(i) = lh.Ex * lh.h; ar(i) = ll.areal;
    end
    yt = hold_i(PC, 'y_top', sg); yb = hold_i(PC, 'y_bot', sg); yc = hold_i(PC, 'y_c', sg);
    N_top =  Ext .* BC.kappa .* (yt - yc);          % + tension
    N_bot = -Ext .* BC.kappa .* (yc - yb);
    N_c = max(max(-N_top, -N_bot), 1e-6);           % compression at ultimate
    load_fac = 1;
    if strcmpi(cfg.rib_check_load, 'limit'), load_fac = Ld.k_limit / Ld.k_ult; end
    if strcmpi(cfg.rib_check_load, 'aero'),  load_fac = Ld.k_speed / Ld.k_ult; end
    a_ult1 = pi * sqrt(D11 ./ N_c);                 % pinned, ultimate (reference)
    a_max = (pi / cfg.rib_K) * sqrt(D11 ./ (load_fac * N_c));

    % strength
    e_top = BC.kappa .* (yt - yc); e_bot = -BC.kappa .* (yc - yb);
    MS_skin = min(M.twill.eps_t(1) / max([e_top; e_bot]), M.twill.eps_c(1) / max([-e_top; -e_bot])) - 1;
    sig_p = al{2}.E * abs(BB.kappa) .* hold_i(PB, 'h_p', BB.s) / 2;
    MS_plate = M.alum.Fty(1) / max(sig_p) - 1;

    % mass (half wing): perimeter x areal density by zone, plus plate
    perim = hold_i(PA, 'm_s', sg) ./ hold_i(PA, 'areal', sg);
    m_skin = trapz(sg, perim .* ar); m_plate = trapz(sg, hold_i(PA, 'm_p', sg));

    % loft corners
    Nt_f = interp1(sg, N_top, Gm.front.s, 'linear', 'extrap');
    Nb_f = interp1(sg, N_bot, Gm.front.s, 'linear', 'extrap');
    K = kink_corners(Gm.front.s, Gm.front.y_top, Gm.front.y_bot, Nt_f, Nb_f, cfg.kink_curv_min);
    fixed = cfg.mount_s;
    for k = 1:numel(K)
        K(k).needs_rib = max(K(k).F_top, K(k).F_bot) >= cfg.kink_rib_min;
        if ~no_loft && K(k).needs_rib && min(abs(K(k).s - fixed)) > cfg.rib_merge
            fixed(end+1) = K(k).s; %#ok<SAGROW>
        end
    end
    % no ribs inside a loft: a rib just outside each end of each loft instead
    % (an end within cfg.rib_merge of an existing fixed rib uses that rib)
    for j = 1:size(lofts, 1)
        for e = [lofts(j,1) - cfg.rib_loft_margin, lofts(j,2) + cfg.rib_loft_margin]
            if e > 0 && e < s_tip - cfg.rib_merge && min(abs(e - fixed)) > cfg.rib_merge
                fixed(end+1) = e; %#ok<SAGROW>
            end
        end
    end

    [ribs, bays] = rib_layout(sg, a_max, fixed, s_tip);
    % ribs that still fell inside a loft: the skin there cannot span the loft at this spacing
    in_loft = false(size(ribs));
    for j = 1:size(lofts, 1)
        in_loft = in_loft | (ribs > lofts(j,1) + 1e-9 & ribs < lofts(j,2) - 1e-9);
    end
    for k = 1:numel(K), K(k).d_rib = min(abs(K(k).s - ribs)); end
    ri = ribs(ribs < s_tip - 1e-6);                  % ribs on this half (endplate excluded)
    m_ribs = sum(hold_i(PA, 'A_enc', ri) * cfg.rib_areal + hold_i(PA, 'm_s', ri) ./ hold_i(PA, 'areal', ri) * cfg.rib_fillet);
    % sensitivity: rib count for other fixity / check-load assumptions
    Ks = [1.0 0.7 0.5]; lf = [1, Ld.k_limit/Ld.k_ult, Ld.k_speed/Ld.k_ult]; nsens = zeros(numel(Ks), 3);
    for i = 1:numel(Ks)
        for j = 1:3
            r_ij = rib_layout(sg, a_ult1 / Ks(i) / sqrt(lf(j)), fixed, s_tip);
            nsens(i,j) = nnz(r_ij < s_tip - 1e-6);
        end
    end
    V(v).nsens = nsens;
    V(v).lofts = lofts; V(v).ribs_in_loft = ribs(in_loft);
    V(v).name = var.name; V(v).ribs = ribs; V(v).bays = bays; V(v).a = a_max; V(v).s = sg;
    V(v).corners = K; V(v).fixed = sort(fixed);
    V(v).m = m_skin + m_plate; V(v).m_ribs = m_ribs; V(v).w = 1e3*BA.w(end); V(v).th = rad2deg(BA.phi(end));
    V(v).MS = min(MS_skin, MS_plate);
    V(v).n_ribs = nnz(ribs < s_tip - 1e-6);          % one half, excludes the tip (endplate)
    V(v).n_full = 2*V(v).n_ribs - any(abs(ribs) < 1e-9);   % whole wing (a centre-plane rib counts once)
end

%% Report
fprintf('Design: base skin %s, %.3f in 6061-T6 plate %.0f mm behind LE, mount at %.2f mm (%s)\n', ...
    mat2str(D.base_layup), D.plate_in, 1e3*D.spar_aft, 1e3*cfg.mount_s, cfg.mount_bc);
fprintf('Loads: limit factor %.3f, ultimate %.3f\n', Ld.k_limit, Ld.k_ult);

K = V(1).corners;
fprintf('\nLoft corners (%s, ultimate):\n', V(1).name);
fprintf('%9s %16s %9s %13s %13s %10s\n', 's [mm]', 'turning [mm]', 'angle', 'F_top [N/mm]', 'F_bot [N/mm]', 'own rib?');
for k = 1:numel(K)
    yn = 'no';
    if K(k).needs_rib
        if no_loft, yn = sprintf('%.0f mm off', 1e3*K(k).d_rib);
        elseif abs(K(k).s - cfg.mount_s) <= cfg.rib_merge, yn = 'mount rib'; else, yn = 'yes'; end
    end
    fprintf('%9.0f %8.0f-%-7.0f %8.1fd %13.2f %13.2f %10s\n', 1e3*K(k).s, 1e3*K(k).s_range, ...
        rad2deg(K(k).dtheta), K(k).F_top/1e3, K(k).F_bot/1e3, yn);
end
fprintf('(F = skin force x angle change: the sideways push per mm of chord at the corner)\n');
if no_loft
    fprintf('No ribs inside lofted regions (cfg.rib_no_loft): %s mm. A rib sits just outside each end.\n', ...
        strjoin(arrayfun(@(j) sprintf('%.0f-%.0f', 1e3*lofts(j,1), 1e3*lofts(j,2)), 1:size(lofts,1), 'UniformOutput', false), ', '));
    fprintf('("own rib?" = distance from a loaded corner to the nearest rib; the skin carries the corner load over it - check in FEA)\n');
end

fprintf('\n%-30s %6s %6s %13s %9s %9s %10s %10s %7s\n', 'Variant (masses: half wing)', 'ribs/', 'ribs/', 'skin+plate[g]', 'ribs[g]', 'total[g]', 'tip w[mm]', 'twist[deg]', 'MS min');
fprintf('%-30s %6s %6s\n', '', 'half', 'wing');
for v = 1:nv
    fprintf('%-30s %6d %6d %13.0f %9.0f %9.0f %10.2f %10.3f %7.2f\n', V(v).name, V(v).n_ribs, V(v).n_full, 1e3*V(v).m, ...
        1e3*V(v).m_ribs, 1e3*(V(v).m + V(v).m_ribs), V(v).w, V(v).th, V(v).MS);
end
fprintf('(rib mass uses cfg.rib_areal and cfg.rib_fillet - placeholders until the rib material is chosen)\n');

% thicker skin vs more ribs: compare each variant with the first
fprintf('\nThicker skin vs more ribs (whole wing, compared with "%s"):\n', V(1).name);
fprintf('%-30s %13s %12s %22s\n', 'Variant', 'extra skin[g]', 'ribs saved', 'break-even rib mass [g]');
for v = 2:nv
    dm = 2*1e3*(V(v).m - V(1).m); dn = V(1).n_full - V(v).n_full;
    if dn > 0, be = sprintf('%.0f', dm/dn); else, be = 'n/a'; end
    fprintf('%-30s %13.0f %12d %22s\n', V(v).name, dm, dn, be);
end
fprintf(['(Thicker skin is lighter overall only if one rib, with its bonds, weighs more than the\n' ...
         ' break-even mass. Current rib estimate: %.0f g each.)\n'], 1e3*V(1).m_ribs/max(V(1).n_ribs,1));
for v = 1:nv
    fprintf('\n%s - ribs at [mm]: %s, tip %.0f\n', V(v).name, ...
        strjoin(arrayfun(@(x) sprintf('%.0f', 1e3*x), V(v).ribs(V(v).ribs < s_tip - 1e-6)', 'UniformOutput', false), ', '), 1e3*s_tip);
    b = V(v).bays;
    fprintf('  bay lengths [mm]: %s\n', strjoin(arrayfun(@(k) sprintf('%.0f', 1e3*b(k,3)), 1:size(b,1), 'UniformOutput', false), ', '));
    if ~isempty(V(v).ribs_in_loft)
        fprintf('  ** NOT VALID with no ribs in the loft: this skin cannot span it (ribs fell inside at %s mm). Thicken the skin over the loft.\n', ...
            strjoin(arrayfun(@(x) sprintf('%.0f', 1e3*x), V(v).ribs_in_loft', 'UniformOutput', false), ', '));
    end
end
fprintf('\nRib count sensitivity (ONE half wing, excluding endplate; whole wing is about double):\n');
fprintf('%-30s %21s %21s %21s\n', '', 'K=1.0 pinned', 'K=0.7 partly fixed', 'K=0.5 clamped');
fprintf('%-30s %7s %6s %6s %7s %6s %6s %7s %6s %6s\n', 'check load ->', 'ult', 'limit', 'aero', 'ult', 'limit', 'aero', 'ult', 'limit', 'aero');
for v = 1:nv
    n = V(v).nsens;
    fprintf('%-30s %7d %6d %6d %7d %6d %6d %7d %6d %6d\n', V(v).name, n(1,1), n(1,2), n(1,3), n(2,1), n(2,2), n(2,3), n(3,1), n(3,2), n(3,3));
end
fprintf('(ult = %.2f x CFD, limit = %.2f x CFD, aero = %.2f x CFD: max-speed aero load with no load or safety factor)\n', ...
    Ld.k_ult, Ld.k_limit, Ld.k_speed);
fprintf('(current settings: K = %.2f, check at %s)\n', cfg.rib_K, cfg.rib_check_load);
fprintf(['\nNotes: rib spacing uses a flat wide-column buckling estimate (skin curvature and rib\n' ...
         'fixity ignored) - conservative. A buckling FEA of the real skin can relax it.\n' ...
         'Rib count excludes the endplate at the tip; a rib on the centre plane counts as one.\n']);

%% Save
fid = fopen(fullfile('results', [cfg.name '_ribs.csv']), 'w');
fprintf(fid, 'variant,skin_plate_g,ribs_g,tip_defl_mm,tip_twist_deg,MS_min,n_ribs,rib_positions_mm\n');
for v = 1:nv
    fprintf(fid, '"%s",%.0f,%.0f,%.3f,%.4f,%.3f,%d,"%s"\n', V(v).name, 1e3*V(v).m, 1e3*V(v).m_ribs, V(v).w, V(v).th, V(v).MS, ...
        V(v).n_ribs, strjoin(arrayfun(@(x) sprintf('%.1f', 1e3*x), V(v).ribs', 'UniformOutput', false), ' '));
end
fclose(fid);
save(fullfile('results', [cfg.name '_ribs.mat']), 'V', '-v7');

figure; hold on; grid on;
mk = {'-', '--', ':', '-.', '-', '--'};
for v = 1:nv
    plot(1e3*V(v).s, 1e3*min(V(v).a, 0.3), mk{1+mod(v-1,6)}, 'LineWidth', 1.5, 'DisplayName', ['allowable: ' V(v).name]);
end
for v = 1:nv
    r = V(v).ribs;
    plot(1e3*r, (70 - 10*v)*ones(size(r)), 'v', 'DisplayName', ['ribs: ' V(v).name]);
end
xlabel('Span from centre plane [mm]'); ylabel('Allowable rib spacing [mm] (capped at 300)');
title('Rib spacing from skin buckling (ultimate) and resulting rib layouts');
legend('Location', 'northwest');
print('-dpng', fullfile('results', [cfg.name '_ribs.png']));
