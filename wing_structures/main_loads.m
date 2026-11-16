% MAIN_LOADS  Spanwise load distribution on the front wing main element
% from the Fluent wall export, plus the lumped flap + endplate load.
%
% Steps:
%   1. Read STL (to CFD frame) and Fluent faces.
%   2. Check: summed face forces reproduce Fluent's Report > Forces.
%   3. Classify faces: main element = within classify_tol of the STL skin.
%      Everything else (flaps, endplate) is lumped at the main element tip.
%   4. Diagnostics: CFD-to-STL distance by span, and area coverage per
%      strip (classified CFD area / STL skin area; ~1 means complete).
%   5. Integrate main element faces into spanwise strips. Strips with
%      coverage below cfg.coverage_min are interpolated from neighbours,
%      and the flap + endplate lump is recomputed as (whole zone - main
%      element) so total force and moment still match the CFD.
%   6. Scale to design speed, limit and ultimate loads.
% Run from the wing_structures folder.
clear; close all;
cfg = front_wing_config();
if ~exist('results', 'dir'), mkdir('results'); end

%% 1. Read
[T, Tskin, ginfo] = load_main_element(cfg);
F = read_fluent_ascii(cfg.cfd);
s_face = cfg.side * (F.xyz(:,3) - cfg.z_sym);
z_of_s = @(s) cfg.z_sym + s / cfg.side;
fprintf('STL: %d triangles (%d end-cap triangles removed), span %.1f-%.1f mm\n', ...
    size(T,1), ginfo.n_cap, 1e3*ginfo.s_min, 1e3*ginfo.s_max);
fprintf('CFD: %d faces, span %.1f-%.1f mm\n', size(F.xyz,1), 1e3*min(s_face), 1e3*max(s_face));

%% 2. Check integration against Fluent
[Ff, Fp, Fv] = face_forces(F);
Pt = sum(Fp, 1); Vt = sum(Fv, 1);
ep = norm(Pt - cfg.fluent_pressure) / norm(cfg.fluent_pressure);
ev = norm(Vt - cfg.fluent_viscous)  / norm(cfg.fluent_viscous);
fprintf('\nCheck vs Fluent Report > Forces (whole zone):\n');
fprintf('  pressure  [%10.4f %10.4f %10.4f]  Fluent [%10.4f %10.4f %10.4f]  rel err %.1e\n', Pt, cfg.fluent_pressure, ep);
fprintf('  viscous   [%10.4f %10.4f %10.4f]  Fluent [%10.4f %10.4f %10.4f]  rel err %.1e\n', Vt, cfg.fluent_viscous, ev);
if ep > 1e-4 || ev > 1e-4
    warning('main_loads:check', 'Integrated forces do not match Fluent. Check export quantities/units before using results.');
end

%% 3. Classify faces
tic;
dist_wide = distance_to_mesh(F.xyz, Tskin, 6e-3);   % also used for diagnostics
dist = dist_wide;
nlat = abs(F.A(:,3)) ./ F.Amag;
main = dist <= cfg.classify_tol & nlat < cfg.max_lateral_normal & ...
       s_face >= ginfo.s_min - cfg.classify_tol & s_face <= ginfo.s_max + cfg.classify_tol;
fprintf('\nClassification (tol %.1f mm, %.1f s): main element %d faces, other %d faces\n', ...
    1e3*cfg.classify_tol, toc, nnz(main), nnz(~main));

%% 4. Diagnostics: distance and coverage by span
edges = (0:cfg.strip_ds:ceil(ginfo.s_max/cfg.strip_ds)*cfg.strip_ds)';
edges(end) = ginfo.s_max;   % last strip ends at the main element tip (per-span values use true width)
nb = numel(edges) - 1;
% Expected main element surface area per strip: STL section perimeter x
% strip width (slice at the middle of the part of the strip the STL covers)
A_stl = zeros(nb, 1);
for k = 1:nb
    a = edges(k); b = min(edges(k+1), ginfo.s_max);
    if b <= a, continue; end
    lp = slice_mesh(T, 3, z_of_s((a + b)/2));
    if isempty(lp), continue; end
    xy = lp{1}(:,1:2);
    A_stl(k) = sum(sqrt(sum(diff([xy; xy(1,:)]).^2, 2))) * (b - a);
end

fprintf('\nCFD-to-STL distance of faces within 6 mm, by span band:\n');
fprintf('  %-14s %8s %10s %10s\n', 'span [mm]', 'faces', 'median', 'p95 [mm]');
bands = 0:0.1:ceil(10*ginfo.s_max)/10;
for k = 1:numel(bands)-1
    m = isfinite(dist_wide) & s_face >= bands(k) & s_face < bands(k+1) & nlat < cfg.max_lateral_normal;
    if any(m)
        dk = sort(dist_wide(m));
        fprintf('  %4.0f - %-7.0f %8d %10.2f %10.2f\n', 1e3*bands(k), 1e3*bands(k+1), nnz(m), ...
            1e3*median(dk), 1e3*dk(max(1, round(0.95*numel(dk)))));
    end
end

%% 5. Integrate
Ldummy = integrate_strips(F.xyz(main,:), Ff(main,:), F.Amag(main), s_face(main), edges, [0 0], z_of_s);
coverage = Ldummy.area ./ max(A_stl, eps);
low = find(coverage < cfg.coverage_min & A_stl > 0);
if ~isempty(low)
    fprintf(['\nNote: main element faces missing in %d strip(s) (coverage < %.2f), first at s = %.0f mm.\n' ...
        '      These strips are interpolated from their neighbours below. Usual cause: coarse STL\n' ...
        '      or a CAD/CFD geometry difference.\n'], numel(low), cfg.coverage_min, 1e3*Ldummy.s_mid(low(1)));
end

% Moment reference line: through the most forward point of the main
% element, parallel to z (any fixed line works; stage 2 moves moments to
% the spar axis).
[xmax, iv] = max(reshape(T(:,[1 4 7]), [], 1));
Yv = reshape(T(:,[2 5 8]), [], 1);
ref_xy = [xmax, Yv(iv)];
L = integrate_strips(F.xyz(main,:), Ff(main,:), F.Amag(main), s_face(main), edges, ref_xy, z_of_s);

% Fill strips where main element faces were missed
L_raw = L;
[L, filled] = fill_low_coverage(L, coverage, cfg.coverage_min);

% Lumped flap + endplate load about the centre of the main element tip:
% whole zone minus the (corrected) main element, so totals are conserved.
tip_loop = slice_mesh(T, 3, z_of_s(ginfo.s_max - 1e-3));
[cx, cy] = polygon_centroid(tip_loop{1}(:,1), tip_loop{1}(:,2));
r_ep = [cx cy z_of_s(ginfo.s_max)];
F_all = sum(Ff, 1);
M_all = sum(cross(F.xyz - r_ep, Ff, 2), 1);          % whole zone about r_ep
r0 = [repmat(ref_xy, nb, 1) z_of_s(L.s_mid)];       % strip reference points
M_main_ep = sum(L.M + cross(r0 - r_ep, L.F, 2), 1);  % main element about r_ep
F_ep = F_all - sum(L.F, 1);
M_ep = M_all - M_main_ep;
% check: conservation (exact by construction, guards against bookkeeping errors)
e_cons = norm(sum(L.F,1) + F_ep - [cfg.fluent_pressure + cfg.fluent_viscous]) / norm(F_all);
if e_cons > 1e-6
    warning('main_loads:conservation', 'Main + lump does not reproduce the Fluent total (rel err %.1e).', e_cons);
end

k_speed = (cfg.V_design / cfg.V_cfd)^2;
k_limit = k_speed * cfg.n_load;
k_ult   = k_limit * cfg.SF;

%% Report
F_main = sum(L.F, 1);
fprintf('\nMain element (half wing) at CFD speed %.0f m/s:\n', cfg.V_cfd);
fprintf('  F = [%.2f %.2f %.2f] N   (downforce %.2f N, drag %.2f N)\n', F_main, -F_main(2), -F_main(1));
fprintf('Flaps + endplate (half wing), lumped at main element tip [%.1f %.1f %.1f] mm:\n', 1e3*r_ep);
fprintf('  F = [%.2f %.2f %.2f] N,  M = [%.3f %.3f %.3f] N m\n', F_ep, M_ep);
fprintf('Main + other = [%.4f %.4f %.4f] N (should equal the Fluent total)\n', F_main + F_ep);
if any(filled)
    fprintf('Interpolated %d strip(s) with coverage < %.2f: s = %s mm\n', nnz(filled), cfg.coverage_min, ...
        mat2str(round(1e3*L.s_mid(filled))'));
    fprintf('  main element downforce before filling %.2f N, after %.2f N\n', -sum(L_raw.F(:,2)), -F_main(2));
end
fprintf('\nFactors: speed %.3f, limit %.3f (n = %.2f), ultimate %.3f (SF = %.2f)\n', ...
    k_speed, k_limit, cfg.n_load, k_ult, cfg.SF);
fprintf('Half-wing downforce: main %.1f N / flaps+endplate %.1f N at CFD speed; limit %.1f / %.1f N\n', ...
    -F_main(2), -F_ep(2), -k_limit*F_main(2), -k_limit*F_ep(2));

fprintf('\n%8s %10s %10s %10s %9s\n', 's [mm]', 'w [N/m]', 'drag[N/m]', 'm_z[Nm/m]', 'coverage');
tag = {'', '  interpolated'};
for k = 1:nb
    fprintf('%8.1f %10.2f %10.2f %10.3f %9.3f%s\n', 1e3*L.s_mid(k), -L.f(k,2), -L.f(k,1), L.m(k,3), coverage(k), tag{filled(k)+1});
end
fprintf('(w = downforce per unit span; m_z = pitching moment per unit span about the reference line,\n');
fprintf(' + = nose-up about +z (right). Reference line x = %.1f mm, y = %.1f mm.)\n', 1e3*ref_xy);

%% Save
loads = struct('cfg', cfg, 'L', L, 'L_raw', L_raw, 'filled', filled, 'F_ep', F_ep, 'M_ep', M_ep, 'r_ep', r_ep, ...
    'k_speed', k_speed, 'k_limit', k_limit, 'k_ult', k_ult, 'coverage', coverage);
save(fullfile('results', [cfg.name '_loads.mat']), 'loads');
fid = fopen(fullfile('results', [cfg.name '_loads.csv']), 'w');
fprintf(fid, '# %s main element, half wing, CFD speed %.1f m/s. Multiply by k_limit=%.4f or k_ult=%.4f.\n', cfg.name, cfg.V_cfd, k_limit, k_ult);
fprintf(fid, '# Moments per unit span about x=%.6f m, y=%.6f m at each strip centre.\n', ref_xy);
fprintf(fid, '# Flap + endplate lump at [%.6f %.6f %.6f] m: F = [%.6f %.6f %.6f] N, M = [%.6f %.6f %.6f] N m\n', r_ep, F_ep, M_ep);
fprintf(fid, 's_mid_m,ds_m,fx_N_per_m,fy_N_per_m,fz_N_per_m,mx_Nm_per_m,my_Nm_per_m,mz_Nm_per_m,coverage,interpolated\n');
fprintf(fid, '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%d\n', [L.s_mid L.ds L.f L.m coverage double(filled)]');
fclose(fid);

%% Plots
figure;
subplot(2,1,1); hold on; grid on;
stairs([edges(1:end-1); edges(end)]*1e3, [-L_raw.f(:,2); -L_raw.f(end,2)], ':', 'DisplayName', 'CFD, before filling');
stairs([edges(1:end-1); edges(end)]*1e3, [-L.f(:,2); -L.f(end,2)], 'DisplayName', sprintf('%.0f m/s (CFD, filled)', cfg.V_cfd));
stairs([edges(1:end-1); edges(end)]*1e3, k_limit*[-L.f(:,2); -L.f(end,2)], 'DisplayName', 'limit load');
ylabel('Downforce per span [N/m]'); legend('Location', 'best');
title('Front wing main element, half span');
subplot(2,1,2); hold on; grid on;
stairs([edges(1:end-1); edges(end)]*1e3, [coverage; coverage(end)]);
plot([0 edges(end)]*1e3, [1 1], 'k--');
xlabel('Span from centreline [mm]'); ylabel('CFD / STL area');
title('Classification coverage (1 = all main element faces found)');
print('-dpng', fullfile('results', [cfg.name '_loads.png']));
