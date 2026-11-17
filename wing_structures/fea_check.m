% FEA_CHECK  Prediction for an independent SolidWorks Simulation check of
% the beam model (halfspan_beam + section properties).
%
% Test case (deliberately simple, so both models see exactly the same thing):
%   - skin only, isotropic 6061-T6 shell, t = 1.0 mm (laminate theory is
%     checked separately by the unit tests; this checks the beam model)
%   - overhang only: skin edge FIXED at the mount plane (s = 206.85 mm)
%   - 5 planar ribs (same 1 mm shell) at s = 275, 375, 475, 575, 655 mm
%   - each rib carries one vertical force (limit load of the main element
%     band around it; the tip rib also carries the flap + endplate force),
%     applied as a total force on the rib face -> acts at the rib's area
%     centroid
% Load case B (torsion): only a 50 N m nose-up torque on the tip rib face
% (pushes the leading edge up, trailing edge down).
% Compare FEA vertical displacement at the leading and trailing edge of
% the ribs with the prediction printed here.
%
% Writes results/fea_check_loads.csv (forces and points in CFD AND
% SolidWorks part coordinates) and results/fea_check_prediction.csv.
clear;
cfg = front_wing_config();
Gm  = load(fullfile('results', [cfg.name '_geometry.mat']));
Ld  = load(fullfile('results', [cfg.name '_loads.mat'])); Ld = Ld.loads;
[T, ~, ginfo] = load_main_element(cfg);
z_of_s = @(s) cfg.z_sym + s / cfg.side;
to_part = @(p) [(cfg.stl_t(1) - p(:,1)) p(:,2) -p(:,3)] * 1e3;   % CFD [m] -> STL part [mm]

E = 69e9; nu = 0.33; t_skin = 1.0e-3;
ply = struct('E1', E, 'E2', E, 'G12', E/(2*(1+nu)), 'nu12', nu, 't', t_skin, 'rho', 2700);
skin = laminate_clt(ply, 0);
noplate = struct('t', 1e-9, 'E', 1, 'G', 1, 'rho', 0);          % skin only

s_rib = [0.275 0.375 0.475 0.575 ginfo.s_max];
edges = [cfg.mount_s, (s_rib(1:end-1) + s_rib(2:end))/2, ginfo.s_max + 1e-6];
L = Ld.L; w_lim = Ld.k_limit;

nr = numel(s_rib);
Fy = zeros(nr,1); C = zeros(nr,2); LE = zeros(nr,2); TE = zeros(nr,2);
for k = 1:nr
    in = L.s_mid > edges(k) & L.s_mid <= edges(k+1);
    Fy(k) = w_lim * sum(L.F(in, 2));
    if k == nr, Fy(k) = Fy(k) + w_lim * Ld.F_ep(2); end
    sl = slice_mesh(T, 3, z_of_s(min(s_rib(k), ginfo.s_max - 1e-3)));
    xy = sl{1}(:,1:2);
    [C(k,1), C(k,2)] = polygon_centroid(xy(:,1), xy(:,2));
    [~, i] = max(xy(:,1)); LE(k,:) = xy(i,:);
    [~, i] = min(xy(:,1)); TE(k,:) = xy(i,:);
end

P  = station_props(Gm.outlines, Gm.G.s, skin, noplate, mean(C(:,1)), 0);
ld = struct('s', s_rib(:), 'r', [C z_of_s(s_rib(:))], 'F', [zeros(nr,1) Fy zeros(nr,1)], 'C', zeros(nr,3));
o  = struct('s_m', cfg.mount_s, 's_tip', ginfo.s_max, 'z_of_s', z_of_s, 'side', cfg.side, ...
            'bc', 'clamped', 'ns', 1500, 'r_mount_xy', [P.x_ea(1) P.y_ea(1)]);
B  = halfspan_beam(P, ld, o);

% downward displacement of a point at x: w - phi * (x - x_ea)
fprintf('Skin-only check case: 6061-T6 shell %.1f mm, fixed at s = %.2f mm, limit load\n', 1e3*t_skin, 1e3*cfg.mount_s);
fprintf('\n%6s %9s %24s %24s\n', 'rib', 'force', 'centroid (part, mm)', 'centroid (CFD, mm)');
for k = 1:nr
    pc = to_part([C(k,:) z_of_s(s_rib(k))]);
    fprintf('%4.0fmm %8.1fN   [%6.1f %6.1f %6.1f]   [%6.1f %6.1f %6.1f]\n', 1e3*s_rib(k), -Fy(k), pc, ...
        1e3*[C(k,:) z_of_s(s_rib(k))]);
end
fprintf('Total %.1f N (downward, -Y in both frames)\n', -sum(Fy));

fprintf('\nPrediction (downward displacement, mm; twist nose-up deg):\n');
fprintf('%6s %10s %10s %10s %10s   %s\n', 'rib', 'LE', 'TE', 'centroid', 'twist', 'LE / TE points (part mm)');
pred = zeros(nr, 5);
for k = 1:nr
    w  = interp1(B.s, B.w, s_rib(k)); ph = interp1(B.s, B.phi, s_rib(k)); xe = interp1(B.s, B.x_ea, s_rib(k));
    d  = @(x) w - ph * (x - xe);
    pred(k,:) = [1e3*s_rib(k) 1e3*d(LE(k,1)) 1e3*d(TE(k,1)) 1e3*d(C(k,1)) rad2deg(ph)];
    pl = to_part([LE(k,:) z_of_s(s_rib(k))]); pt = to_part([TE(k,:) z_of_s(s_rib(k))]);
    fprintf('%4.0fmm %10.3f %10.3f %10.3f %10.4f   LE [%.1f %.1f %.1f]  TE [%.1f %.1f %.1f]\n', pred(k,:), pl, pt);
end
fprintf('\nTwist from FEA: (d_TE - d_LE) / (x_LE - x_TE) in radians, nose-up positive.\n');

% Case B: pure torque at the tip
Tq = 50;
ldB = struct('s', ginfo.s_max, 'r', [C(end,:) z_of_s(ginfo.s_max)], 'F', [0 0 0], 'C', [0 0 Tq]);
BB = halfspan_beam(P, ldB, o);
fprintf('\nCase B: %.0f N m nose-up torque on the tip rib only (leading edge pushed UP)\n', Tq);
fprintf('%6s %10s %10s %10s\n', 'rib', 'LE', 'TE', 'twist');
predB = zeros(nr, 4);
for k = 1:nr
    w  = interp1(BB.s, BB.w, s_rib(k)); ph = interp1(BB.s, BB.phi, s_rib(k)); xe = interp1(BB.s, BB.x_ea, s_rib(k));
    d  = @(x) w - ph * (x - xe);
    predB(k,:) = [1e3*s_rib(k) 1e3*d(LE(k,1)) 1e3*d(TE(k,1)) rad2deg(ph)];
    fprintf('%4.0fmm %10.3f %10.3f %10.4f\n', predB(k,:));
end
fprintf('(downward displacement in mm: negative = up)\n');
GJm = interp1(P.s, P.GJ, linspace(cfg.mount_s, ginfo.s_max, 200), 'linear', 'extrap');
fprintf('Hand check: twist = T * integral(ds/GJ) = %.4f deg\n', rad2deg(Tq*trapz(linspace(cfg.mount_s, ginfo.s_max, 200), 1./GJm)));

fid = fopen(fullfile('results', 'fea_check_loads.csv'), 'w');
fprintf(fid, 'rib_s_mm,force_down_N,cx_part_mm,cy_part_mm,cz_part_mm,cx_cfd_mm,cy_cfd_mm,cz_cfd_mm\n');
for k = 1:nr
    fprintf(fid, '%.1f,%.3f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', 1e3*s_rib(k), -Fy(k), ...
        to_part([C(k,:) z_of_s(s_rib(k))]), 1e3*[C(k,:) z_of_s(s_rib(k))]);
end
fclose(fid);
fid = fopen(fullfile('results', 'fea_check_prediction.csv'), 'w');
fprintf(fid, 'rib_s_mm,LE_down_mm,TE_down_mm,centroid_down_mm,twist_noseup_deg\n');
fprintf(fid, '%.1f,%.4f,%.4f,%.4f,%.5f\n', pred');
fprintf(fid, '\ncase_B_torque_%gNm\nrib_s_mm,LE_down_mm,TE_down_mm,twist_noseup_deg\n', Tq);
fprintf(fid, '%.1f,%.4f,%.4f,%.5f\n', predB');
fclose(fid);
