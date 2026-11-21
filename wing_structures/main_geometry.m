% MAIN_GEOMETRY  Slice the main element STL along the span and map the
% depth available for a spar, plus skin section properties per station.
%
% Output per station (s = distance from centreline):
%   chord, chord-line angle, leading-edge height, max thickness and where
%   it is, skin EIxx_eff / GJ / mass from the skin layup in the config.
% Depth map: vertical thickness vs x/c at each station. The spar must fit
% inside this envelope (minus skin and bond allowances) along its path.
% Run from the wing_structures folder.
clear; close all;
cfg = wing_config();                  % active wing: see wing_config.m
M   = materials_input();
if ~exist('results', 'dir'), mkdir('results'); end
fprintf('Wing: %s\n', cfg.name);

[T, ~, ginfo] = load_main_element(cfg);
z_of_s = @(s) cfg.z_sym + s / cfg.side;
st = cfg.stations(cfg.stations > ginfo.s_min & cfg.stations < ginfo.s_max);
ns = numel(st);
nq = 101;

tw   = bound_props(M.twill, 1);                 % low stiffness bound
skin = laminate_clt(tw, cfg.skin_layup);

G = struct('s', st(:), 'chord', zeros(ns,1), 'angle', zeros(ns,1), 'y_LE', zeros(ns,1), ...
    'x_LE', zeros(ns,1), 't_max', zeros(ns,1), 'xc_tmax', zeros(ns,1), 'EI', zeros(ns,1), ...
    'GJ', zeros(ns,1), 'mass', zeros(ns,1), 'nloops', zeros(ns,1), 't_wall_cad', nan(ns,1));
depth = zeros(nq, ns); xc = []; outlines = cell(ns, 1);

for k = 1:ns
    loops = slice_mesh(T, 3, z_of_s(st(k)));
    G.nloops(k) = numel(loops);
    xy = loops{1}(:, 1:2);                       % [x forward, y up], outer mould line
    if numel(loops) >= 2
        % A shelled CAD body gives an inner offset loop; record its wall thickness
        xi = loops{2}(:, 1:2);
        [~, ~, Ao] = polygon_centroid(xy(:,1), xy(:,2));
        [cxi, cyi, Ai] = polygon_centroid(xi(:,1), xi(:,2));
        Po = sum(sqrt(sum(diff([xy; xy(1,:)]).^2, 2)));
        Pi = sum(sqrt(sum(diff([xi; xi(1,:)]).^2, 2)));
        if Pi > 0.9*Po && inpolygon(cxi, cyi, xy(:,1), xy(:,2))
            G.t_wall_cad(k) = (Ao - Ai) / ((Po + Pi)/2);
            G.nloops(k) = numel(loops) - 1;
        end
    end
    outlines{k} = xy;
    S = section_depth(xy, nq);
    xc = S.xc; depth(:, k) = S.t;
    G.chord(k) = S.chord; G.angle(k) = S.angle_deg;
    G.x_LE(k) = S.LE(1); G.y_LE(k) = S.LE(2);
    G.t_max(k) = S.t_max; G.xc_tmax(k) = S.xc_tmax;
    Sk = thinwall_closed_section(xy, skin.h, skin.Ex, skin.Gxy, skin.rho);
    G.EI(k) = Sk.EIxx_eff; G.GJ(k) = Sk.GJ; G.mass(k) = Sk.mass_per_len;
end
if any(G.nloops > 1)
    warning('main_geometry:loops', 'Some slices returned extra loops besides the skin; the longest was used.');
end
if any(isfinite(G.t_wall_cad))
    fprintf('CAD body is shelled: wall thickness %.2f-%.2f mm (outer mould line used).\n', ...
        1e3*min(G.t_wall_cad), 1e3*max(G.t_wall_cad));
end

fprintf('Skin [%s], t = %.2f mm, Ex = %.1f GPa (low bound)\n', num2str(cfg.skin_layup), 1e3*skin.h, skin.Ex/1e9);
fprintf('\n%7s %8s %8s %8s %8s %8s %10s %9s %8s\n', 's[mm]', 'chord', 'angle', 'y_LE', 't_max', 'x/c@t', 'EI[Nm2]', 'GJ[Nm2]', 'kg/m');
for k = 1:ns
    fprintf('%7.1f %8.1f %8.2f %8.1f %8.1f %8.3f %10.0f %9.0f %8.3f\n', 1e3*G.s(k), 1e3*G.chord(k), ...
        G.angle(k), 1e3*G.y_LE(k), 1e3*G.t_max(k), G.xc_tmax(k), G.EI(k), G.GJ(k), G.mass(k));
end
fprintf('(angle = chord line vs horizontal, + nose up; y_LE = leading edge height)\n');

% Front-view skin lines (highest and lowest point of each section) on a
% fine spanwise grid, used to find the loft corners (main_ribs)
sf = (0:0.004:ginfo.s_max)';
front = struct('s', sf, 'y_top', zeros(size(sf)), 'y_bot', zeros(size(sf)));
for i = 1:numel(sf)
    lp = slice_mesh(T, 3, z_of_s(min(max(sf(i), ginfo.s_min + 1e-4), ginfo.s_max - 1e-4)));
    front.y_top(i) = max(lp{1}(:,2)); front.y_bot(i) = min(lp{1}(:,2));
end

save(fullfile('results', [cfg.name '_geometry.mat']), 'G', 'depth', 'xc', 'outlines', 'skin', 'front', '-v7');

%% Plots
figure;
contourf(1e3*st, xc, 1e3*depth, 20); colorbar;
xlabel('Span from centreline [mm]'); ylabel('x/c (0 = leading edge)');
title('Available depth [mm] (outer mould line, before skin/bond allowance)');
set(gca, 'YDir', 'reverse');
print('-dpng', fullfile('results', [cfg.name '_depth_map.png']));

figure; hold on; grid on; axis equal;
for k = 1:ns
    plot(1e3*outlines{k}([1:end 1],1), 1e3*outlines{k}([1:end 1],2));
end
xlabel('x (forward) [mm]'); ylabel('y (up) [mm]');
title('Main element sections at each station (side view)');
print('-dpng', fullfile('results', [cfg.name '_sections.png']));

figure;
subplot(3,1,1); plot(1e3*G.s, 1e3*G.t_max, 'o-'); grid on; ylabel('t_{max} [mm]');
title('Main element along span');
subplot(3,1,2); plot(1e3*G.s, G.xc_tmax, 'o-'); grid on; ylabel('x/c at t_{max}');
subplot(3,1,3); plot(1e3*G.s, G.EI, 'o-'); grid on; ylabel('skin EI [N m^2]');
xlabel('Span from centreline [mm]');
print('-dpng', fullfile('results', [cfg.name '_spanwise.png']));
