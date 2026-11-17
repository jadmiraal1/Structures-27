% REVIEW_FIGURES  Design-review figures for one rib variant.
%
% Reads the saved results of main_geometry, main_loads and main_ribs and
% writes three PNG figures to results/review:
%   fig1_loads.png       spanwise downforce at limit load, bending moment at ultimate
%   fig2_section.png     section at the mount: skin, spar plate, neutral axis
%   fig3_rib_layout.png  rib stations on the planform, bay length vs allowable spacing
%
% Run main_ribs first (and again after any change to front_wing_config.m),
% then run this from wing_structures. No sizing is done here: every number
% comes from the results files and the config.
clear; close all;
variant = 'inboard of 350';              % any text contained in the variant name (cfg.rib_variants)
out_dir = fullfile('results', 'review');

cfg = front_wing_config();
M   = materials_input();
Gm  = load(fullfile('results', [cfg.name '_geometry.mat']));
Ld  = load(fullfile('results', [cfg.name '_loads.mat'])); Ld = Ld.loads;
Rb  = load(fullfile('results', [cfg.name '_ribs.mat'])); V = Rb.V;
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

iv = find(arrayfun(@(x) ~isempty(strfind(lower(x.name), lower(variant))), V), 1);
if isempty(iv)
    error('review_figures:variant', 'No variant in %s_ribs.mat contains "%s".', cfg.name, variant);
end
v  = V(iv);
ic = find(cellfun(@(c) strcmp(c.name, v.name), cfg.rib_variants), 1);
if isempty(ic)
    error('review_figures:stale', 'Variant "%s" is not in front_wing_config any more. Re-run main_ribs.', v.name);
end
var = cfg.rib_variants{ic};
D   = cfg.rib_design;

ink = [0.13 0.13 0.13]; blue = [0.12 0.31 0.47]; orange = [0.77 0.35 0.07]; grey = [0.50 0.50 0.50];
tint_thick = [0.86 0.90 0.95]; tint_base = [0.95 0.95 0.95]; tint_loft = [0.78 0.78 0.78];
fs = 12;

%% Common quantities
L = Ld.L; s_mid = L.s_mid(:); ds = L.ds(:);
Fy = L.F(:,2);                                        % strip force, + up [N] at the CFD speed
s_tip = cfg.side * (Ld.r_ep(3) - cfg.z_sym);
ribs = v.ribs(:); bays = v.bays; fixed = v.fixed(:);
lofts = zeros(0, 2);
if isfield(v, 'lofts'), lofts = reshape(v.lofts, [], 2); end

% skin plies along the span for this variant
nb = numel(D.base_layup);
sz = linspace(0, s_tip, 2001)'; npl = nb * ones(size(sz));
for z = 1:numel(var.zones)
    in = sz >= var.zones(z).s(1) & sz <= var.zones(z).s(2);
    npl(in) = numel(var.zones(z).layup);
end
n_mount = interp1(sz, npl, cfg.mount_s, 'nearest');

%% Figure 1: loads
w_lim = -Fy ./ ds * Ld.k_limit;                       % downforce per span at limit [N/m]
edges = [s_mid - ds/2; s_mid(end) + ds(end)/2];
% Bending moment at ultimate from statics (pinned mount, symmetric centre
% plane): the same quantity as halfspan_beam's M_hog for cfg.mount_bc =
% 'pinned'. The small F_z * height term is left out.
sg = linspace(0, s_tip, 600)';
F_down = -(sum(Fy) + Ld.F_ep(2));                     % total downforce on the half wing
Mh = zeros(size(sg));
for i = 1:numel(sg)
    o = s_mid > sg(i);
    Mh(i) = sum(-Fy(o) .* (s_mid(o) - sg(i))) - Ld.F_ep(2) * (s_tip - sg(i)) ...
          + cfg.side * (sum(L.M(o,1)) + Ld.M_ep(1));
    if sg(i) < cfg.mount_s, Mh(i) = Mh(i) - F_down * (cfg.mount_s - sg(i)); end
end
Mh = Mh * Ld.k_ult;
[M_peak, ip] = max(Mh);

f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 640]);
subplot(2, 1, 1); hold on; grid on; box off; set(gca, 'FontSize', fs);
yl = 1.15 * max(w_lim);
stairs(1e3*edges, [w_lim; w_lim(end)], 'Color', blue, 'LineWidth', 2);
plot(1e3*cfg.mount_s*[1 1], [0 yl], '--', 'Color', grey, 'LineWidth', 1.2);
text(1e3*cfg.mount_s + 6, 0.94*yl, 'mount', 'Color', grey, 'FontSize', fs);
text(15, 0.12*yl, sprintf('Main element: %.0f N per half wing', -sum(Fy)*Ld.k_limit), 'Color', blue, 'FontSize', fs);
text(1e3*s_tip - 10, 0.30*yl, sprintf('Flaps + endplate: %.0f N at the tip', -Ld.F_ep(2)*Ld.k_limit), ...
    'Color', orange, 'FontSize', fs, 'HorizontalAlignment', 'right');
xlim([0 1e3*s_tip + 5]); ylim([0 yl]);
ylabel('Downforce per span [N/m]');
title(sprintf('Limit load (CFD x %.2f)', Ld.k_limit));
subplot(2, 1, 2); hold on; grid on; box off; set(gca, 'FontSize', fs);
yl2 = [min(Mh) - 5, 1.15*M_peak];
plot(1e3*sg, Mh, 'Color', blue, 'LineWidth', 2.5);
plot([0 1e3*s_tip + 5], [0 0], 'Color', ink, 'LineWidth', 0.8);
plot(1e3*cfg.mount_s*[1 1], yl2, '--', 'Color', grey, 'LineWidth', 1.2);
text(1e3*sg(ip) + 15, M_peak, sprintf('%.0f N m at the mount', M_peak), 'Color', ink, 'FontSize', fs, ...
    'VerticalAlignment', 'bottom');
xlim([0 1e3*s_tip + 5]); ylim(yl2);
xlabel('Span from centre plane [mm]'); ylabel('Bending moment [N m]');
title(sprintf('Ultimate load (CFD x %.2f); positive = top skin in tension', Ld.k_ult));
print(f1, '-dpng', '-r200', fullfile(out_dir, 'fig1_loads.png'));

%% Figure 2: section at the mount
[~, k] = min(abs(Gm.G.s(:) - cfg.mount_s));
xy = Gm.outlines{k};
p2 = xy([2:end 1], :); Lg = sqrt(sum((p2 - xy).^2, 2)); mid = (xy + p2) / 2;
yc = sum(Lg .* mid(:,2)) / sum(Lg);                   % skin centroid height = neutral axis
xa = 1e3 * (max(xy(:,1)) - xy(:,1));                  % mm behind the leading edge
ya = 1e3 * (xy(:,2) - yc);
chord = max(xa);
x_LE = max(cellfun(@(o) max(o(:,1)), Gm.outlines));   % as in main_ribs
xp = x_LE - D.spar_aft;
Dp = section_depth(xy, 2001);
y_hi = 1e3 * (interp1(Dp.x, Dp.y_hi, xp) - yc); y_lo = 1e3 * (interp1(Dp.x, Dp.y_lo, xp) - yc);
clr = 1e3 * (n_mount * M.twill.t + cfg.bond_gap);     % skin thickness + bond gap
h_plate = (y_hi - y_lo) - 2*clr;
xs = 1e3 * (max(xy(:,1)) - xp);
[y_bot, ib] = min(ya);
it = find(ya > 0 & abs(xa - 0.5*chord) < 4);
y_top = mean(ya(it));

f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 380]);
hold on; grid on; box off; set(gca, 'FontSize', fs);
plot([xa; xa(1)], [ya; ya(1)], 'Color', ink, 'LineWidth', 2.5);
plot([-8 chord + 8], [0 0], '--', 'Color', grey, 'LineWidth', 1.2);
text(0.55*chord, 2.5, 'neutral axis', 'Color', grey, 'FontSize', fs);
patch('XData', xs + [-1.2 1.2 1.2 -1.2], 'YData', [y_lo + clr, y_lo + clr, y_hi - clr, y_hi - clr], ...
    'FaceColor', orange, 'EdgeColor', 'none');
text(xs + 6, -9, {sprintf('Spar plate, %.0f mm behind LE', 1e3*D.spar_aft), sprintf('%.0f mm tall', h_plate)}, ...
    'Color', orange, 'FontSize', fs, 'VerticalAlignment', 'middle');
plot([xa(ib), xa(ib) + 60], [y_bot, y_bot - 9], 'Color', blue, 'LineWidth', 1);
text(xa(ib) + 63, y_bot - 6, {sprintf('Bottom skin: %.1f mm below the neutral axis', -y_bot), '(compression under downforce)'}, ...
    'Color', blue, 'FontSize', fs, 'VerticalAlignment', 'top');
plot([0.5*chord, 0.5*chord + 12], [y_top, y_top + 9], 'Color', blue, 'LineWidth', 1);
text(0.5*chord + 15, y_top + 10, sprintf('Top skin: about %.0f mm above (tension)', y_top), 'Color', blue, 'FontSize', fs);
xlim([-8 chord + 8]); ylim([min(ya) - 22, max(ya) + 20]);
daspect([1 1 1]);
xlabel('Distance behind leading edge [mm]'); ylabel('Height [mm]');
title(sprintf('Section at the mount (%.0f mm station): chord %.0f mm', 1e3*Gm.G.s(k), chord));
print(f2, '-dpng', '-r200', fullfile(out_dir, 'fig2_section.png'));

%% Figure 3: rib layout
others = fixed(abs(fixed - cfg.mount_s) > 1e-6 & abs(fixed - s_tip) > 1e-6);
n_half = nnz(ribs < s_tip - 1e-6);
n_full = 2*n_half - any(abs(ribs) < 1e-9);

f3 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 680]);
subplot(2, 1, 1); hold on; box off; set(gca, 'FontSize', fs, 'YDir', 'reverse');
brk = [1; find(diff(npl) ~= 0) + 1; numel(sz) + 1];
for j = 1:numel(brk) - 1
    z0 = 1e3 * sz(brk(j)); z1 = 1e3 * sz(min(brk(j+1), numel(sz))); n = npl(brk(j));
    col = tint_base; if n > nb, col = tint_thick; end
    patch('XData', [z0 z1 z1 z0], 'YData', [0 0 chord chord], 'FaceColor', col, 'EdgeColor', 'none');
    text((z0 + z1)/2, chord + 50, sprintf('%d-ply skin, %.2f mm', n, 1e3*n*M.twill.t), 'Color', ink, ...
        'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'BackgroundColor', col);
end
for j = 1:size(lofts, 1)
    patch('XData', 1e3*[lofts(j,1) lofts(j,2) lofts(j,2) lofts(j,1)], 'YData', [0 0 chord chord], ...
        'FaceColor', tint_loft, 'EdgeColor', 'none');
    text(1e3*mean(lofts(j,:)), 0.62*chord, 'loft', 'Color', ink, 'FontSize', fs, 'HorizontalAlignment', 'center');
end
plot([0 1e3*s_tip], 1e3*D.spar_aft*[1 1], 'Color', orange, 'LineWidth', 3);
text(1e3*mean(ribs(1:2)), 1e3*D.spar_aft + 14, 'spar plate', 'Color', orange, 'FontSize', fs, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
for i = 1:numel(ribs)
    r = ribs(i); col = grey; lw = 2;
    if any(abs(r - others) < 1e-6), col = blue; lw = 4; end
    if abs(r - cfg.mount_s) < 1e-6, col = orange; lw = 4; end
    if abs(r - s_tip) < 1e-6, col = ink; lw = 4; end
    plot(1e3*r*[1 1], [0 chord], 'Color', col, 'LineWidth', lw);
    text(1e3*r, -10, sprintf('%.0f', 1e3*r), 'Color', col, 'FontSize', fs - 1, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
text(1e3*cfg.mount_s, chord + 12, 'mount rib', 'Color', orange, 'FontSize', fs, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
if ~isempty(others)
    lab = 'corner rib'; if ~isempty(lofts), lab = 'loft-end ribs (blue)'; end
    text(1e3*others(1) + 4, chord + 12, lab, 'Color', blue, 'FontSize', fs, 'VerticalAlignment', 'top');
end
text(1e3*s_tip, chord + 12, 'endplate', 'Color', ink, 'FontSize', fs, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
xlim([-5 1e3*s_tip + 5]); ylim([-55 chord + 95]);
set(gca, 'YTick', [0 chord], 'YTickLabel', {'LE', 'TE'}, 'XTick', []);
title(sprintf('Half wing, top view: %d ribs per side (%d in the wing), stations in mm', n_half, n_full));

subplot(2, 1, 2); hold on; grid on; box off; set(gca, 'FontSize', fs);
h1 = plot(1e3*v.s, 1e3*min(v.a, 0.5), 'Color', blue, 'LineWidth', 2.5);
h2 = [];
for i = 1:size(bays, 1)
    h2 = plot(1e3*[max(bays(i,1), 0) bays(i,2)], 1e3*bays(i,3)*[1 1], 'Color', orange, 'LineWidth', 5);
end
xlim([-5 1e3*s_tip + 5]); ylim([0 200]);
xlabel('Span from centre plane [mm]'); ylabel('Rib spacing [mm]');
legend([h1 h2], {sprintf('Allowable spacing (skin buckling, %s, K = %g)', cfg.rib_check_load, cfg.rib_K), ...
    'Actual bay length'}, 'Location', 'northwest');
print(f3, '-dpng', '-r200', fullfile(out_dir, 'fig3_rib_layout.png'));

%% Summary
fprintf('Variant: %s  (K = %g, checked at %s)\n', v.name, cfg.rib_K, cfg.rib_check_load);
fprintf('Limit load per half wing: main element %.0f N, flaps + endplate %.0f N\n', ...
    -sum(Fy)*Ld.k_limit, -Ld.F_ep(2)*Ld.k_limit);
fprintf('Ultimate bending moment at the mount: %.1f N m\n', M_peak);
fprintf('Rib stations [mm]: %s\n', strjoin(arrayfun(@(x) sprintf('%.0f', 1e3*x), ribs', 'UniformOutput', false), ', '));
fprintf('Bay / allowable [mm]: %s\n', strjoin(arrayfun(@(i) sprintf('%.0f/%.0f', 1e3*bays(i,3), 1e3*bays(i,4)), ...
    1:size(bays,1), 'UniformOutput', false), ', '));
fprintf('Tightest bay uses %.0f%% of its allowable spacing\n', 100*max(bays(:,3) ./ bays(:,4)));
fprintf('Skin + plate %.0f g per half, tip deflection %.2f mm, twist %.2f deg (limit), min strength margin %.1f\n', ...
    1e3*v.m, v.w, v.th, v.MS);
fprintf('Figures written to %s\n', out_dir);
