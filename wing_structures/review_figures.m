% REVIEW_FIGURES  Design-review figures for one rib variant.
%
% Reads the saved results of main_geometry, main_loads and main_ribs and
% writes three PNG figures to results/review:
%   fig1_loads.png       spanwise downforce at limit load, bending moment at ultimate
%   fig2_section.png     section at the mount: skin, spar plate, neutral axis
%   fig3_rib_layout.png  rib stations on the planform, bay length vs allowable spacing
%
% Run main_ribs first (and again after any change to the wing's config),
% then run this from wing_structures. The wing is the active one (wing_config). No sizing is done here: every number
% comes from the results files and the config.
clear; close all;
variant = '';                            % any text contained in the variant name (cfg.rib_variants);
                                         % '' = the wing's cfg.review_variant

cfg = wing_config();                  % active wing: see wing_config.m
M   = materials_input();
Gm  = load(fullfile('results', [cfg.name '_geometry.mat']));
Ld  = load(fullfile('results', [cfg.name '_loads.mat'])); Ld = Ld.loads;
Rb  = load(fullfile('results', [cfg.name '_ribs.mat'])); V = Rb.V;

% Refuse rib results that were not made with the current code and config
no_loft = isfield(cfg, 'rib_no_loft') && cfg.rib_no_loft;
if ~isfield(Rb, 'meta')
    error('review_figures:stale', ['%s_ribs.mat was written by an older main_ribs, so it does not match ' ...
        '%s_config.m.\nRun main_ribs, then review_figures.'], cfg.name, cfg.name);
end
meta = Rb.meta;
if meta.rib_K ~= cfg.rib_K || ~strcmpi(meta.rib_check_load, cfg.rib_check_load) ...
        || logical(meta.rib_no_loft) ~= no_loft || abs(meta.mount_s - cfg.mount_s) > 1e-9 ...
        || ~strcmpi(meta.mount_bc, cfg.mount_bc) ...
        || (isfield(meta, 'strain_model') && isfield(cfg, 'strain_model') && ~strcmpi(meta.strain_model, cfg.strain_model))
    error('review_figures:stale', ['%s_ribs.mat does not match its config (it was saved with K = %g, ' ...
        'check at %s, no ribs in lofts = %d, mount %s).\nRun main_ribs, then review_figures.'], ...
        cfg.name, meta.rib_K, meta.rib_check_load, meta.rib_no_loft, meta.mount_bc);
end
if isempty(variant), variant = cfg.review_variant; end
out_dir = fullfile('results', ['review_' cfg.name]);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

iv = find(arrayfun(@(x) ~isempty(strfind(lower(x.name), lower(variant))), V), 1);
if isempty(iv)
    error('review_figures:variant', 'No variant in %s_ribs.mat contains "%s".', cfg.name, variant);
end
v  = V(iv);
ic = find(cellfun(@(c) strcmp(c.name, v.name), cfg.rib_variants), 1);
if isempty(ic)
    error('review_figures:stale', 'Variant "%s" is not in %s_config any more. Re-run main_ribs.', v.name, cfg.name);
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
bad = [];                                             % ribs that fell inside a loft (skin too thin to span it)
if isfield(v, 'ribs_in_loft'), bad = v.ribs_in_loft(:); end
if ~isempty(bad)
    warning('review_figures:loft', 'Variant "%s" has ribs inside a loft (%s mm): its skin cannot span the loft. Pick another variant.', ...
        v.name, strjoin(arrayfun(@(x) sprintf('%.0f', 1e3*x), bad', 'UniformOutput', false), ', '));
end

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
% Bending moment at ultimate. Fallback from statics (pinned mount, symmetric
% centre plane, F_z * height term left out); replaced below by the beam
% model's M_hog when main_ribs saved it.
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
if isfield(v, 'M_ult') && ~isempty(v.M_ult)          % the beam model's own moment, when main_ribs saved it
    sg = v.s(:); Mh = v.M_ult(:);
end
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
% Neutral axis under a vertical load. 'horizontal': through the skin centroid.
% 'unsymmetric': tilted by Ixz/Izz (uniform skin, so the modulus cancels), the
% same assumption beam_skin_strains makes. ea = vertical distance from it [mm].
strain_model = 'horizontal';
if isfield(cfg, 'strain_model'), strain_model = cfg.strain_model; end
xc = sum(Lg .* mid(:,1)) / sum(Lg);
dsg = p2 - xy; dxm = mid(:,1) - xc; dym = mid(:,2) - yc;
Izz = sum(Lg .* (dxm.^2 + dsg(:,1).^2/12)); Ixz = sum(Lg .* (dxm.*dym + dsg(:,1).*dsg(:,2)/12));
na_slope = 0;
if strcmpi(strain_model, 'unsymmetric'), na_slope = Ixz / Izz; end
ea = ya - 1e3 * na_slope * (xy(:,1) - xc);
na_y = @(xq) -na_slope * (xq - 1e3*(max(xy(:,1)) - xc));   % neutral axis height at xq mm behind the LE
[e_bot, ib] = min(ea);                                % most compressed fibre when the wing hogs
[e_top, it] = max(ea);

f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 520]);
axes('Position', [0.08 0.42 0.88 0.48]); hold on; grid on; box off; set(gca, 'FontSize', fs);
plot([xa; xa(1)], [ya; ya(1)], 'Color', ink, 'LineWidth', 2.5);
plot([-8 chord + 8], na_y([-8 chord + 8]), '--', 'Color', grey, 'LineWidth', 1.2);
plot([xs xs], [y_lo + clr, y_hi - clr], 'Color', orange, 'LineWidth', 4);
plot(xa(ib), ya(ib), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', blue, 'MarkerFaceColor', blue);
plot(xa(it), ya(it), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', blue, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
xlim([-8 chord + 8]); ylim([min(ya) - 10, max(ya) + 10]);
daspect([1 1 1]);
xlabel('Distance behind leading edge [mm]'); ylabel('Height [mm]');
title(sprintf('Section at the mount (%.0f mm station): chord %.0f mm', 1e3*Gm.G.s(k), chord));
% key, drawn by hand in its own axes so it never covers the section
axes('Position', [0.08 0.02 0.88 0.27]); hold on; axis off; xlim([0 1]); ylim([0 1]);
key = {'Skin (outer mould line)', 'Neutral axis under a vertical load', ...
    sprintf('Spar plate: %.0f mm behind LE, %.0f mm tall', 1e3*D.spar_aft, h_plate), ...
    sprintf('Most compressed skin under downforce: %.1f mm from the neutral axis', -e_bot), ...
    sprintf('Most stretched skin under downforce: %.1f mm from the neutral axis', e_top)};
ky = 0.92 - 0.2*(0:4);
plot([0.01 0.05], ky(1)*[1 1], 'Color', ink, 'LineWidth', 2.5);
plot([0.01 0.05], ky(2)*[1 1], '--', 'Color', grey, 'LineWidth', 1.2);
plot([0.01 0.05], ky(3)*[1 1], 'Color', orange, 'LineWidth', 4);
plot(0.03, ky(4), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', blue, 'MarkerFaceColor', blue);
plot(0.03, ky(5), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', blue, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
for r = 1:numel(key)
    text(0.07, ky(r), key{r}, 'Color', ink, 'FontSize', fs, 'VerticalAlignment', 'middle');
end
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
text(1e3*ribs(1) + 4, 1e3*D.spar_aft + 14, 'spar plate', 'Color', orange, 'FontSize', fs, 'VerticalAlignment', 'top');
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
xlim([-5 1e3*s_tip + 5]); ylim([-55 chord + 110]);
set(gca, 'YTick', [0 chord], 'YTickLabel', {'LE', 'TE'}, 'XTick', []);
ttl = sprintf('Half wing, top view: %d ribs per side (%d in the wing), stations in mm', n_half, n_full);
if ~isempty(bad), ttl = [ttl ' - NOT VALID: ribs inside a loft']; end
title(ttl);

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
