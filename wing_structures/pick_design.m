% PICK_DESIGN  Choose the lightest design from the main_sizing results that
% meets strength (margins >= 0 at ultimate) and your stiffness limits at
% limit load (cfg.max_tip_deflection_mm, cfg.max_tip_twist_deg).
% If the limits are not set, lists the lightest strength-passing design
% for each skin layup, so you can see the trade-off and set the limits.
% Needs results/<name>_sizing.mat from main_sizing. Run from wing_structures.
clear;
cfg = front_wing_config();
S = load(fullfile('results', [cfg.name '_sizing.mat'])); R = S.R;
R = R([R.pass]);                                   % strength at ultimate
m = [R.m_tot]'; w = [R.w_tip_mm]'; th = abs([R.twist_tip_deg])';
hdr = sprintf('%-14s %8s %8s %9s %10s %10s %8s %8s', 'skin', 'spar aft', 'plate', 'mass [g]', ...
              'tip w[mm]', 'twist[deg]', 'MS min', 'rib max');
row = @(r) sprintf('%-14s %6.0fmm %6.3fin %9.0f %10.2f %10.3f %8.2f %6.0fmm', r.layup, r.spar_aft_mm, ...
               r.t_in, 1e3*r.m_tot, r.w_tip_mm, abs(r.twist_tip_deg), min(r.MS_plate, r.MS_skin), r.rib_max_mm);

limw = cfg.max_tip_deflection_mm; limt = cfg.max_tip_twist_deg;
if ~isnan(limw) || ~isnan(limt)
    ok = true(size(m));
    if ~isnan(limw), ok = ok & w <= limw; end
    if ~isnan(limt), ok = ok & th <= limt; end
    fprintf('Limits: tip deflection <= %g mm, tip twist <= %g deg (limit load)\n', limw, limt);
    if ~any(ok)
        fprintf('No design in the sweep meets these limits. Add stiffer skins to cfg.skin_options.\n');
        return;
    end
    idx = find(ok); [~, o] = sort(m(idx)); idx = idx(o);
    fprintf('\nLightest design meeting strength and stiffness:\n%s\n%s\n', hdr, row(R(idx(1))));
    fprintf('\nNext lightest:\n');
    for k = 2:min(6, numel(idx)), fprintf('%s\n', row(R(idx(k)))); end
else
    fprintf('Stiffness limits not set (cfg.max_tip_deflection_mm / cfg.max_tip_twist_deg).\n');
    fprintf('Lightest strength-passing design for each skin:\n\n%s\n', hdr);
    lays = unique({R.layup}, 'stable');
    for k = 1:numel(lays)
        idx = find(strcmp({R.layup}, lays{k}));
        [~, b] = min(m(idx)); fprintf('%s\n', row(R(idx(b))));
    end
end
fprintf(['\nNotes: spar positions within a few grams are equivalent - choose the position for the\n' ...
         'mount hardware and rib layout. The plate gauge here is the strength minimum only; the\n' ...
         'mount (bolt bearing, rib attachment) may need it thicker.\n']);
