function [ribs, bays] = rib_layout(s, a_max, fixed, s_tip)
% RIB_LAYOUT  Place ribs so that no bay is longer than the allowable rib
% spacing anywhere inside it.
%   s, a_max : spanwise grid [m] and allowable rib spacing there [m]
%   fixed    : ribs that must exist (mount, loft corners) [m]
%   s_tip    : tip position [m] (the endplate closes the tip; listed as a rib)
%
% Between consecutive fixed ribs: first the smallest number of bays is
% found by marching outward and making each bay as long as allowed
% (length <= min(a_max over the bay)). The ribs are then spread so every
% bay uses the same fraction of its allowable length, which keeps the
% count minimal while following the local allowable spacing.
% The centre plane is a plane of symmetry, so the bays between the first
% rib on each side (-r1 ... +r1) are spaced evenly and symmetrically.
%   ribs : rib positions on one half wing [m] (0 = a rib on the centre plane)
%   bays : rows [start end length a_allow] for the half wing; a bay that
%          crosses the centre plane is listed as [-r r 2r a_allow]
s = s(:); a_max = a_max(:);
fixed = unique([fixed(:); s_tip]);
% allowable spacing for a bay between x0 and x1 (uses |s|, symmetric)
amin = @(x0, x1) min([a_max(s >= bmin(x0, x1) - 1e-12 & s <= max(abs([x0 x1])) + 1e-12); ...
                      interp1(s, a_max, [bmin(x0, x1); max(abs([x0 x1]))], 'linear', 'extrap')]);

% centre: symmetric layout on [-r1, r1]
r1 = fixed(1);
for n = 1:500
    p = linspace(-r1, r1, n + 1);
    if all(arrayfun(@(q) p(q+1) - p(q) <= amin(p(q), p(q+1)) + 1e-12, 1:n)), break; end
end
ribs = p(p >= -1e-12)';
% outer bays between fixed ribs
for k = 1:numel(fixed) - 1
    x0 = fixed(k); x1 = fixed(k+1);
    n = numel(march(x0, x1, 1, amin));
    f_lo = 0; f_hi = 1;                    % smallest scale on a_max that still needs only n bays
    for it = 1:40
        f = (f_lo + f_hi)/2;
        if numel(march(x0, x1, f, amin)) <= n, f_hi = f; else, f_lo = f; end
    end
    p = march(x0, x1, f_hi, amin);
    ribs = [ribs; p(:)]; %#ok<AGROW>
end
ribs = unique(round(ribs * 1e7) / 1e7);
ribs(abs(ribs) < 1e-9) = 0;

full = unique([-ribs(ribs > 0); ribs]);
bays = [];
for k = 1:numel(full) - 1
    x0 = full(k); x1 = full(k+1);
    if x1 <= 0, continue; end
    bays = [bays; x0 x1 x1 - x0 amin(x0, x1)]; %#ok<AGROW>
end
end

function m = bmin(x0, x1)
% smallest |s| covered by the bay
if x0 < 0 && x1 > 0, m = 0; else, m = min(abs([x0 x1])); end
end

function p = march(x0, x1, f, amin)
% rib positions after x0 up to and including x1, each bay as long as
% f * (allowable spacing over the bay)
p = []; x = x0;
while x < x1 - 1e-9
    lo = x; hi = x1;
    if hi - x <= f*amin(x, hi) + 1e-12
        x = hi;
    else
        for it = 1:40
            m = (lo + hi)/2;
            if m - x <= f*amin(x, m) + 1e-12, lo = m; else, hi = m; end
        end
        if lo - x < 1e-6, error('rib_layout:spacing', 'Allowable spacing too small near s = %.3f m', x); end
        x = lo;
    end
    p(end+1) = x; %#ok<AGROW>
end
end
