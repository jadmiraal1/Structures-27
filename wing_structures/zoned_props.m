function [P, lams] = zoned_props(outlines, s_st, ply, base_layup, zones, pl, x_plate, gap)
% ZONED_PROPS  station_props with a skin layup that can change along the
% span (e.g. extra plies inboard).
%   ply        : ply properties (bound_props(M.twill, k))
%   base_layup : layup used where no zone applies, e.g. [0 45 0]
%   zones      : struct array with fields s ([start end], m) and layup;
%                later zones override earlier ones; [] for none
% Returns P (as station_props) plus per-station skin data used for
% buckling: P.D11, P.Ex_skin, P.t_skin, and the laminate for each station.
ns = numel(s_st);
lams = cell(ns, 1);
for k = 1:ns
    lay = base_layup;
    for z = 1:numel(zones)
        if s_st(k) >= zones(z).s(1) && s_st(k) <= zones(z).s(2), lay = zones(z).layup; end
    end
    lams{k} = laminate_clt(ply, lay);
end
for k = 1:ns
    Pk = station_props(outlines(k), s_st(k), lams{k}, pl, x_plate, gap);
    if k == 1, P = Pk; f = fieldnames(Pk); continue; end
    for i = 1:numel(f)
        v = Pk.(f{i});
        if isnumeric(v) && numel(v) == 1 && ~strcmp(f{i}, 'x_plate'), P.(f{i})(k,1) = v; end
    end
end
P.D11 = cellfun(@(l) l.D(1,1), lams);
P.Ex_skin = cellfun(@(l) l.Ex, lams);
P.t_skin = cellfun(@(l) l.h, lams);
P.areal = cellfun(@(l) l.areal, lams);
end
