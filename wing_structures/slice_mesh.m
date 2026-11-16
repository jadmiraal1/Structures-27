function loops = slice_mesh(T, axis_n, z0)
% SLICE_MESH  Cut a triangle mesh with the plane {x : x(axis_n) = z0}.
%   T      : nT x 9 triangle array
%   axis_n : 1, 2 or 3 - coordinate axis normal to the cutting plane
%   z0     : plane position
% Returns a cell array of closed loops, each an M x 3 array of points in
% order (last point not repeated), sorted by enclosed length, longest first.

V = reshape(T', 3, 3, [])  ;     % V(:, vertex, tri)
zz = squeeze(V(axis_n, :, :))';  % nT x 3 plane-normal coordinate
span = max(abs(zz(:))) + 1;
% nudge the plane off any vertex lying exactly on it
eps_z = 1e-9 * span;
while any(abs(zz(:) - z0) < eps_z)
    z0 = z0 + 3.1*eps_z;
end
s = zz - z0;
cross_t = find(min(s, [], 2) < 0 & max(s, [], 2) > 0);
nseg = numel(cross_t);
if nseg == 0, loops = {}; return; end

seg = zeros(nseg, 6);
edges = [1 2; 2 3; 3 1];
for q = 1:nseg
    i = cross_t(q);
    pts = zeros(2, 3); k = 0;
    for e = 1:3
        sa = s(i, edges(e,1)); sb = s(i, edges(e,2));
        if (sa < 0) ~= (sb < 0)
            pa = V(:, edges(e,1), i)'; pb = V(:, edges(e,2), i)';
            k = k + 1;
            pts(k,:) = pa + (sa/(sa - sb)) * (pb - pa);
        end
    end
    seg(q,:) = [pts(1,:) pts(2,:)];
end

% Merge coincident endpoints (shared edges give identical points up to
% round-off) and chain segments into loops.
P = [seg(:,1:3); seg(:,4:6)];
tolm = 1e-9 * span;
[~, ~, id] = unique(round(P / tolm), 'rows');
e1 = id(1:nseg); e2 = id(nseg+1:end);
npt = max(id);
pt = zeros(npt, 3); pt(id, :) = P;
adj = cell(npt, 1);
for q = 1:nseg
    adj{e1(q)}(end+1) = q; adj{e2(q)}(end+1) = q;
end
used = false(nseg, 1);
loops = {};
for q0 = 1:nseg
    if used(q0), continue; end
    used(q0) = true;
    chain = [e1(q0) e2(q0)];
    cur = e2(q0);
    while true
        nxt = adj{cur}(~used(adj{cur}));
        if isempty(nxt), break; end
        q = nxt(1); used(q) = true;
        if e1(q) == cur, cur = e2(q); else, cur = e1(q); end
        if cur == chain(1), break; end
        chain(end+1) = cur; %#ok<AGROW>
    end
    loops{end+1} = pt(chain, :); %#ok<AGROW>
end
len = cellfun(@(L) sum(sqrt(sum(diff([L; L(1,:)]).^2, 2))), loops);
[~, o] = sort(len, 'descend');
loops = loops(o);
end
