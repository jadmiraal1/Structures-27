function [d, itri] = distance_to_mesh(P, T, tol)
% DISTANCE_TO_MESH  Exact distance from points P (N x 3) to a triangle
% mesh T (nT x 9), for points within TOL of the mesh. Points farther than
% TOL get d = Inf. itri is the index of the nearest triangle (0 if none).
%
% Points and triangles are bucketed on a uniform grid (cell size h). Each
% small triangle is registered in every cell its bounding box (grown by
% TOL) overlaps; then for each cell all its points are tested against all
% its triangles in one vectorised call. Triangles spanning more than 3
% cells in any direction (e.g. long facets on a straight extrusion) are
% handled one at a time against the points near them.

N = size(P, 1);
d = inf(N, 1); itri = zeros(N, 1);
h  = max(tol, 4e-3);
lo = min(P, [], 1) - h;
ijk  = floor((P - lo) / h) + 1;
dims = max(ijk, [], 1) + 1;
pkey = sub2ind(dims, ijk(:,1), ijk(:,2), ijk(:,3));
[pkey_s, porder] = sort(pkey);
ncell  = prod(dims);
pcnt   = accumarray(pkey_s, 1, [ncell 1]);
pstart = [1; cumsum(pcnt) + 1];

% triangle bounding boxes in cell indices
X = T(:, [1 4 7]); Y = T(:, [2 5 8]); Z = T(:, [3 6 9]);
bmin = [min(X,[],2) min(Y,[],2) min(Z,[],2)] - tol;
bmax = [max(X,[],2) max(Y,[],2) max(Z,[],2)] + tol;
i0 = max(floor((bmin - lo)/h) + 1, 1);
i1 = min(floor((bmax - lo)/h) + 1, dims);
valid = all(i1 >= i0, 2);
span  = i1 - i0 + 1;
small = valid & all(span <= 3, 2);
big   = find(valid & ~small);

%% Small triangles: (cell, triangle) pairs, then per-cell vectorised test
ts = find(small);
ck = []; tk = [];
for ox = 0:2
    for oy = 0:2
        for oz = 0:2
            m = span(ts,1) > ox & span(ts,2) > oy & span(ts,3) > oz;
            t = ts(m);
            c = sub2ind(dims, i0(t,1) + ox, i0(t,2) + oy, i0(t,3) + oz);
            keep = pcnt(c) > 0;               % only cells that hold points
            ck = [ck; c(keep)]; tk = [tk; t(keep)]; %#ok<AGROW>
        end
    end
end
[ck, o] = sort(ck); tk = tk(o);
[cells, first] = unique(ck, 'first');
last = [first(2:end) - 1; numel(ck)];
for q = 1:numel(cells)
    tri = tk(first(q):last(q));
    pid = porder(pstart(cells(q)):pstart(cells(q)+1) - 1);
    np = numel(pid); nt = numel(tri);
    ip = repmat((1:np)', nt, 1);
    it = reshape(repmat(tri', np, 1), [], 1);
    dd = point_triangle_distance(P(pid(ip),:), T(it,1:3), T(it,4:6), T(it,7:9));
    dd = reshape(dd, np, nt);
    [dmin, j] = min(dd, [], 2);
    better = dmin < d(pid) & dmin <= tol;
    d(pid(better)) = dmin(better);
    itri(pid(better)) = tri(j(better));
end

%% Large triangles: one at a time
for i = big'
    [I, J, K] = ndgrid(i0(i,1):i1(i,1), i0(i,2):i1(i,2), i0(i,3):i1(i,3));
    c = sub2ind(dims, I(:), J(:), K(:));
    c = c(pcnt(c) > 0);
    if isempty(c), continue; end
    idx = cell2mat(arrayfun(@(q) (pstart(q):pstart(q+1)-1)', c, 'UniformOutput', false));
    pid = porder(idx);
    di = point_triangle_distance(P(pid,:), T(i,1:3), T(i,4:6), T(i,7:9));
    better = di < d(pid) & di <= tol;
    d(pid(better)) = di(better);
    itri(pid(better)) = i;
end
end
