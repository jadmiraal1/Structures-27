% RUN_TESTS_GEOMETRY  Checks for the STL, distance, slicing and load
% integration tools against exact results.
% Run from the wing_structures folder:  addpath('tests'); run_tests_geometry
nfail = 0;

%% 1. STL round trip: binary and ASCII
T0 = box_tris([0 0 0], [1 2 3]);
fb = [tempname '.stl']; write_stl_binary(fb, T0);
Tb = read_stl(fb, 1);
nfail = nfail + check('STL binary round trip', max(abs(Tb(:) - T0(:))) + 1, 1, 1e-6);
fa = [tempname '.stl']; fid = fopen(fa, 'w');
fprintf(fid, 'solid test\n');
for i = 1:size(T0,1)
    fprintf(fid, ' facet normal 0 0 0\n  outer loop\n');
    fprintf(fid, '   vertex %.9g %.9g %.9g\n', reshape(T0(i,:), 3, 3));
    fprintf(fid, '  endloop\n endfacet\n');
end
fprintf(fid, 'endsolid test\n'); fclose(fid);
Ta = read_stl(fa, 1e-3);
nfail = nfail + check('STL ASCII round trip + scale', max(abs(Ta(:) - 1e-3*T0(:))) + 1, 1, 1e-12);
delete(fb); delete(fa);

%% 2. Point-triangle distance, exact cases
a = [0 0 0]; b = [1 0 0]; c = [0 1 0];
P = [0.2 0.2 0.5;    % above interior -> 0.5
     -1 -1 0;        % beyond vertex a -> sqrt(2)
     0.5 -2 0;       % beyond edge ab -> 2
     1 1 0;          % beyond edge bc -> sqrt(2)/2
     2 0 1];         % beyond vertex b -> sqrt(2)
dexp = [0.5; sqrt(2); 2; sqrt(2)/2; sqrt(2)];
d = point_triangle_distance(P, a, b, c);
nfail = nfail + check('Point-triangle exact cases', max(abs(d - dexp)) + 1, 1, 1e-12);

%% 3. Point-triangle distance vs dense sampling (random points)
rand('seed', 1); %#ok<RAND>
a = [0.1 0.3 -0.2]; b = [1.2 0.1 0.4]; c = [0.4 1.1 0.9];
P = 3*rand(400, 3) - 1;
d = point_triangle_distance(P, a, b, c);
[u, v] = meshgrid(linspace(0, 1, 301));
keep = u + v <= 1; u = u(keep); v = v(keep);
S = a + u.*(b - a) + v.*(c - a);
ds = zeros(size(P,1), 1);
for i = 1:size(P,1), ds(i) = sqrt(min(sum((S - P(i,:)).^2, 2))); end
h = norm(b - a)/300 + norm(c - a)/300;
nfail = nfail + check('Point-triangle never above sampled min', double(all(d <= ds + 1e-12)), 1, 0);
nfail = nfail + check('Point-triangle within sampling spacing', double(all(ds - d <= h)), 1, 0);

%% 4. Distance to mesh with tolerance cut-off
T = box_tris([0 0 0], [0.1 0.05 0.2]);
P = [0.05 0.025 0.201;    % 1 mm above top face
     0.05 0.0505 0.1;     % 0.5 mm off a side
     0.05 0.025 0.21];    % 10 mm away -> Inf at tol 2 mm
d = distance_to_mesh(P, T, 2e-3);
nfail = nfail + check('Mesh distance (1 mm)',   d(1), 1e-3, 1e-9);
nfail = nfail + check('Mesh distance (0.5 mm)', d(2), 5e-4, 1e-9);
nfail = nfail + check('Mesh distance beyond tol is Inf', double(isinf(d(3))), 1, 0);

%% 4b. Distance to mesh with many small triangles vs brute force
Tc = cylinder_tris(0.05, 0.003, 400);                 % small facets -> grid path
Tc = [Tc; transform_tris(Tc, eye(3), [0 0 0.003])];   % stack two rings
rand('seed', 3); %#ok<RAND>
th = 2*pi*rand(300,1); rr = 0.05 + 0.006*(rand(300,1) - 0.5);
Pc = [rr.*cos(th) rr.*sin(th) 0.006*rand(300,1)];
dg = distance_to_mesh(Pc, Tc, 2e-3);
db = inf(300,1);
for i = 1:size(Tc,1)
    db = min(db, point_triangle_distance(Pc, Tc(i,1:3), Tc(i,4:6), Tc(i,7:9)));
end
db(db > 2e-3) = inf;
same = (isinf(dg) & isinf(db)) | abs(dg - db) < 1e-12;
nfail = nfail + check('Mesh distance (grid) = brute force', double(all(same)), 1, 0);

%% 5. Slice a faceted cylinder: perimeter and area of an n-gon
r = 0.05; n = 64;
L = slice_mesh(cylinder_tris(r, 0.2, n), 3, 0.1);
xy = L{1}(:,1:2);
perim = sum(sqrt(sum(diff([xy; xy(1,:)]).^2, 2)));
[~, ~, Ap] = polygon_centroid(xy(:,1), xy(:,2));
nfail = nfail + check('Slice: one loop', numel(L), 1, 0);
nfail = nfail + check('Slice: n-gon perimeter', perim, 2*n*r*sin(pi/n), 1e-9);
nfail = nfail + check('Slice: n-gon area', Ap, 0.5*n*r^2*sin(2*pi/n), 1e-9);

%% 6. Integration: closed surface, p = a*y  ->  sum(p A) = a * V * y_hat
lo = [0.1 0.2 -0.3]; hi = [0.4 0.35 -0.05];
T = box_tris(lo, hi);
Avec = 0.5 * cross(T(:,4:6) - T(:,1:3), T(:,7:9) - T(:,1:3), 2);  % outward
xc = (T(:,1:3) + T(:,4:6) + T(:,7:9)) / 3;
Fs.p = 1000 * xc(:,2); Fs.A = Avec; Fs.Amag = sqrt(sum(Avec.^2, 2)); Fs.tau = zeros(size(Avec));
Ff = face_forces(Fs);
V = prod(hi - lo);
nfail = nfail + check('Hydrostatic box: Fy = a V', sum(Ff(:,2)), 1000*V, 1e-12);
nfail = nfail + check('Hydrostatic box: Fx ~ 0', sum(Ff(:,1)) + 1, 1, 1e-12);
Fs.p = 500 + 0*xc(:,2);
Ff0 = face_forces(Fs);
nfail = nfail + check('Uniform pressure: net 0', norm(sum(Ff0, 1)) + 1, 1, 1e-12);

%% 7. Strip integration conserves totals and moments
rand('seed', 2); %#ok<RAND>
N = 5000; xyz = [rand(N,1) rand(N,1) -0.6*rand(N,1)];
Fq = randn(N, 3); Aq = rand(N, 1);
s = -xyz(:,3); edges = (0:0.05:0.6)';
Ls = integrate_strips(xyz, Fq, Aq, s, edges, [0.3 0.1], @(s) -s);
nfail = nfail + check('Strips: force total', norm(sum(Ls.F,1) - sum(Fq,1)) + 1, 1, 1e-10);
% moment of each strip about its centre, moved to the origin, sums to total
Mtot = sum(cross(xyz, Fq, 2), 1);
Msum = zeros(1,3);
for k = 1:numel(Ls.s_mid)
    r0 = [0.3 0.1 -Ls.s_mid(k)];
    Msum = Msum + Ls.M(k,:) + cross(r0, Ls.F(k,:));
end
nfail = nfail + check('Strips: moment transfer', norm(Msum - Mtot) + 1, 1, 1e-10);

%% 8. Section depth on a NACA 0012 loop (x forward, LE at max x)
xc8 = (1 - cos(linspace(0, pi, 200)'))/2;
yt = 0.6*(0.2969*sqrt(xc8) - 0.1260*xc8 - 0.3516*xc8.^2 + 0.2843*xc8.^3 - 0.1036*xc8.^4);
c = 0.3;
up = [c*(1 - xc8)  c*yt]; lo8 = [c*(1 - xc8) -c*yt];
loop = [up; flipud(lo8(2:end-1,:))];
S = section_depth(loop, 1001);
nfail = nfail + check('Section depth: chord', S.chord, c, 1e-12);
nfail = nfail + check('Section depth: t_max/c', S.t_max/c, 0.12, 2e-3);
nfail = nfail + check('Section depth: x/c at t_max', S.xc_tmax, 0.30, 0.03);

%% 9. Gap filling: good strips untouched, bad strips linear in s
Lf.s_mid = (0.01:0.02:0.19)'; Lf.ds = 0.02*ones(10,1);
Lf.f = [zeros(10,1) 100 - 50*Lf.s_mid zeros(10,1)];   % linear in s
Lf.m = [zeros(10,2) 3*Lf.s_mid];
Lf.F = Lf.f .* Lf.ds; Lf.M = Lf.m .* Lf.ds;
cov = ones(10,1); cov([4 5 6]) = [0.4 0.05 0.6];
Lt = Lf; Lt.f(4:6,:) = 0; Lt.m(4:6,:) = 0;          % corrupted strips
[Lg, filled] = fill_low_coverage(Lt, cov, 0.9);
nfail = nfail + check('Fill: flags bad strips', double(isequal(find(filled)', [4 5 6])), 1, 0);
nfail = nfail + check('Fill: recovers linear f', max(abs(Lg.f(:) - Lf.f(:))) + 1, 1, 1e-12);
nfail = nfail + check('Fill: recovers linear m', max(abs(Lg.m(:) - Lf.m(:))) + 1, 1, 1e-12);
nfail = nfail + check('Fill: F = f*ds rebuilt', max(abs(Lg.F(:) - Lf.F(:))) + 1, 1, 1e-12);

fprintf('\n%d test(s) failed.\n', nfail);
