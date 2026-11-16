function T = transform_tris(T, R, t)
% TRANSFORM_TRIS  Apply x_new = R*x + t to every vertex of an nT x 9
% triangle array (already in metres).
t = t(:)';
for k = 0:2
    c = 3*k + (1:3);
    T(:, c) = T(:, c) * R' + t;
end
end
