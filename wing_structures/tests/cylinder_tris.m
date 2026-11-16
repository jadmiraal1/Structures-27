function T = cylinder_tris(r, h, n)
% CYLINDER_TRIS  Open cylinder along z (side wall only), n facets.
th = 2*pi*(0:n-1)'/n; th2 = [th(2:end); th(1)];
a = [r*cos(th) r*sin(th) zeros(n,1)];  b = [r*cos(th2) r*sin(th2) zeros(n,1)];
c = [b(:,1:2) h*ones(n,1)];            d = [a(:,1:2) h*ones(n,1)];
T = [a b c; a c d];
end
