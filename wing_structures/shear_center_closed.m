function [x_sc, info] = shear_center_closed(xy, t)
% SHEAR_CENTER_CLOSED  Chordwise shear-centre position of a single-cell
% closed thin-walled section (uniform material), for vertical shear.
%   xy : M x 2 wall midline points [x horizontal, y vertical], in order
%   t  : wall thickness, scalar or M x 1 (segment i -> i+1)
% Method (Megson, Aircraft Structures, closed-section shear): open-section
% shear flow q_b from a cut, closing flow q_s0 from the zero-twist
% condition, then moment equilibrium about the centroid gives x_sc.
% A vertical load through x_sc bends the section without twisting it.
M  = size(xy, 1);
if isscalar(t), t = t*ones(M,1); else, t = t(:); end
p1 = xy; p2 = xy([2:M 1], :);
dx = p2(:,1) - p1(:,1); dy = p2(:,2) - p1(:,2);
L  = sqrt(dx.^2 + dy.^2);
xm = (p1(:,1) + p2(:,1))/2; ym = (p1(:,2) + p2(:,2))/2;
w  = t .* L;
xc = sum(w.*xm)/sum(w); yc = sum(w.*ym)/sum(w);
X  = xm - xc; Y = ym - yc;
Ixx = sum(w .* (Y.^2 + dy.^2/12));
Iyy = sum(w .* (X.^2 + dx.^2/12));
Ixy = sum(w .* (X.*Y + dx.*dy/12));
D   = Ixx*Iyy - Ixy^2;

Sy = 1;
dq = (Sy*Ixy/D) * (t.*L.*X) - (Sy*Iyy/D) * (t.*L.*Y);   % increment over each segment
qn = [0; cumsum(dq)];                                   % at nodes, cut at node 1
qb = (qn(1:M) + qn(2:M+1)) / 2;                          % segment average
qs0 = -sum(qb .* L ./ t) / sum(L ./ t);                  % zero twist
q = qb + qs0;

Fy = sum(q .* dy);                    % resultant (should be +-Sy)
Fx = sum(q .* dx);
Mz = sum((X .* dy - Y .* dx) .* q);   % about the centroid
x_sc = xc + Mz / Fy;
info = struct('xc', xc, 'yc', yc, 'Fy', Fy, 'Fx', Fx, 'q', q);
end
