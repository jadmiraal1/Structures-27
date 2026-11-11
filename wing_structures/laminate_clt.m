function lam = laminate_clt(ply, angles_deg)
% LAMINATE_CLT  Classical laminate theory for a stack of identical plies.
%
%   ply        : struct with scalar E1, E2, G12, nu12, t, rho
%                (use bound_props first if fields are [low high])
%   angles_deg : ply angles, top to bottom, measured from the laminate
%                x-axis (= spanwise for skins/spar faces). For a woven ply,
%                0 means warp along x.
%
%   Returns A, B, D and effective membrane engineering constants Ex, Ey,
%   Gxy, nuxy, plus thickness h, areal mass, and density.
%
%   Effective constants use the full ABD inverse (laminate free to bend),
%   which equals the A-inverse result for symmetric laminates.

n   = numel(angles_deg);
t   = ply.t;
h   = n * t;
z   = linspace(-h/2, h/2, n+1);

nu21 = ply.nu12 * ply.E2 / ply.E1;
den  = 1 - ply.nu12 * nu21;
Q11 = ply.E1 / den;
Q22 = ply.E2 / den;
Q12 = ply.nu12 * ply.E2 / den;
Q66 = ply.G12;

A = zeros(3); B = zeros(3); D = zeros(3);
for k = 1:n
    th = angles_deg(k) * pi/180;
    m = cos(th); s = sin(th);
    Qb = zeros(3);
    Qb(1,1) = Q11*m^4 + 2*(Q12 + 2*Q66)*m^2*s^2 + Q22*s^4;
    Qb(2,2) = Q11*s^4 + 2*(Q12 + 2*Q66)*m^2*s^2 + Q22*m^4;
    Qb(1,2) = (Q11 + Q22 - 4*Q66)*m^2*s^2 + Q12*(m^4 + s^4);
    Qb(3,3) = (Q11 + Q22 - 2*Q12 - 2*Q66)*m^2*s^2 + Q66*(m^4 + s^4);
    Qb(1,3) = (Q11 - Q12 - 2*Q66)*m^3*s + (Q12 - Q22 + 2*Q66)*m*s^3;
    Qb(2,3) = (Q11 - Q12 - 2*Q66)*m*s^3 + (Q12 - Q22 + 2*Q66)*m^3*s;
    Qb(2,1) = Qb(1,2); Qb(3,1) = Qb(1,3); Qb(3,2) = Qb(2,3);

    A = A + Qb * (z(k+1) - z(k));
    B = B + Qb * (z(k+1)^2 - z(k)^2) / 2;
    D = D + Qb * (z(k+1)^3 - z(k)^3) / 3;
end

ABD = [A B; B D];
abd = inv(ABD);

lam.A = A; lam.B = B; lam.D = D;
lam.h    = h;
lam.Ex   = 1 / (h * abd(1,1));
lam.Ey   = 1 / (h * abd(2,2));
lam.Gxy  = 1 / (h * abd(3,3));
lam.nuxy = -abd(1,2) / abd(1,1);
lam.areal = ply.rho * h;   % kg/m^2 (no resin-rich/bond allowance)
lam.rho   = ply.rho;
lam.angles = angles_deg;

if norm(B, 'fro') > 1e-6 * norm(A, 'fro') * h
    warning('laminate_clt:unsymmetric', ...
        'Laminate is unsymmetric (B ~= 0): bend-extension coupling present.');
end
end
