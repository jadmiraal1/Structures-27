function S = sandwich_web_section(h_total, face, c, rho_core, cap, bond_areal)
% SANDWICH_WEB_SECTION  Vertical sandwich-panel spar web, optional twill caps.
%
%   h_total    : total spar depth [m] (web + caps)
%   face       : laminate struct from laminate_clt for ONE face sheet
%                (x-axis = spanwise)
%   c          : core thickness [m]
%   rho_core   : core density [kg/m^3]
%   cap        : [] for no caps, or struct with fields
%                  .lam : laminate struct for one cap (x = spanwise)
%                  .w   : cap width [m] (normally the panel thickness)
%   bond_areal : adhesive/resin mass per face per unit area [kg/m^2]
%
%   Bending (EIxx, about the horizontal axis, resists downforce):
%     faces act as two thin vertical plates of height h_w:
%       EI_faces = 2 * Ex_f * t_f * h_w^3 / 12
%     caps sit at +/-(h_total/2 - t_c/2).
%     Core axial stiffness is neglected (honeycomb carries ~no axial load).
%
%   Torsion:
%     with caps   -> closed cell, Bredt-Batho with faces + caps as walls.
%     without caps-> sandwich-plate result GJ = 2 * Gxy_f * t_f * d^2 * h_w
%                    (d = face centroid spacing). This assumes the core is
%                    shear-stiff enough to close the cell; treat as an
%                    upper bound.

tf = face.h;
d  = c + tf;                 % face centroid spacing
t_panel = c + 2*tf;

if isempty(cap)
    tc = 0;
else
    tc = cap.lam.h;
end
h_w = h_total - 2*tc;
if h_w <= 0
    error('sandwich_web_section:depth', 'Caps thicker than available depth.');
end

EI = 2 * face.Ex * tf * h_w^3 / 12;
EA = 2 * face.Ex * tf * h_w;
m  = (2*face.areal + 2*bond_areal) * h_w + rho_core * c * h_w;

if ~isempty(cap)
    Ac  = cap.w * tc;
    arm = h_total/2 - tc/2;
    EI  = EI + 2 * cap.lam.Ex * (Ac*arm^2 + cap.w*tc^3/12);
    EA  = EA + 2 * cap.lam.Ex * Ac;
    m   = m + 2 * cap.lam.areal * cap.w;
    % Closed cell: faces (height h_w, thickness tf) + caps (width d, thickness tc)
    Aenc = (h_w + tc) * d;   % midline box
    GJ   = 4 * Aenc^2 / (2*(h_w + tc)/(face.Gxy*tf) + 2*d/(cap.lam.Gxy*tc));
else
    GJ = 2 * face.Gxy * tf * d^2 * h_w;
end

S.EIxx = EI;
S.EIxx_eff = EI;   % symmetric section, no product term
S.EA   = EA;
S.GJ   = GJ;
S.mass_per_len = m;
S.t_panel = t_panel;
S.h_web = h_w;
end
