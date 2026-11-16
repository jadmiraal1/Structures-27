function M = materials_input()
% MATERIALS_INPUT  Material property table. SI units throughout (Pa, m, kg/m^3).
%
% Stiffness/strength properties are given as [low high] bounds.
%   - Stiffness checks (deflection, twist): run with the LOW bound.
%   - Load-sharing/strength checks: run BOTH bounds, check each member against
%     its LOW strength in whichever run loads it more.
% Thickness and density are single nominal values.
%
% EVERY VALUE BELOW IS A PLACEHOLDER. Replace with your vendor datasheets,
% coupon data, or MMPDS (aluminum). Fill in the 'source' field for each.

%% Carbon twill ply: Hexcel AS4 3K 2x2 twill, 203 gsm, wet layup + vacuum bag,
%% Composite Envisions thin epoxy 3:1. Laminate x-axis = spanwise.
% Bounds span fibre volume fraction Vf = 0.40 (low) to 0.50 (high), typical for
% wet layup under vacuum. Estimated by micromechanics (rule of mixtures for the
% fibre direction, Halpin-Tsai transverse and shear, 0.95 crimp knockdown),
% cross-checked against Hexcel woven-fabric data scaled to the same Vf (<3%).
% Inputs: AS4 E = 231 GPa, rho = 1.79 g/cc (Hexcel AS4 datasheet); epoxy
% E = 293,000 psi (2.02 GPa), mixed 9.06 lb/gal (1.086 g/cc): Composite
% Envisions Thin Epoxy 1159 + medium 3:1 hardener 1161, TDS rev 03/11/2022.
% AS4 transverse/shear moduli 15 GPa (assumed, small effect).
% ROOM TEMPERATURE values. Resin Tg (DMA onset) is only 46 C after the TDS
% post-cure (3 h at 275 F); matrix-dominated properties (G12, compression,
% bonds) drop sharply as a sun-heated black part approaches that.
% Replace with coupon data when available.
M.twill.E1    = [46.5e9 57.8e9];   % warp modulus, Vf 0.40-0.50
M.twill.E2    = [46.5e9 57.8e9];   % fill modulus (balanced twill)
M.twill.G12   = [1.6e9 2.0e9];     % in-plane shear (matrix dominated), Vf 0.40-0.50, room temperature
M.twill.nu12  = [0.05 0.05];
M.twill.t     = 0.25e-3;           % cured ply thickness at Vf 0.45 (0.23-0.28 mm over the range)
M.twill.rho   = 1400;              % kg/m^3 at Vf 0.45
M.twill.eps_t = [0.006 0.008];     % design tensile strain (wet layup, knocked down from ~1.2% failure)
M.twill.eps_c = [0.004 0.006];     % design compressive strain
M.twill.source = 'Micromechanics: Hexcel AS4 datasheet + CE Thin Epoxy 1159/1161 TDS, Vf 0.40-0.50, RT';

%% Aluminum (spar tube option)
M.alum.E      = [68e9 70e9];     % PLACEHOLDER
M.alum.G      = [25.5e9 26.5e9]; % PLACEHOLDER
M.alum.rho    = 2700;            % PLACEHOLDER
M.alum.Fty    = [240e6 240e6];   % PLACEHOLDER - verify for alloy/temper/form (later stage)
M.alum.Ftu    = [290e6 290e6];   % PLACEHOLDER (later stage)
M.alum.source = 'PLACEHOLDER';

%% Nomex honeycomb core (spar web option)
M.nomex.rho   = 48;              % PLACEHOLDER - core density [kg/m^3]
M.nomex.G_L   = [35e6 45e6];     % PLACEHOLDER - ribbon-direction shear modulus (later stage)
M.nomex.G_W   = [20e6 28e6];     % PLACEHOLDER - transverse shear modulus (later stage)
M.nomex.tau_L = [1.0e6 1.2e6];   % PLACEHOLDER - core shear strength (later stage)
M.nomex.source = 'PLACEHOLDER';

%% Balsa (ribs - used in a later stage)
M.balsa.rho   = 150;             % PLACEHOLDER - weigh your actual sheets
M.balsa.source = 'PLACEHOLDER';

%% Adhesive / resin added mass on sandwich faces (areal, kg/m^2 per face)
M.bond_areal = 0.15;             % PLACEHOLDER
end
