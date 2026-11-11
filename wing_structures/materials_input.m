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

%% Carbon twill ply (woven, so E1 = E2). Laminate x-axis = spanwise.
M.twill.E1    = [45e9 55e9];     % PLACEHOLDER - warp modulus
M.twill.E2    = [45e9 55e9];     % PLACEHOLDER - fill modulus
M.twill.G12   = [3.5e9 4.5e9];   % PLACEHOLDER - in-plane shear modulus
M.twill.nu12  = [0.05 0.05];     % PLACEHOLDER
M.twill.t     = 0.25e-3;         % PLACEHOLDER - cured ply thickness [m]
M.twill.rho   = 1450;            % PLACEHOLDER - cured ply density [kg/m^3]
M.twill.eps_t = [0.008 0.010];   % PLACEHOLDER - allowable tensile strain (used in later stage)
M.twill.eps_c = [0.006 0.008];   % PLACEHOLDER - allowable compressive strain (later stage)
M.twill.source = 'PLACEHOLDER';

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
