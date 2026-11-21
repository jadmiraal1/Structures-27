function cfg = front_wing_config()
% FRONT_WING_CONFIG  Inputs for the front wing main element.
% Coordinates follow the CFD: x forward, y up, z right (driver's view).
% Paths are relative to the wing_structures folder.

cfg.name = 'front_wing';

%% Geometry (SolidWorks STL of the main element only)
cfg.stl       = fullfile('data', 'front_wing', 'geometry', 'main_element.STL');
cfg.stl_scale = 1e-3;                 % STL exported in mm
% STL part frame -> CFD frame:  x_cfd = R * x_stl + t   [m]
% Found by matching the inboard leading edge; inboard agreement ~0.15 mm.
% If you re-export the STL in the CFD coordinate system, set R = eye(3), t = [0 0 0].
cfg.stl_R = [-1 0 0; 0 1 0; 0 0 -1];  % 180 deg about vertical
cfg.stl_t = [0.921962875 0 0];

%% CFD (Fluent ASCII export, whole front wing zone, half car)
cfg.cfd     = fullfile('data', 'front_wing', 'cfd', 'front_wing');
cfg.z_sym   = 0;                      % symmetry plane z [m]
cfg.side    = -1;                     % modelled half is z < z_sym (left)
% Fluent Report > Forces for the same zone [N], used to check the integration
cfg.fluent_pressure = [-8.9581374  -74.684269    0.39330499];
cfg.fluent_viscous  = [-0.69002864 -0.0032465998 -0.02527921];

%% Load case
cfg.V_cfd    = 15;                    % m/s, speed of the CFD run
cfg.V_design = 25;                    % m/s, predicted max speed
cfg.n_load   = 2.0;                   % PLACEHOLDER limit load factor (braking dive, bumps)
cfg.SF       = 1.5;                   % PLACEHOLDER safety factor, ultimate = SF * limit

%% Face classification
cfg.classify_tol       = 2e-3;        % m, max distance from main element STL surface
cfg.max_lateral_normal = 0.8;         % faces with |n_z| above this are endplate-like
% Strips where less than this fraction of the main element surface was
% found are replaced by interpolation from neighbouring strips (the CFD
% geometry's loft transitions sit ~15 mm inboard of the current CAD).
cfg.coverage_min = 0.90;
% How the flap loads reach the main element tip (main_loads):
%   'clamped' force plus the full couple about the tip: flaps cantilevered from the endplate
%   'pinned'  force plus the pitch couple only: flaps simply supported at the endplates, and
%             the endplate's tipping moment carried by axial forces in the elements (frame)
cfg.flap_tip_fixity = 'clamped';

%% Spanwise discretisation
cfg.strip_ds = 0.020;                 % m, load strip width
cfg.stations = 0.010:0.025:0.640;     % m from centreline, geometry slices

%% Skin layup for section properties (materials from materials_input.m)
cfg.skin_layup = [0 45 0];            % PLACEHOLDER (used by main_geometry)

%% Structure / sizing (main_sizing)
cfg.mount_s   = 0.20685;              % m from centre plane, one mount per side
cfg.mount_bc  = 'pinned';             % 'pinned' (conservative) or 'clamped'
cfg.bond_gap  = 0.5e-3;               % m, spar-to-skin clearance each side
% Candidate skin layups (twill plies, deg from spanwise)
cfg.skin_options = {[0], [0 0], [45 45], [0 45 0], [45 0 45], [0 45 45 0]};
% Candidate 6061-T6 plate thicknesses [in] (common sheet/plate gauges)
cfg.plate_gauges_in = [0.063 0.080 0.090 0.100 0.125 0.160 0.190 0.250];
% Candidate spar plane positions, distance aft of the leading edge [m]
% (flat plate -> constant x plane; chord is 321-331 mm)
cfg.spar_aft_LE = 0.050:0.015:0.125;

%% Rib analysis (main_ribs): the chosen design and skin variants to compare
cfg.rib_design = struct('spar_aft', 0.080, 'plate_in', 0.063, 'base_layup', [0 45 0]);
% Each variant: a name and zones of extra plies ([] = base layup everywhere).
% zone fields: s = [start end] in m from centre plane, layup = local layup.
cfg.rib_variants = { ...
    struct('name', '[0 45 0] everywhere',          'zones', []), ...
    struct('name', '[0 0 45 0 0] inboard of 275',  'zones', struct('s', [0 0.275], 'layup', [0 0 45 0 0])), ...
    struct('name', '[0 0 45 0 0] inboard of 350',  'zones', struct('s', [0 0.350], 'layup', [0 0 45 0 0])), ...
    struct('name', '[0 0 45 0 0] inboard of 450',  'zones', struct('s', [0 0.450], 'layup', [0 0 45 0 0])), ...
    struct('name', '[0 45 45 0] everywhere',       'zones', struct('s', [0 1], 'layup', [0 45 45 0])), ...
    struct('name', '[0 0 45 0 0] everywhere',      'zones', struct('s', [0 1], 'layup', [0 0 45 0 0])) };
cfg.kink_curv_min = 3;       % rad/m: skin counts as turning above this
cfg.kink_rib_min  = 1000;    % N per m of chord at ultimate: corner needs its own rib above this
cfg.rib_merge     = 0.015;   % m: a corner this close to the mount rib uses the mount rib
cfg.rib_no_loft     = true;        % true: no ribs inside lofted regions (where the wing rises or drops in
                                   %   front view). A rib is placed just outside each end of each loft instead,
                                   %   and corners no longer get their own rib. false: previous behaviour.
cfg.loft_slope_min  = deg2rad(1);  % rad: skin lines steeper than this in front view count as lofted
cfg.rib_loft_margin = 0.003;       % m: rib centre this far outside the loft end (rib on constant section)
cfg.review_variant = 'inboard of 350';   % variant shown by review_figures (text contained in its name)
cfg.rib_K         = 0.7;     % skin buckling effective length factor: 1 = ribs act as pins (conservative),
                             %   ~0.7 = filleted/tabbed ribs partly clamp the skin, 0.5 = fully clamped
cfg.rib_check_load = 'ultimate';  % 'ultimate' (no skin buckling below ultimate), 'limit', or
                                  % 'aero' (max-speed aero load only, no factors - for comparing with past wings)
cfg.rib_areal     = 0.45;    % kg/m^2 of rib face. PLACEHOLDER: 3 mm balsa at 150 kg/m^3, no lightening holes
cfg.rib_fillet    = 0.010;   % kg per m of rib perimeter (fillets both sides). PLACEHOLDER

%% Stiffness requirements at LIMIT load (used by pick_design). YOUR CALL.
% NaN = not set; pick_design then lists the trade-off designs instead.
cfg.max_tip_deflection_mm = NaN;
cfg.max_tip_twist_deg     = NaN;
end
