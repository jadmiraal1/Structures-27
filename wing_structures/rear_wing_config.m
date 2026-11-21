function cfg = rear_wing_config()
% REAR_WING_CONFIG  Inputs for the rear wing main element.
% Coordinates follow the CFD: x forward, y up, z right (driver's view).
% Paths are relative to the wing_structures folder.
% Make this the active wing with:  wing_config('rear_wing')

cfg.name = 'rear_wing';

%% Geometry (SolidWorks STL of the main element only, left half)
cfg.stl       = fullfile('data', 'rear_wing', 'geometry', 'RW_mainelement.STL');
cfg.stl_scale = 1e-3;                 % STL exported in mm
% STL part frame -> CFD frame:  x_cfd = R * x_stl + t   [m]
% Found by fitting the STL to the CFD wall faces (translation only, same
% rotation as the front wing): 94% of the STL triangles lie within 2 mm of a
% CFD face, median distance 0.012 mm, so CAD and CFD are the same revision.
% If you re-export the STL in the CFD coordinate system, set R = eye(3), t = [0 0 0].
cfg.stl_R = [-1 0 0; 0 1 0; 0 0 -1];  % 180 deg about vertical
cfg.stl_t = [-1.322203 0.063738 0.000730];

%% CFD (Fluent ASCII export, whole rear wing zone, half car)
cfg.cfd     = fullfile('data', 'rear_wing', 'cfd', 'rear_wing');
cfg.z_sym   = 0;                      % symmetry plane z [m]
cfg.side    = -1;                     % modelled half is z < z_sym (left)
% Fluent Report > Forces for the same zone [N], used to check the integration.
% NOT PROVIDED YET: main_loads prints the integrated totals and warns. Enter
% Fluent's numbers here (pressure and viscous, [Fx Fy Fz]) to turn the check on.
cfg.fluent_pressure = [NaN NaN NaN];
cfg.fluent_viscous  = [NaN NaN NaN];

%% Load case
cfg.V_cfd    = 15;                    % m/s, speed of the CFD run (peak pressure 144 Pa matches q at 15 m/s)
cfg.V_design = 25;                    % m/s, predicted max speed
cfg.n_load   = 2.0;                   % limit load factor (same basis as the front wing)
cfg.SF       = 1.5;                   % safety factor, ultimate = SF * limit

%% Face classification
cfg.classify_tol       = 2e-3;        % m, max distance from main element STL surface
cfg.max_lateral_normal = 0.8;         % faces with |n_z| above this are endplate-like
cfg.coverage_min = 0.90;
% How the flap loads reach the main element tip (main_loads):
%   'clamped' force plus the full couple about the tip: flaps cantilevered from the endplate
%   'pinned'  force plus the pitch couple only: flaps simply supported at the endplates, and
%             the endplate's tipping moment carried by axial forces in the elements (frame)
% ASSUMPTION: the rear flaps run endplate to endplate and are bolted to both, so 'pinned'.
cfg.flap_tip_fixity = 'pinned';

%% Spanwise discretisation
cfg.strip_ds = 0.020;                 % m, load strip width
cfg.stations = 0.010:0.025:0.385;     % m from centreline, geometry slices (half span 0.400)

%% Skin layup for section properties (materials from materials_input.m)
cfg.skin_layup = [0 45 0];            % used by main_geometry

%% Skin strain model (beam_skin_strains)
% The rear section sits at 6 deg (centre) to 13 deg (tip) to the horizontal, so its
% neutral axis is far from horizontal: use unsymmetric bending. Drag is 40% of the
% downforce on this wing (80% for the flap + endplate load), so include its bending.
cfg.strain_model = 'unsymmetric';     % 'horizontal' (front wing method) or 'unsymmetric'
cfg.drag_bending = true;

%% Structure / sizing (main_sizing)
cfg.mount_s   = 7 * 0.0254;           % m from centre plane: mount face 7 in from the midplane, one per side
cfg.mount_bc  = 'pinned';             % 'pinned' (conservative) or 'clamped'
cfg.bond_gap  = 0.5e-3;               % m, spar-to-skin clearance each side
% Candidate skin layups (twill plies, deg from spanwise)
cfg.skin_options = {[0], [0 0], [45 45], [0 45 0], [45 0 45], [0 45 45 0]};
% Candidate 6061-T6 plate thicknesses [in] (common sheet/plate gauges)
cfg.plate_gauges_in = [0.063 0.080 0.090 0.100 0.125 0.160 0.190 0.250];
% Candidate spar plane positions, distance aft of the leading edge [m]
% (flat plate -> constant x plane; chord is 358 mm, max thickness 44 mm at 70-76 mm)
cfg.spar_aft_LE = 0.045:0.015:0.120;

%% Rib analysis (main_ribs): the chosen design and skin variants to compare
cfg.rib_design = struct('spar_aft', 0.075, 'plate_in', 0.063, 'base_layup', [0 45 0]);
% Each variant: a name and zones of extra plies ([] = base layup everywhere).
% zone fields: s = [start end] in m from centre plane, layup = local layup.
cfg.rib_variants = { ...
    struct('name', '[0 45 0] everywhere',          'zones', []), ...
    struct('name', '[0 0 45 0 0] inboard of 250',  'zones', struct('s', [0 0.250], 'layup', [0 0 45 0 0])), ...
    struct('name', '[0 45 45 0] everywhere',       'zones', struct('s', [0 1], 'layup', [0 45 45 0])), ...
    struct('name', '[0 0 45 0 0] everywhere',      'zones', struct('s', [0 1], 'layup', [0 0 45 0 0])) };
cfg.kink_curv_min = 3;       % rad/m: skin counts as turning above this
cfg.kink_rib_min  = 1000;    % N per m of chord at ultimate: corner needs its own rib above this
cfg.rib_merge     = 0.015;   % m: a corner this close to the mount rib uses the mount rib
cfg.rib_no_loft     = true;        % no ribs inside lofted regions (none are found on this wing: it twists
                                   %   along the span, the section does not rise or drop as a whole)
cfg.loft_slope_min  = deg2rad(1);  % rad: both skin lines steeper than this in front view = lofted
cfg.rib_loft_margin = 0.003;       % m: rib centre this far outside the loft end
cfg.review_variant = 'inboard of 250';        % variant shown by review_figures (text contained in its name)
cfg.rib_K         = 0.7;     % skin buckling effective length factor: 1 = ribs act as pins (conservative),
                             %   ~0.7 = filleted/tabbed ribs partly clamp the skin, 0.5 = fully clamped
cfg.rib_check_load = 'ultimate';  % 'ultimate' (no skin buckling below ultimate), 'limit', or 'aero'
cfg.rib_areal     = 0.45;    % kg/m^2 of rib face. PLACEHOLDER: 3 mm balsa at 150 kg/m^3, no lightening holes
cfg.rib_fillet    = 0.010;   % kg per m of rib perimeter (fillets both sides). PLACEHOLDER

%% Stiffness requirements at LIMIT load (used by pick_design). YOUR CALL.
cfg.max_tip_deflection_mm = NaN;
cfg.max_tip_twist_deg     = NaN;
end
