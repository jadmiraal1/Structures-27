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
end
