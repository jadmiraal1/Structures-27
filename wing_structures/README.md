# Wing structures calculator (FSAE main elements)

Hand-calc-level sizing tool for rib-and-spar wing main elements with a
bonded carbon twill skin. Written for MATLAB; also runs in GNU Octave.

## Status

| Stage | Scope | State |
|---|---|---|
| 1 | Section properties: skin shell vs spar options (Al rectangular tube, Nomex sandwich web, sandwich web + twill caps) | Done |
| 1b | Geometry: slice the CAD STL along the span, depth map for spar placement, skin section properties per station | Done |
| 1c | Loads: spanwise load distribution from the Fluent wall export, flap + endplate load lumped at the main element tip | Done (needs fine STL, see below) |
| 2 | Sizing sweep: skin layup x 6061-T6 plate spar gauge x spar position. Half-span beam (symmetry at centre, mount at 206.85 mm): bending, torsion, deflection and twist at limit; plate stress, skin strain, rib spacing for skin buckling and mount loads at ultimate | Done (placeholder materials) |
| 3 | Detailed checks: mount and rib joints (bond shear, bearing), curved-panel skin buckling, plate stability, FEA verification | Planned |

## Running

```matlab
cd wing_structures
addpath('tests'); run_tests; run_tests_geometry   % validation against exact results
main_stage1      % skin vs spar options (generic section)
main_geometry    % slice CAD along span: depth map, skin properties per station
main_loads       % CFD spanwise loads, checks, limit/ultimate factors
main_sizing      % stage 2 sweep (needs the two results files above)
pick_design      % lightest design meeting strength + your stiffness limits
main_ribs        % loft corners, rib spacing from skin buckling, rib layout per skin variant
fea_check        % prediction for the SolidWorks FEA check of the beam model
```

Inputs: `front_wing_config.m` (paths, coordinate transform, load case,
load factor, safety factor, stations) and `materials_input.m`.
Values marked PLACEHOLDER must be replaced before results mean anything.
Outputs go to `results/` (not committed).

## Input data

```
data/front_wing/geometry/main_element.STL   SolidWorks STL of the main element (mm)
data/front_wing/cfd/front_wing              Fluent ASCII wall export of the front wing zone
data/front_wing/cfd/case_notes.txt          axes, symmetry plane, speeds, Fluent force totals
```

Fluent export: File > Export > Solution Data, ASCII, Cell Center, comma
delimited; Static Pressure (gauge, not absolute), X/Y/Z Wall Shear Stress,
X/Y/Z Face Area. STL and CFD exports are git-ignored (large, regenerable).

## Checks built into the pipeline

- Summed face forces must reproduce Fluent Report > Forces (pressure and
  viscous separately). Verifies the export, units and sign convention.
- CFD-to-STL distance by span band: shows whether the CAD and the CFD
  geometry agree.
- Coverage per strip: CFD area classified as main element / STL surface
  area. Values well below 1 mean main element faces were missed and their
  load went into the flap + endplate lump.

## Known issues

- The CFD was run on a slightly different revision of the main element
  than the current CAD: the loft transitions (where the profile rises and
  drops) sit ~15 mm further inboard in the CFD geometry. Sections, root and
  tip match. Between ~210-280 mm and ~500-580 mm span the surfaces differ by
  up to ~6 mm, which is about the main element-to-flap slot gap, so those
  main element faces are not classified (coverage ~0.06). Strips with
  coverage below `cfg.coverage_min` (0.90) are interpolated from their
  neighbours (`fill_low_coverage.m`), and the flap + endplate lump is
  recomputed as whole zone minus main element, so total force and moment
  still match the CFD. Acceptable here because the load is nearly flat
  across those regions; re-running CFD on the current CAD removes it.

## Modelling assumptions

- Skin is bonded to the ribs only, not along the spar. Ribs enforce matching
  vertical deflection and twist at rib stations but transfer no spanwise
  shear, so skin and spar each bend about their own neutral axis:
  EI_total ~ EI_skin + EI_spar.
- Thin-walled sections; airfoil outer mould line used as the skin midline.
- Loads: face force = p*A + tau*|A| (A = Fluent face area vector). Main
  element faces are those within `classify_tol` of the STL skin (end caps
  excluded, endplate-like faces with |n_z| > 0.8 excluded). Everything else
  in the zone (flaps, endplate) is lumped at the main element tip section,
  since the flaps attach through the endplates.
- Design loads: CFD loads x (V_design/V_cfd)^2 x n (limit), x SF (ultimate).
  Deflection checks use limit loads, strength checks ultimate loads.
- Stage 2 beam: projected half span, slope = 0 and torque = 0 at the centre
  plane (symmetric loading), vertical + pitch reaction at the mount
  (`pinned` or `clamped`). Overhang moments are statically determinate.
  The spar is a flat vertical 6061-T6 plate at constant x following the
  profile; skin and plate share bending by EI (ribs only, each about its own
  axis); the plate has negligible GJ, so the skin carries torsion. Elastic
  axis = EI-weighted mean of skin shear centre and plate.
- Flap + endplate load enters the main element tip as one force + couple
  (flaps cantilevered from the endplates).
- Rib spacing from skin buckling uses a flat wide column between ribs,
  a_max = pi*sqrt(D11/N), ignoring skin curvature and rib fixity:
  conservative.
- Loft corners: a skin force N turning through an angle dtheta pushes
  sideways with N*dtheta per unit chord; corners above cfg.kink_rib_min
  get their own rib.
- Material properties as [low high] bounds: low for stiffness checks, both
  bounds for load sharing, low strength for allowables.
- Honeycomb core carries no axial load. Uncapped sandwich-web torsion uses
  the sandwich-plate result (upper bound, assumes a shear-stiff core).

## Files

| File | Purpose |
|---|---|
| `materials_input.m` | Material table ([low high] bounds, sources) |
| `bound_props.m` | Pick low or high bound |
| `laminate_clt.m` | Classical laminate theory, effective Ex, Gxy |
| `thinwall_closed_section.m` | Modulus-weighted EI, GJ (Bredt-Batho), mass for closed thin-wall sections |
| `sandwich_web_section.m` | Sandwich-panel spar web, optional caps |
| `rect_tube_midline.m` | Rectangular tube contour |
| `load_airfoil.m` | Selig .dat loader (or NACA 4-digit placeholder) |
| `main_stage1.m` | Stage 1 driver |
| `front_wing_config.m` | Front wing inputs |
| `read_stl.m`, `read_fluent_ascii.m` | Readers |
| `transform_tris.m`, `load_main_element.m` | STL to CFD frame, end-cap removal |
| `point_triangle_distance.m`, `distance_to_mesh.m` | Exact point-to-mesh distance |
| `slice_mesh.m`, `section_depth.m`, `polygon_centroid.m` | Spanwise slicing, thickness distribution |
| `face_forces.m`, `integrate_strips.m` | Face forces, spanwise strip integration |
| `fill_low_coverage.m` | Interpolate strips where main element faces were missed |
| `shear_center_closed.m` | Shear centre of a closed thin-walled cell |
| `station_props.m` | Skin + plate properties per station |
| `halfspan_beam.m` | Half-span bending/torsion with centre symmetry and mount |
| `main_sizing.m` | Stage 2 sizing sweep driver |
| `pick_design.m` | Picks the lightest design meeting strength and the stiffness limits in the config |
| `main_ribs.m` | Rib analysis for the chosen design: loft corner loads, allowable rib spacing, rib layout, skin variants |
| `zoned_props.m` | Section properties with a skin layup that changes along the span |
| `kink_corners.m` | Finds loft corners and the sideways skin load at each |
| `rib_layout.m` | Places ribs so no bay exceeds its allowable spacing |
| `fea_check.m` | Skin-only check case and predictions for an independent FEA comparison |
| `main_geometry.m`, `main_loads.m` | Geometry and loads drivers |
| `tests/` | Validation tests |

## References

- Megson, *Aircraft Structures for Engineering Students* (thin-wall sections, shear flow, boom idealization)
- Bruhn, *Analysis and Design of Flight Vehicle Structures* (buckling)
- Hexcel, *HexWeb Honeycomb Sandwich Design Technology*; Zenkert, *An Introduction to Sandwich Construction*

## Design-review figures (`review_figures.m`)

MATLAB script that reads the saved results and writes three PNG figures for slides to
`results/review`: loads and bending moment, the section at the mount, and the rib layout
against the allowable spacing. It does no sizing of its own.

1. From `wing_structures`, run `main_ribs` (it needs the saved results of `main_geometry`
   and `main_loads`). Re-run it whenever `front_wing_config.m` changes.
2. Run `review_figures`.
3. To show a different skin variant, change `variant` at the top of the script to any text
   contained in the variant name (default `'inboard of 350'`).

The script also prints the rib stations, each bay against its allowable spacing, mass,
deflection and twist. Lofted regions are shaded from what `main_ribs` saved. With
`cfg.rib_no_loft = true` (see `front_wing_config.m`) no rib is placed inside a loft;
`main_ribs` flags any skin variant that cannot span a loft without one.

`review_figures.py` is an earlier Python version of the same figures (needs numpy, scipy,
matplotlib; run `python review_figures.py`). It is not needed if you use the MATLAB script.
