# Wing structures calculator (FSAE main elements)

Hand-calc-level sizing tool for rib-and-spar wing main elements with a
bonded carbon twill skin. Written for MATLAB; also runs in GNU Octave.

## Status

| Stage | Scope | State |
|---|---|---|
| 1 | Section properties: skin shell vs spar options (Al rectangular tube, Nomex sandwich web, sandwich web + twill caps) | Done |
| 2 | Coupled two-beam FE (skin + spar linked at ribs), CFD loads, deflection, twist, rib transfer loads | Planned |
| 3 | Strength and buckling checks (strain allowables, skin panel buckling, sandwich failure modes, rib/mount loads) | Planned |

## Running

```matlab
cd wing_structures
addpath('tests'); run_tests   % validation against closed-form results
main_stage1                   % skin vs spar comparison
```

All inputs are at the top of `main_stage1.m` and in `materials_input.m`.
Values marked PLACEHOLDER must be replaced before results mean anything.

## Modelling assumptions

- Skin is bonded to the ribs only, not along the spar. Ribs enforce matching
  vertical deflection and twist at rib stations but transfer no spanwise
  shear, so skin and spar each bend about their own neutral axis:
  EI_total ~ EI_skin + EI_spar.
- Thin-walled sections; airfoil outer mould line used as the skin midline.
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
| `tests/` | Validation tests |

## References

- Megson, *Aircraft Structures for Engineering Students* (thin-wall sections, shear flow, boom idealization)
- Bruhn, *Analysis and Design of Flight Vehicle Structures* (buckling)
- Hexcel, *HexWeb Honeycomb Sandwich Design Technology*; Zenkert, *An Introduction to Sandwich Construction*
