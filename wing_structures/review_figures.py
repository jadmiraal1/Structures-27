#!/usr/bin/env python3
"""REVIEW_FIGURES  Design-review figures for the front wing main element.

Reads the results written by the MATLAB pipeline (main_geometry, main_loads,
main_ribs) and writes three PNG figures:

  fig1_loads.png       spanwise downforce at limit load, bending moment at ultimate
  fig2_section.png     section at the mount: skin, spar plate, neutral axis
  fig3_rib_layout.png  rib stations on the planform, and bay length vs allowable

Usage (from wing_structures):
  python review_figures.py
  python review_figures.py --variant "inboard of 350"
  python review_figures.py --variant "inboard of 350" --saved-K 1.0 --relayout-K 0.7

Nothing here sizes the wing: every number comes from the results files.
The one exception is --relayout-K, for when front_wing_ribs.mat was saved
with a different rib fixity K than the one you want to show. It rescales the
saved allowable spacing (a ~ 1/K) and re-places the ribs with a line-for-line
port of rib_layout.m. Prefer re-running main_ribs and omitting the flag.

Needs: numpy, scipy, matplotlib.
"""
import argparse
import io
import os
import re

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

INK, GRID = "#222222", "#D0D0D0"
BLUE, ORANGE, GREY, TINT_A, TINT_B = "#1F4E79", "#C55A11", "#7F7F7F", "#DCE6F1", "#F2F2F2"
plt.rcParams.update({
    "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "axes.edgecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.8, "legend.frameon": False,
    "figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight", "font.family": "DejaVu Sans",
})


# ----------------------------------------------------------------- inputs
def load_loads_csv(path):
    """Spanwise strip loads at the CFD speed, plus the header's factors and tip lump."""
    head = [r for r in open(path) if r.startswith("#")]
    body = [r for r in open(path) if not r.startswith("#")]
    d = np.genfromtxt(io.StringIO("".join(body)), delimiter=",", names=True)
    txt = " ".join(head)
    k_lim = float(re.search(r"k_limit=(\d+\.\d+)", txt).group(1))
    k_ult = float(re.search(r"k_ult=(\d+\.\d+)", txt).group(1))
    m = re.search(r"lump at \[([^\]]+)\] m: F = \[([^\]]+)\] N, M = \[([^\]]+)\]", txt)
    r_ep, F_ep, M_ep = (np.array(g.split(), float) for g in m.groups())
    return d, k_lim, k_ult, r_ep, F_ep, M_ep


def hogging_moment(s, d, F_ep, M_ep, s_tip, s_mount):
    """Bending moment (+ = top skin in tension) from statics: pinned mount,
    symmetric centre plane. Same quantity as halfspan_beam's M_hog."""
    s_i, Fy = d["s_mid_m"], d["fy_N_per_m"] * d["ds_m"]
    total_down = -(Fy.sum() + F_ep[1])
    M = np.zeros_like(s)
    for j, s0 in enumerate(s):
        o = s_i > s0
        M[j] = np.sum(-Fy[o] * (s_i[o] - s0)) + (-F_ep[1]) * (s_tip - s0) - M_ep[0]
        if s0 < s_mount:
            M[j] -= total_down * (s_mount - s0)
    return M


# ------------------------------------------------- port of rib_layout.m
def rib_layout(s, a_max, fixed, s_tip):
    fixed = np.unique(np.r_[np.atleast_1d(fixed), s_tip])

    def bmin(x0, x1):
        return 0.0 if (x0 < 0 < x1) else min(abs(x0), abs(x1))

    def amin(x0, x1):
        lo, hi = bmin(x0, x1), max(abs(x0), abs(x1))
        inside = a_max[(s >= lo - 1e-12) & (s <= hi + 1e-12)]
        return min(np.r_[inside, np.interp([lo, hi], s, a_max)])

    def march(x0, x1, f):
        p, x = [], x0
        while x < x1 - 1e-9:
            lo, hi = x, x1
            if hi - x <= f * amin(x, hi) + 1e-12:
                x = hi
            else:
                for _ in range(40):
                    mid = (lo + hi) / 2
                    if mid - x <= f * amin(x, mid) + 1e-12:
                        lo = mid
                    else:
                        hi = mid
                if lo - x < 1e-6:
                    raise ValueError("allowable spacing too small")
                x = lo
            p.append(x)
        return p

    r1 = fixed[0]
    for n in range(1, 501):
        p = np.linspace(-r1, r1, n + 1)
        if all(p[q + 1] - p[q] <= amin(p[q], p[q + 1]) + 1e-12 for q in range(n)):
            break
    ribs = list(p[p >= -1e-12])
    for x0, x1 in zip(fixed[:-1], fixed[1:]):
        n = len(march(x0, x1, 1.0))
        f_lo, f_hi = 0.0, 1.0
        for _ in range(40):
            f = (f_lo + f_hi) / 2
            try:
                ok = len(march(x0, x1, f)) <= n
            except ValueError:
                ok = False
            f_lo, f_hi = (f_lo, f) if ok else (f, f_hi)
        ribs += march(x0, x1, f_hi)
    ribs = np.unique(np.round(np.array(ribs) * 1e7) / 1e7)
    ribs[np.abs(ribs) < 1e-9] = 0
    full = np.unique(np.r_[-ribs[ribs > 0], ribs])
    bays = [(x0, x1, x1 - x0, amin(x0, x1)) for x0, x1 in zip(full[:-1], full[1:]) if x1 > 0]
    return ribs, np.array(bays)


# ---------------------------------------------------------------- figures
def fig_loads(d, k_lim, k_ult, F_ep, M_ep, s_tip, s_mount, out):
    s = np.linspace(0, s_tip, 600)
    M = hogging_moment(s, d, F_ep, M_ep, s_tip, s_mount) * k_ult
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.4), sharex=True)
    edges = np.r_[d["s_mid_m"] - d["ds_m"] / 2, d["s_mid_m"][-1] + d["ds_m"][-1] / 2] * 1e3
    w = -d["fy_N_per_m"] * k_lim
    a1.stairs(w, edges, color=BLUE, lw=2, fill=True, alpha=0.18)
    a1.stairs(w, edges, color=BLUE, lw=2)
    main_total = -(d["fy_N_per_m"] * d["ds_m"]).sum() * k_lim
    a1.annotate(f"Flaps + endplate\n{-F_ep[1] * k_lim:.0f} N at the tip", xy=(s_tip * 1e3, 0.55 * w.max()),
                xytext=(s_tip * 1e3 - 150, 0.22 * w.max()), color=ORANGE, ha="center",
                arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2))
    a1.text(15, 0.12 * w.max(), f"Main element: {main_total:.0f} N per half wing", color=BLUE, va="bottom")
    a1.set_ylabel("Downforce per span [N/m]")
    a1.set_ylim(0, 1.15 * w.max())
    a1.set_title(f"Limit load (CFD x {k_lim:.2f})", loc="left")

    a2.plot(s * 1e3, M, color=BLUE, lw=2.5)
    a2.axhline(0, color=INK, lw=0.8)
    i = np.argmax(M)
    a2.annotate(f"{M[i]:.0f} N m at the mount", xy=(s[i] * 1e3, M[i]), xytext=(s[i] * 1e3 + 60, M[i] * 0.93),
                color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=1))
    a2.set_ylabel("Bending moment [N m]")
    a2.set_xlabel("Span from centre plane [mm]")
    a2.set_title(f"Ultimate load (CFD x {k_ult:.2f}); positive = top skin in tension", loc="left")
    for a in (a1, a2):
        a.axvline(s_mount * 1e3, color=GREY, ls="--", lw=1.2)
    a1.text(s_mount * 1e3 + 6, 1.06 * w.max(), "mount", color=GREY)
    a2.set_xlim(0, s_tip * 1e3 + 5)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return M[i], main_total, -F_ep[1] * k_lim


def fig_section(outline, spar_aft, plate_t, skin_t, gap, label, out):
    xy = np.asarray(outline, float)
    if np.allclose(xy[0], xy[-1]):
        xy = xy[:-1]
    p2 = np.roll(xy, -1, axis=0)
    L = np.hypot(*(p2 - xy).T)
    mid = (xy + p2) / 2
    yc = np.sum(L * mid[:, 1]) / L.sum()                # skin centroid = neutral axis height
    x = (xy[:, 0].max() - xy[:, 0]) * 1e3               # mm behind the leading edge
    y = (xy[:, 1] - yc) * 1e3
    near = np.abs(x - spar_aft * 1e3) < 2.0
    y_hi, y_lo = y[near & (y > 0)].mean(), y[near & (y < 0)].mean()
    clr = (skin_t + gap) * 1e3

    fig, ax = plt.subplots(figsize=(10, 3.4))
    ax.plot(np.r_[x, x[0]], np.r_[y, y[0]], color=INK, lw=2.5)
    ax.axhline(0, color=GREY, ls="--", lw=1.2)
    ax.text(0.55 * x.max(), 1.0, "neutral axis", color=GREY, ha="left", va="bottom")
    ax.add_patch(Rectangle((spar_aft * 1e3 - 1.2, y_lo + clr), 2.4, (y_hi - y_lo) - 2 * clr, color=ORANGE, zorder=3))
    ax.annotate(f"Spar plate, {spar_aft * 1e3:.0f} mm behind LE\n{(y_hi - y_lo) - 2 * clr:.0f} mm tall",
                xy=(spar_aft * 1e3 + 1.5, (y_hi + y_lo) / 2), xytext=(spar_aft * 1e3 + 24, -9), color=ORANGE,
                va="center", arrowprops=dict(arrowstyle="-", color=ORANGE, lw=1))
    ib = np.argmin(y)
    ax.annotate(f"Bottom skin: {-y[ib]:.1f} mm below the neutral axis\n(compression under downforce)",
                xy=(x[ib], y[ib]), xytext=(x[ib] + 70, y[ib] - 13), color=BLUE,
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.5))
    top = (y > 0) & (x > 0.35 * x.max()) & (x < 0.65 * x.max())
    xt = x[top][np.argmin(np.abs(x[top] - 0.5 * x.max()))]
    yt = y[top][np.argmin(np.abs(x[top] - 0.5 * x.max()))]
    ax.annotate(f"Top skin: about {yt:.0f} mm above (tension)", xy=(xt, yt), xytext=(xt + 10, yt + 14),
                color=BLUE, arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.5))
    ax.set_aspect("equal")
    ax.set_xlim(-8, x.max() + 8)
    ax.set_ylim(y.min() - 22, y.max() + 20)
    ax.set_xlabel("Distance behind leading edge [mm]")
    ax.set_ylabel("Height [mm]")
    ax.set_title(f"Section at {label}: chord {x.max():.0f} mm", loc="left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_ribs(s, a, ribs, bays, s_mount, fixed, zones, chord, spar_aft, K, title, out, lofts=None):
    s_tip = ribs.max()
    others = [x for x in fixed if abs(x - s_mount) > 1e-6 and abs(x - s_tip) > 1e-6]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.8), sharex=True, gridspec_kw={"height_ratios": [1, 1.25]})
    # planform
    a1.grid(False)
    for (z0, z1, lab, col) in zones:
        a1.add_patch(Rectangle((z0, 0), z1 - z0, chord, color=col, zorder=0))
        a1.text((z0 + z1) / 2, chord + 48, lab, ha="center", va="top", color=INK, zorder=6,
                bbox=dict(facecolor=col, edgecolor="none", pad=4))
    for (l0, l1) in ([] if lofts is None else lofts):
        a1.add_patch(Rectangle((l0 * 1e3, 0), (l1 - l0) * 1e3, chord, facecolor="none", edgecolor=GREY,
                               hatch="///", lw=0, zorder=1))
        a1.text(0.5 * (l0 + l1) * 1e3, chord * 0.62, "loft", ha="center", va="center", color=INK, zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", pad=3))
    a1.plot([0, s_tip * 1e3], [spar_aft, spar_aft], color=ORANGE, lw=3)
    a1.text(0.5 * (ribs[0] + ribs[1]) * 1e3, spar_aft + 16, "spar plate", color=ORANGE, va="top", ha="center")
    for r in ribs:
        tip, mnt = abs(r - s_tip) < 1e-6, abs(r - s_mount) < 1e-6
        crn = any(abs(r - x) < 1e-6 for x in others)
        col = INK if tip else (ORANGE if mnt else (BLUE if crn else GREY))
        a1.plot([r * 1e3] * 2, [0, chord], color=col, lw=4 if (tip or mnt or crn) else 2)
        a1.text(r * 1e3, -10, f"{r * 1e3:.0f}", ha="center", va="bottom", color=col, fontsize=11.5)
    a1.text(s_mount * 1e3, chord + 12, "mount rib", ha="center", va="top", color=ORANGE)
    if others:
        a1.text(others[0] * 1e3 + 4, chord + 12, "loft-end ribs (blue)" if lofts is not None else "corner rib",
                ha="left", va="top", color=BLUE)
    a1.text(s_tip * 1e3, chord + 12, "endplate", ha="right", va="top", color=INK)
    a1.set_ylim(chord + 95, -55)
    a1.set_yticks([0, chord])
    a1.set_yticklabels(["LE", "TE"])
    a1.set_title(title, loc="left")
    for sp in ("left", "bottom"):
        a1.spines[sp].set_visible(False)
    a1.tick_params(length=0)
    # spacing check
    a2.plot(s * 1e3, np.minimum(a, 0.5) * 1e3, color=BLUE, lw=2.5, label=f"Allowable spacing (skin buckling, ultimate, K = {K:g})")
    first = True
    for x0, x1, L, _ in bays:
        a2.plot([max(x0, 0) * 1e3, x1 * 1e3], [L * 1e3] * 2, color=ORANGE, lw=5, solid_capstyle="butt",
                label="Actual bay length" if first else None)
        first = False
    for r in ribs:
        a2.axvline(r * 1e3, color=GRID, lw=0.8, zorder=0)
    a2.set_ylim(0, 200)
    a2.set_xlim(-5, s_tip * 1e3 + 5)
    a2.set_ylabel("Rib spacing [mm]")
    a2.set_xlabel("Span from centre plane [mm]")
    a2.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


# ------------------------------------------------------------------- main
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=os.path.join(here, "results"))
    ap.add_argument("--out", default=None, help="output folder (default: <results>/review)")
    ap.add_argument("--name", default="front_wing")
    ap.add_argument("--variant", default="inboard of 350", help="text contained in the variant name")
    ap.add_argument("--mount", type=float, default=0.20685, help="mount position [m] (cfg.mount_s)")
    ap.add_argument("--spar-aft", type=float, default=0.080, help="spar plane behind LE [m]")
    ap.add_argument("--plate-in", type=float, default=0.063)
    ap.add_argument("--bond-gap", type=float, default=0.5e-3)
    ap.add_argument("--ply-t", type=float, default=0.25e-3)
    ap.add_argument("--saved-K", type=float, default=None, help="rib_K the ribs .mat was saved with")
    ap.add_argument("--relayout-K", type=float, default=None, help="rib_K to show instead (needs --saved-K)")
    ap.add_argument("--K", type=float, default=0.7, help="rib_K label when not re-laying out")
    ap.add_argument("--fixed", type=float, nargs="+", default=None,
                    help="with --relayout-K: forced rib stations [m] instead of the saved ones (first = mount)")
    ap.add_argument("--loft", type=float, nargs="+", default=None,
                    help="lofted regions to shade, start end [start end ...] in m (default: read from the ribs .mat)")
    a = ap.parse_args()
    out = a.out or os.path.join(a.results, "review")
    os.makedirs(out, exist_ok=True)
    f = lambda ext: os.path.join(a.results, f"{a.name}_{ext}")

    geo = sio.loadmat(f("geometry.mat"), squeeze_me=True, struct_as_record=False)
    V = np.atleast_1d(sio.loadmat(f("ribs.mat"), squeeze_me=True, struct_as_record=False)["V"])
    d, k_lim, k_ult, r_ep, F_ep, M_ep = load_loads_csv(f("loads.csv"))
    v = next(x for x in V if a.variant.lower() in str(x.name).lower())
    s, allow = np.asarray(v.s, float), np.asarray(v.a, float)
    ribs, bays = np.atleast_1d(v.ribs).astype(float), np.atleast_2d(v.bays).astype(float)
    fixed = np.atleast_1d(v.fixed).astype(float) if a.fixed is None else np.array(a.fixed, float)
    s_tip, K = ribs.max(), a.K
    if a.relayout_K:
        assert a.saved_K, "--relayout-K needs --saved-K"
        allow, K = allow * a.saved_K / a.relayout_K, a.relayout_K
        ribs, bays = rib_layout(s, allow, fixed, s_tip)
    lofts = None
    if a.loft:
        lofts = np.array(a.loft, float).reshape(-1, 2)
    elif hasattr(v, "lofts") and np.size(v.lofts):
        lofts = np.asarray(v.lofts, float).reshape(-1, 2)

    M_peak, F_main, F_tip = fig_loads(d, k_lim, k_ult, F_ep, M_ep, s_tip, a.mount, os.path.join(out, "fig1_loads.png"))

    st = np.atleast_1d(geo["G"].s)
    k = int(np.argmin(np.abs(st - a.mount)))
    n_base = len(str(V[0].name).split("]")[0].strip("[ ").split())
    name = str(v.name)
    zone_end = float(re.search(r"inboard of (\d+)", name).group(1)) if "inboard of" in name else None
    n_in = len(name.split("]")[0].strip("[ ").split())
    skin_t = (n_in if (zone_end is None or st[k] * 1e3 <= zone_end) else n_base) * a.ply_t
    fig_section(geo["outlines"][k], a.spar_aft, a.plate_in * 0.0254, skin_t, a.bond_gap,
                f"the mount ({st[k] * 1e3:.0f} mm station)", os.path.join(out, "fig2_section.png"))

    xy = np.asarray(geo["outlines"][k], float)
    chord = (xy[:, 0].max() - xy[:, 0].min()) * 1e3
    lay = lambda n: f"{n}-ply skin, {n * a.ply_t * 1e3:.2f} mm"
    if zone_end:
        zones = [(0, zone_end, lay(n_in), TINT_A), (zone_end, s_tip * 1e3, lay(n_base), TINT_B)]
    else:
        zones = [(0, s_tip * 1e3, lay(n_in), TINT_B)]
    n_half = int(np.sum(ribs < s_tip - 1e-6))
    n_full = 2 * n_half - int(np.any(np.abs(ribs) < 1e-9))
    fig_ribs(s, allow, ribs, bays, a.mount, fixed, zones, chord, a.spar_aft * 1e3, K,
             f"Half wing, top view: {n_half} ribs per side ({n_full} in the wing), stations in mm",
             os.path.join(out, "fig3_rib_layout.png"), lofts)

    print(f"Variant: {name}   (rib K shown = {K:g})")
    print(f"Limit load per half wing: main element {F_main:.0f} N, flaps + endplate {F_tip:.0f} N")
    print(f"Ultimate bending moment at the mount: {M_peak:.1f} N m")
    print("Rib stations [mm]: " + ", ".join(f"{r * 1e3:.0f}" for r in ribs))
    print("Bay / allowable [mm]: " + ", ".join(f"{b[2] * 1e3:.0f}/{b[3] * 1e3:.0f}" for b in bays))
    print(f"Tightest bay uses {100 * max(b[2] / b[3] for b in bays):.0f}% of its allowable spacing")
    print(f"From the saved run: skin + plate {v.m * 1e3:.0f} g per half, tip deflection {v.w:.2f} mm, "
          f"twist {v.th:.2f} deg (limit), min strength margin {v.MS:.1f}")
    print(f"Figures written to {out}")


if __name__ == "__main__":
    main()
