"""Minimal semiconducting hidden-altermagnet bilayer (no SOC).

Run: python semiconductor_bilayer.py
Outputs figures in assets/ and plain-text numerical data in data/.
All energies are dimensionless model units; lattice constants are set to one.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ASSETS = ROOT / "assets"
DATA = ROOT / "data"
TX, TY, TAB, DELTA = 0.20, -0.05, 0.25, 0.60
TPERP, V = 0.08, 0.06  # onsite energies: upper +V, lower -V
MU = 0.16             # illustrative electron-doped constant-energy cut
NGRID, NSTEPS = 301, 100
P = np.eye(4)[[3, 2, 1, 0]]  # Mu<->Nd, Nu<->Md
MXY = np.eye(4)[[1, 0, 3, 2]]
BLUE, RED = "#1764ab", "#cf4141"


def terms(kx, ky):
    cx, cy = np.cos(kx), np.cos(ky)
    e0 = -(TX + TY) * (cx + cy)
    d = -(TX - TY) * (cx - cy)
    g = -4 * TAB * np.cos(kx / 2) * np.cos(ky / 2)
    return e0, d, g


def monolayer(kx, ky, spin):
    """2x2 block in (M,N). spin=+1 (up), -1 (down)."""
    e0, d, g = terms(kx, ky)
    h = np.zeros(np.shape(e0) + (2, 2))
    h[..., 0, 0] = e0 + d - spin * DELTA
    h[..., 1, 1] = e0 - d + spin * DELTA
    h[..., 0, 1] = h[..., 1, 0] = g
    return h


def bilayer(kx, ky, spin, potential):
    """4x4 block in (Mu,Nu,Md,Nd); lower layer is inversion partner."""
    e0, d, g = terms(kx, ky)
    h = np.zeros(np.shape(e0) + (4, 4))
    h[..., 0, 0] = e0 + d - spin * DELTA + potential
    h[..., 1, 1] = e0 - d + spin * DELTA + potential
    h[..., 2, 2] = e0 - d - spin * DELTA - potential
    h[..., 3, 3] = e0 + d + spin * DELTA - potential
    h[..., 0, 1] = h[..., 1, 0] = g
    h[..., 2, 3] = h[..., 3, 2] = g
    h[..., 0, 2] = h[..., 2, 0] = -TPERP
    h[..., 1, 3] = h[..., 3, 1] = -TPERP
    return h


def path():
    nodes = [("Y", (0, np.pi)), (r"$\Gamma$", (0, 0)),
             ("X", (np.pi, 0)), ("M", (np.pi, np.pi)),
             (r"$\Gamma$", (0, 0))]
    segments = []
    for (_, a), (_, b) in zip(nodes[:-1], nodes[1:]):
        f = np.linspace(0, 1, NSTEPS, endpoint=False)[:, None]
        segments.append(np.asarray(a) + f * (np.asarray(b) - np.asarray(a)))
    k = np.vstack(segments + [np.asarray(nodes[-1][1])[None, :]])
    dist = np.r_[0, np.cumsum(np.linalg.norm(np.diff(k, axis=0), axis=1))]
    ticks = np.arange(0, len(k), NSTEPS)
    return k, dist, ticks, [n for n, _ in nodes]


def band_edges(up, down):
    e = np.sort(np.concatenate([up, down], axis=-1), axis=-1)
    half = e.shape[-1] // 2
    vmax = float(e[..., half - 1].max())
    cmin = float(e[..., half].min())
    return vmax, cmin, cmin - vmax


def save_txt(name, columns, array, note=""):
    np.savetxt(DATA / name, array, fmt="%.10e",
               header=f"{columns}\n{note}")


def axes_bz(ax):
    ax.set_aspect("equal")
    ax.set(xlim=(-np.pi, np.pi), ylim=(-np.pi, np.pi),
           xlabel=r"$k_x$", ylabel=r"$k_y$")
    ax.set_xticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])
    ax.set_yticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])


def draw_bands(ax, distance, ticks, labels, up, down, title, edges):
    for energies, color, label in [(up, BLUE, "up"), (down, RED, "down")]:
        for i in range(energies.shape[-1]):
            ax.plot(distance, energies[:, i], color=color, lw=1.15,
                    label=label if i == 0 else None)
    for t in ticks[1:-1]:
        ax.axvline(distance[t], color="0.88", lw=.7)
    ax.axhline(edges[0], color="0.35", ls=":", lw=.7)
    ax.axhline(edges[1], color="0.35", ls=":", lw=.7)
    ax.axhline(0, color="0.65", lw=.5)
    ax.set_xticks(distance[ticks], labels)
    ax.set(xlim=(distance[0], distance[-1]), ylim=(-1.75, 1.75),
           ylabel="Energy (model units)", title=title)
    ax.legend(frameon=False, fontsize=8, loc="upper right")


def main():
    ASSETS.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    k = np.linspace(-np.pi, np.pi, NGRID)
    kx, ky = np.meshgrid(k, k, indexing="xy")
    pk, distance, ticks, labels = path()
    grids, paths, edges = {}, {}, {}
    for key, hfun, potential in [
        ("mono", monolayer, None), ("zero", bilayer, 0.0),
        ("field", bilayer, V),
    ]:
        args_grid = (kx, ky)
        args_path = (pk[:, 0], pk[:, 1])
        grids[key] = tuple(np.linalg.eigvalsh(hfun(*args_grid, s) if potential is None
                            else hfun(*args_grid, s, potential)) for s in (+1, -1))
        paths[key] = tuple(np.linalg.eigvalsh(hfun(*args_path, s) if potential is None
                            else hfun(*args_path, s, potential)) for s in (+1, -1))
        edges[key] = band_edges(*grids[key])
        if edges[key][2] <= 0:
            raise RuntimeError(f"{key}: indirect gap closed; change parameters")
        up, down = grids[key]
        flat = np.column_stack((kx.ravel(), ky.ravel(),
                                up.reshape(-1, up.shape[-1]),
                                down.reshape(-1, down.shape[-1])))
        columns = "kx ky " + " ".join(f"up_{i+1}" for i in range(up.shape[-1]))
        columns += " " + " ".join(f"down_{i+1}" for i in range(down.shape[-1]))
        save_txt(f"{key}_grid.txt", columns, flat,
                 f"Energy in model units; {NGRID}x{NGRID}; ky outer, kx inner")
        pu, pd = paths[key]
        save_txt(f"{key}_bands.txt", "distance kx ky " +
                 " ".join(f"up_{i+1}" for i in range(pu.shape[-1])) + " " +
                 " ".join(f"down_{i+1}" for i in range(pd.shape[-1])),
                 np.column_stack((distance, pk, pu, pd)),
                 "Path: Y-Gamma-X-M-Gamma; energy in model units")

    # Algebraic symmetry checks at a generic momentum, not only along a path.
    qx, qy = .37, 1.11
    pt_zero = float(np.max(np.abs(P @ bilayer(qx, qy, +1, 0) @ P.T
                                  - bilayer(qx, qy, -1, 0))))
    pt_field = float(np.max(np.abs(P @ bilayer(qx, qy, +1, V) @ P.T
                                   - bilayer(qx, qy, -1, V))))
    mirror = float(np.max(np.abs(MXY @ bilayer(qx, qy, +1, V) @ MXY.T
                                 - bilayer(qy, qx, -1, V))))
    zero_split = float(np.max(np.abs(grids["zero"][0] - grids["zero"][1])))
    if max(pt_zero, mirror, zero_split) > 1e-10 or abs(pt_field - 2*V) > 1e-10:
        raise RuntimeError("Symmetry test failed")

    # Figure 1: lattice and isolated-layer bands.
    fig, axs = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)
    ax = axs[0]
    for x in range(-1, 2):
        for y in range(-1, 2):
            ax.scatter([x], [y+.5], c=BLUE, s=90, marker="o")
            ax.scatter([x+.5], [y], c=RED, s=90, marker="s")
    ax.annotate("M ↑", (0, .5), (.10, .69), color=BLUE)
    ax.annotate("N ↓", (.5, 0), (.62, -.19), color=RED)
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=".5", lw=1)
    ax.annotate(r"x bonds: M $t_x$, N $t_y$", (.1, -1.2), fontsize=10)
    ax.annotate(r"y bonds exchanged by $m_{xy}$", (.1, -1.4), fontsize=10)
    ax.set(xlim=(-.6, 1.55), ylim=(-1.5, 1.55), xlabel="x / a", ylabel="y / a",
           title="(a) One magnetic layer")
    ax.set_aspect("equal")
    draw_bands(axs[1], distance, ticks, labels, *paths["mono"],
               "(b) Semiconducting monolayer", edges["mono"])
    fig.savefig(ASSETS / "fig1_monolayer.png", dpi=230)
    plt.close(fig)

    # Figure 2: stacking and electric-field symmetry breaking.
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.2), constrained_layout=True)
    ax = axs[0]
    for z, suffix in [(.6, "u"), (.4, "d")]:
        ax.scatter([0, .5], [z, z], c=[BLUE, RED], s=180, marker="o", zorder=3)
        ax.text(0, z+.023, f"M{suffix} ↑", ha="center", fontsize=9)
        ax.text(.5, z+.023, f"N{suffix} ↓", ha="center", fontsize=9)
    ax.plot([0, 0], [.4, .6], color="0.4", lw=2)
    ax.plot([.5, .5], [.4, .6], color="0.4", lw=2)
    ax.scatter([.25], [.5], marker="x", c="k", s=45)
    ax.annotate("P center", (.27, .47), fontsize=9)
    ax.annotate(r"$+V$", (.77, .6), color="darkgreen")
    ax.annotate(r"$-V$", (.77, .4), color="darkgreen")
    ax.arrow(.75, .44, 0, .11, width=.004, color="darkgreen")
    ax.set(xlim=(-.18, .95), ylim=(.34, .69), xlabel="in-plane position (schematic)",
           ylabel="z / c", title="(a) Inversion-partner layers")
    draw_bands(axs[1], distance, ticks, labels, *paths["zero"],
               "(b) No field: PT degenerate", edges["zero"])
    draw_bands(axs[2], distance, ticks, labels, *paths["field"],
               "(c) Field: spin splitting", edges["field"])
    fig.savefig(ASSETS / "fig2_bilayer_bands.png", dpi=230)
    plt.close(fig)

    # Figure 3: electron-doped isoenergy curves. An insulator has no Fermi
    # surface at neutrality; this is deliberately a cut above the CBM.
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.5), constrained_layout=True)
    for ax, key, title in zip(axs, ("zero", "field"),
                              ("(a) No field", "(b) With field")):
        up, down = grids[key]
        for energies, color, spin_name in [(up, "#3c5266" if key == "zero" else BLUE,
                                            "up = down" if key == "zero" else "up"),
                                           (down, "#3c5266" if key == "zero" else RED,
                                            None if key == "zero" else "down")]:
            for i in range(energies.shape[-1]):
                values = energies[..., i]
                if values.min() < MU < values.max():
                    ax.contour(kx, ky, values, levels=[MU], colors=[color],
                               linewidths=1.35)
            if spin_name is not None:
                ax.plot([], [], color=color, lw=1.4, label=spin_name)
        axes_bz(ax)
        ax.set_title(title + rf": $E={MU:.2f}$")
        ax.legend(frameon=False, fontsize=8)
    fig.savefig(ASSETS / "fig3_constant_energy.png", dpi=230)
    plt.close(fig)

    # Figure 4: layer-resolved spin spectral weight near the same energy.
    fig, axs = plt.subplots(2, 2, figsize=(9.7, 8.4), constrained_layout=True)
    cmap = plt.colormaps["RdBu_r"].copy()
    cmap.set_bad("white")
    for row, (key, potential) in enumerate((("zero", 0.0), ("field", V))):
        spin_data = []
        for spin in (+1, -1):
            e, vec = np.linalg.eigh(bilayer(kx, ky, spin, potential))
            spin_data.append((e, vec))
        for col, (layer, orbitals) in enumerate((("upper", slice(0, 2)),
                                                  ("lower", slice(2, 4)))):
            a = []
            for energy, vectors in spin_data:
                weight = np.sum(np.abs(vectors[..., orbitals, :])**2, axis=-2)
                lorentz = .02 / np.pi / ((MU - energy)**2 + .02**2)
                a.append(np.sum(weight * lorentz, axis=-1))
            pol = (a[0] - a[1]) / (a[0] + a[1] + 1e-15)
            intensity = a[0] + a[1]
            save_txt(f"{key}_{layer}_spectral.txt",
                     "kx ky A_up A_down polarization total_intensity",
                     np.column_stack((kx.ravel(), ky.ravel(), a[0].ravel(),
                                      a[1].ravel(), pol.ravel(), intensity.ravel())),
                     f"Lorentzian width eta=0.02; E={MU}; {NGRID}x{NGRID}; ky outer")
            masked = np.ma.masked_where(intensity < .08 * intensity.max(), pol)
            ax = axs[row, col]
            im = ax.imshow(masked, origin="lower", extent=(-np.pi, np.pi)*2,
                           cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
            axes_bz(ax)
            ax.set_title(f"({chr(97+2*row+col)}) {layer}, V={potential:.2f}")
    fig.colorbar(im, ax=axs, shrink=.84, label=r"$(A_{\uparrow}-A_{\downarrow})/(A_{\uparrow}+A_{\downarrow})$")
    fig.savefig(ASSETS / "fig4_layer_spectral.png", dpi=220)
    plt.close(fig)

    with (DATA / "summary.txt").open("w", encoding="utf-8") as f:
        f.write("Energies in model units; coordinates in lattice units.\n")
        f.write(f"tx={TX}, ty={TY}, tab={TAB}, delta={DELTA}, tperp={TPERP}, V={V}, mu={MU}\n")
        for key in ("mono", "zero", "field"):
            vmax, cmin, gap = edges[key]
            f.write(f"{key}: VBM={vmax:.10f}, CBM={cmin:.10f}, indirect_gap={gap:.10f}\n")
        f.write(f"PT matrix error (V=0): {pt_zero:.3e}\n")
        f.write(f"PT matrix error (V={V}): {pt_field:.3e} (expected 2V={2*V:.3e})\n")
        f.write(f"Zero-field spin-band difference: {zero_split:.3e}\n")
        f.write(f"Field spin-mirror matrix error: {mirror:.3e}\n")
    print((DATA / "summary.txt").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
