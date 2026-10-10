# -*- coding: utf-8 -*-
"""2-D flow figures (meridional cut y = 0 and downstream cut z = 60 mm) from the OpenFOAM VTP slices.

The slices are exported from the last time step with ParaView (meridien.vtp, avalZ.vtp).
Pass their folder as first argument or set FLOW_VTK (default: data/flow_vtk).
"""
import os, sys
from pathlib import Path
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.interpolate import griddata

REPO = Path(__file__).resolve().parents[2]
VTK_DIR = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FLOW_VTK", REPO / "data" / "flow_vtk"))
OUT = REPO / "figures" / "flow"
PAPER = REPO / "paper" / "figures"


def save(name):
    """PNG for the README; 200 dpi PNG for the paper (filled contours on ~10^5 triangles are too heavy as PDF)."""
    plt.savefig(OUT / f"{name}.png", dpi=130)
    PAPER.mkdir(parents=True, exist_ok=True)
    plt.savefig(PAPER / f"{name}.png", dpi=200)

if not (VTK_DIR / "meridien.vtp").exists() or not (VTK_DIR / "avalZ.vtp").exists():
    sys.exit(f"meridien.vtp and avalZ.vtp not found in {VTK_DIR}.\n"
             "These ParaView slices of the last time step are not published; export them from the case\n"
             "and pass their folder: python render_flow_2d.py <folder>")


def read(path):
    r = vtk.vtkXMLPolyDataReader(); r.SetFileName(str(path)); r.Update()
    tf = vtk.vtkTriangleFilter(); tf.SetInputData(r.GetOutput()); tf.Update()
    c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(tf.GetOutput()); c2p.Update()
    pd = c2p.GetOutput()
    pts = vtk_to_numpy(pd.GetPoints().GetData())
    polys = vtk_to_numpy(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
    U = vtk_to_numpy(pd.GetPointData().GetArray("U"))
    p = vtk_to_numpy(pd.GetPointData().GetArray("p"))
    return pts, polys, U, p


def masked_tri(a, b, tri, max_mm=12):
    """Triangulation with the large triangles (solid holes: casing, motor, blades) masked out."""
    T = mtri.Triangulation(a, b, tri)
    at, bt = a[tri], b[tri]
    T.set_mask(np.maximum(at.max(1) - at.min(1), bt.max(1) - bt.min(1)) > max_mm)
    return T


def finish(ax, title, name):
    ax.set_xlabel("z, fan axis (mm)   [upstream < 0 < downstream]"); ax.set_ylabel("radius x (mm)")
    ax.set_title(title)
    ax.set_aspect("equal"); ax.set_xlim(-120, 250); ax.set_ylim(-110, 110)
    plt.tight_layout(); save(name.removesuffix(".png")); plt.close()


# ---------- meridional cut (y = 0): axial z horizontal, radius x vertical ----------
pts, tri, U, p = read(VTK_DIR / "meridien.vtp")
z = pts[:, 2] * 1000; x = pts[:, 0] * 1000
Umag = np.linalg.norm(U, axis=1)
T = masked_tri(z, x, tri)
lv = np.linspace(0, 7.0, 29)   # fixed 0-7 m/s scale, round ticks

fig, ax = plt.subplots(figsize=(11, 5.2))
c = ax.tricontourf(T, Umag, levels=lv, cmap="turbo", extend="max")
plt.colorbar(c, ax=ax, label="|U| (m/s)", shrink=0.85, ticks=range(0, 8))
ax.axvline(0, color="white", lw=0.6, ls="--", alpha=0.5)
finish(ax, "Velocity magnitude, meridional cut: annular jet and hub wake", "meridional_velocity.png")

fig, ax = plt.subplots(figsize=(11, 5.2))
c = ax.tricontourf(T, Umag, levels=lv, cmap="turbo", extend="max", alpha=0.9)
plt.colorbar(c, ax=ax, label="|U| (m/s)", shrink=0.85, ticks=range(0, 8))
zi, xi = np.meshgrid(np.linspace(-50, 230, 34), np.linspace(-95, 95, 22))
uu = griddata(np.column_stack([z, x]), U[:, 2], (zi, xi), method="linear")
vv = griddata(np.column_stack([z, x]), U[:, 0], (zi, xi), method="linear")
ax.quiver(zi, xi, uu, vv, color="white", scale=130, width=0.0022, alpha=0.85)
finish(ax, "In-plane velocity vectors: axial jet and recirculation behind the hub", "wake_vectors.png")

fig, ax = plt.subplots(figsize=(11, 5.2))
pPa = p * 1.2   # kinematic pressure -> Pa
vmax = np.percentile(np.abs(pPa), 98)
c = ax.tricontourf(T, pPa, levels=np.linspace(-vmax, vmax, 25), cmap="RdBu_r", extend="both")
plt.colorbar(c, ax=ax, label="p (Pa)", shrink=0.85)
finish(ax, "Static pressure, meridional cut: suction upstream of the rotor, overpressure downstream",
       "meridional_pressure.png")

# ---------- downstream cut (z = 60 mm): the wake seen from the front ----------
pts2, tri2, U2, _ = read(VTK_DIR / "avalZ.vtp")
xx = pts2[:, 0] * 1000; yy = pts2[:, 1] * 1000
Umag2 = np.linalg.norm(U2, axis=1)
T2 = masked_tri(xx, yy, tri2)
fig, ax = plt.subplots(figsize=(6.6, 5.6))
c = ax.tricontourf(T2, Umag2, levels=np.linspace(0, 6.0, 25), cmap="turbo", extend="max")
plt.colorbar(c, ax=ax, label="|U| (m/s)", ticks=range(0, 7))
ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
ax.set_title("Wake 60 mm downstream of the rotor")
ax.set_aspect("equal"); ax.set_xlim(-110, 110); ax.set_ylim(-110, 110)
plt.tight_layout(); save("wake_cut_z60"); plt.close()

print(f"|U| max (meridional) = {Umag.max():.1f} m/s, jet p90 (z > 50 mm) = {np.percentile(Umag[z > 50], 90):.1f} m/s")
print(f"Figures written to {OUT}")
