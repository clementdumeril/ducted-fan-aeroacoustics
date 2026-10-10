# -*- coding: utf-8 -*-
"""Blade loads versus time and established-regime performance (thrust, torque, power, flow rate).

Loads come from the OpenFOAM `forces` function object (data/forces_data/). The flow rate is
integrated on the downstream cut z = 60 mm, which needs the exported VTP slices
(pass the folder as first argument, or set FLOW_VTK; default: data/flow_vtk).
"""
import os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
FORCES = REPO / "data" / "forces_data"
VTK_DIR = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FLOW_VTK", REPO / "data" / "flow_vtk"))
OUT = REPO / "figures" / "acoustics"
PAPER = REPO / "paper" / "figures"


def save(name):
    """PNG for the README, vector PDF for the paper."""
    plt.savefig(OUT / f"{name}.png", dpi=130)
    PAPER.mkdir(parents=True, exist_ok=True)
    plt.savefig(PAPER / f"{name}.pdf")

OMEGA = 158.0  # rad/s


def load():
    t, Fz, Mz = [], [], []
    for tag in ["0.15", "0.19"]:
        f = np.loadtxt(FORCES / f"force_{tag}.dat")
        m = np.loadtxt(FORCES / f"moment_{tag}.dat")
        t.append(f[:, 0]); Fz.append(f[:, 3]); Mz.append(m[:, 3])
    t, Fz, Mz = map(np.concatenate, (t, Fz, Mz))
    o = np.argsort(t)
    return t[o], Fz[o], Mz[o]


t, Fz, Mz = load()

# --- loads versus time ---
fig, ax1 = plt.subplots(figsize=(11, 4.8))
ax1.plot(t * 1000, -Fz * 1000, color="#185FA5", lw=1.0)
ax1.set_xlabel("Time (ms)"); ax1.set_ylabel("Thrust -Fz (mN)", color="#185FA5")
ax1.tick_params(axis="y", labelcolor="#185FA5"); ax1.grid(alpha=.3)
ax2 = ax1.twinx()
ax2.plot(t * 1000, -Mz * 1000, color="#D85A30", lw=1.0)
ax2.set_ylabel("Torque -Mz (mN·m)", color="#D85A30"); ax2.tick_params(axis="y", labelcolor="#D85A30")
ax1.set_title("Blade loads: start-up transient, then an established periodic regime")
plt.tight_layout(); save("blade_loads_time"); plt.close()

# --- established regime: mean over the last 50 ms ---
mask = t > t.max() - 0.05
Fz_m, Mz_m = Fz[mask].mean(), Mz[mask].mean()
P = abs(Mz_m) * OMEGA
print(f"Thrust            = {-Fz_m*1000:.1f} mN")
print(f"Torque            = {-Mz_m*1000:.3f} mN.m")
print(f"Aero power C*w    = {P*1000:.1f} mW")

# --- volume flow rate through the downstream cut (z = 60 mm) ---
aval = VTK_DIR / "avalZ.vtp"
if aval.exists():
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    r = vtk.vtkXMLPolyDataReader(); r.SetFileName(str(aval)); r.Update()
    tf = vtk.vtkTriangleFilter(); tf.SetInputData(r.GetOutput()); tf.Update()
    pd = tf.GetOutput()
    pts = vtk_to_numpy(pd.GetPoints().GetData())
    tri = vtk_to_numpy(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
    U = vtk_to_numpy(pd.GetCellData().GetArray("U"))
    area = 0.5 * np.linalg.norm(np.cross(pts[tri[:, 1]] - pts[tri[:, 0]], pts[tri[:, 2]] - pts[tri[:, 0]]), axis=1)
    Q = np.sum(area * U[:, 2])
    print(f"Net flow rate     = {Q*1000:.1f} L/s ({Q*3600:.0f} m3/h)")
else:
    print(f"(flow rate skipped: {aval} not found)")
