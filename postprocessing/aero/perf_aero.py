# -*- coding: utf-8 -*-
"""Efforts vs temps (regime etabli) + performances aero (poussee, couple, puissance, debit)."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
OMEGA = 158.0  # rad/s

def load(prefix):
    t, Fz, Mz = [], [], []
    for tag in ["0.15", "0.19"]:
        try:
            f = np.loadtxt(f"{BASE}\\forces_data\\force_{tag}.dat")
            m = np.loadtxt(f"{BASE}\\forces_data\\moment_{tag}.dat")
            t.append(f[:, 0]); Fz.append(f[:, 3]); Mz.append(m[:, 3])
        except Exception as e: print("skip", tag, e)
    return np.concatenate(t), np.concatenate(Fz), np.concatenate(Mz)

t, Fz, Mz = load("force")
order = np.argsort(t); t, Fz, Mz = t[order], Fz[order], Mz[order]

# === FIG efforts vs temps ===
fig, ax1 = plt.subplots(figsize=(11, 4.8))
ax1.plot(t*1000, -Fz*1000, color="#185FA5", lw=1.0, label="poussee -Fz (mN)")
ax1.set_xlabel("temps (ms)"); ax1.set_ylabel("poussee -Fz (mN)", color="#185FA5")
ax1.tick_params(axis="y", labelcolor="#185FA5"); ax1.grid(alpha=.3)
ax2 = ax1.twinx()
ax2.plot(t*1000, -Mz*1000, color="#D85A30", lw=1.0, label="couple -Mz (mN.m)")
ax2.set_ylabel("couple -Mz (mN.m)", color="#D85A30"); ax2.tick_params(axis="y", labelcolor="#D85A30")
ax1.axvspan(0, 5, color="gray", alpha=.15)
ax1.text(8, ax1.get_ylim()[1]*0.95, "rampe omega", fontsize=8, color="gray")
ax1.set_title("Efforts sur les pales en fonction du temps : transitoire puis regime etabli")
plt.tight_layout(); plt.savefig(BASE + r"\perf_efforts_temps.png", dpi=130); plt.close()

# === valeurs de regime (moyenne sur la fin) ===
mask = t > t.max() - 0.05
Fz_m, Mz_m = Fz[mask].mean(), Mz[mask].mean()
P = abs(Mz_m) * OMEGA
print(f"Poussee  = {-Fz_m*1000:.1f} mN")
print(f"Couple   = {-Mz_m*1000:.3f} mN.m")
print(f"Puissance aero = C*omega = {P*1000:.1f} mW  ({P:.3f} W)")

# === debit volumique (integration U_z sur la coupe aval z=60mm) ===
r = vtk.vtkXMLPolyDataReader(); r.SetFileName(BASE + r"\flow_vtk\avalZ.vtp"); r.Update()
tf = vtk.vtkTriangleFilter(); tf.SetInputData(r.GetOutput()); tf.Update()
pd = tf.GetOutput()
pts = vtk_to_numpy(pd.GetPoints().GetData())
polys = vtk_to_numpy(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
U = vtk_to_numpy(pd.GetCellData().GetArray("U"))
# aire de chaque triangle + flux
Q = 0.0; Apos = 0.0
for i, (a, b, c) in enumerate(polys):
    v1 = pts[b]-pts[a]; v2 = pts[c]-pts[a]
    area = 0.5*np.linalg.norm(np.cross(v1, v2))
    uz = U[i, 2]
    Q += area*uz
    if uz > 0: Apos += area
print(f"Debit volumique net (z=60mm) = {Q*1000:.2f} L/s  ({Q*3600:.1f} m3/h)")
Vmean = Q/Apos if Apos > 0 else 0
print(f"Vitesse axiale moyenne (zones sortantes) = {Vmean:.2f} m/s")
print("Figure: perf_efforts_temps.png")
