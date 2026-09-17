# -*- coding: utf-8 -*-
"""Rendu des coupes d'ecoulement (slice meridienne + transverse aval) depuis les VTP OpenFOAM."""
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.tri as mtri

BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
VTK = BASE + r"\flow_vtk"

def read(path):
    r = vtk.vtkXMLPolyDataReader(); r.SetFileName(path); r.Update()
    tf = vtk.vtkTriangleFilter(); tf.SetInputData(r.GetOutput()); tf.Update()
    c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(tf.GetOutput()); c2p.Update()
    pd = c2p.GetOutput()
    pts = vtk_to_numpy(pd.GetPoints().GetData())
    polys = vtk_to_numpy(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
    U = vtk_to_numpy(pd.GetPointData().GetArray("U"))
    p = vtk_to_numpy(pd.GetPointData().GetArray("p"))
    return pts, polys, U, p

# ========== COUPE MERIDIENNE (y=0) : axe z horizontal (flux ->), x vertical (rayon) ==========
pts, tri, U, p = read(VTK + r"\meridien.vtp")
z = pts[:, 2] * 1000   # axial (mm)
x = pts[:, 0] * 1000   # radial (mm)
Umag = np.linalg.norm(U, axis=1)
T = mtri.Triangulation(z, x, tri)
# masque les triangles trop grands (trous solides : carter/moteur/pales)
zt = z[tri]; xt = x[tri]
size = np.maximum(zt.max(1)-zt.min(1), xt.max(1)-xt.min(1))
T.set_mask(size > 12)   # mm

# --- vitesse ---
fig, ax = plt.subplots(figsize=(11, 5.2))
lv = np.linspace(0, max(Umag.max()*0.9, 6), 25)
c = ax.tricontourf(T, Umag, levels=lv, cmap="turbo", extend="max")
plt.colorbar(c, ax=ax, label="|U| (m/s)", shrink=0.85)
ax.axvline(0, color="white", lw=0.6, ls="--", alpha=0.5)
ax.set_xlabel("z — axe du ventilateur (mm)   [amont < 0 < aval]")
ax.set_ylabel("rayon x (mm)")
ax.set_title("Champ de vitesse — coupe meridienne : le JET et le SILLAGE derriere les pales")
ax.set_aspect("equal"); ax.set_xlim(-120, 250); ax.set_ylim(-110, 110)
plt.tight_layout(); plt.savefig(BASE + r"\flow_1_vitesse_meridienne.png", dpi=130); plt.close()

# --- vitesse + vecteurs (sillage) ---
fig, ax = plt.subplots(figsize=(11, 5.2))
c = ax.tricontourf(T, Umag, levels=lv, cmap="turbo", extend="max", alpha=0.9)
plt.colorbar(c, ax=ax, label="|U| (m/s)", shrink=0.85)
# vecteurs sur grille reguliere (interpolation robuste griddata)
from scipy.interpolate import griddata
zi, xi = np.meshgrid(np.linspace(-50, 230, 34), np.linspace(-95, 95, 22))
pts2d = np.column_stack([z, x])
uu = griddata(pts2d, U[:, 2], (zi, xi), method="linear")
vv = griddata(pts2d, U[:, 0], (zi, xi), method="linear")
ax.quiver(zi, xi, uu, vv, color="white", scale=130, width=0.0022, alpha=0.85)
ax.set_xlabel("z — axe (mm)"); ax.set_ylabel("rayon x (mm)")
ax.set_title("Ecoulement derriere les pales : jet axial + recirculation dans le sillage du moyeu")
ax.set_aspect("equal"); ax.set_xlim(-120, 250); ax.set_ylim(-110, 110)
plt.tight_layout(); plt.savefig(BASE + r"\flow_2_sillage_vecteurs.png", dpi=130); plt.close()

# --- pression ---
fig, ax = plt.subplots(figsize=(11, 5.2))
pPa = p * 1.2   # cinematique -> Pa
vmax = np.percentile(np.abs(pPa), 98)
c = ax.tricontourf(T, pPa, levels=np.linspace(-vmax, vmax, 25), cmap="RdBu_r", extend="both")
plt.colorbar(c, ax=ax, label="pression (Pa)", shrink=0.85)
ax.set_xlabel("z — axe (mm)"); ax.set_ylabel("rayon x (mm)")
ax.set_title("Champ de pression — coupe meridienne (depression a l'aspiration, surpression au refoulement)")
ax.set_aspect("equal"); ax.set_xlim(-120, 250); ax.set_ylim(-110, 110)
plt.tight_layout(); plt.savefig(BASE + r"\flow_3_pression_meridienne.png", dpi=130); plt.close()

# ========== COUPE TRANSVERSE AVAL (z=60mm) : le sillage vu de face ==========
pts2, tri2, U2, p2 = read(VTK + r"\avalZ.vtp")
xx = pts2[:, 0] * 1000; yy = pts2[:, 1] * 1000
Umag2 = np.linalg.norm(U2, axis=1)
T2 = mtri.Triangulation(xx, yy, tri2)
xt2 = xx[tri2]; yt2 = yy[tri2]
size2 = np.maximum(xt2.max(1)-xt2.min(1), yt2.max(1)-yt2.min(1))
T2.set_mask(size2 > 12)
fig, ax = plt.subplots(figsize=(6.6, 5.6))
c = ax.tricontourf(T2, Umag2, levels=np.linspace(0, max(Umag2.max()*0.9, 6), 25), cmap="turbo", extend="max")
plt.colorbar(c, ax=ax, label="|U| (m/s)")
ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
ax.set_title("Sillage vu de face (coupe a 60 mm en aval)\nempreinte des 5 pales")
ax.set_aspect("equal")
plt.tight_layout(); plt.savefig(BASE + r"\flow_4_sillage_transverse.png", dpi=130); plt.close()

print("OK. Figures: flow_1_vitesse_meridienne, flow_2_sillage_vecteurs, flow_3_pression_meridienne, flow_4_sillage_transverse")
print(f"|U| max meridien = {Umag.max():.1f} m/s,  jet ~ {np.percentile(Umag[z>50],90):.1f} m/s")
