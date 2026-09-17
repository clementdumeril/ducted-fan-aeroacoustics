# -*- coding: utf-8 -*-
"""Pack de coupes couleur (meridienne + aval) x plusieurs champs, depuis les VTK locaux."""
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os

BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
OUT = os.path.join(BASE, "pack_visuels"); os.makedirs(OUT, exist_ok=True)

def read(path):
    r = vtk.vtkXMLPolyDataReader(); r.SetFileName(path); r.Update()
    tf = vtk.vtkTriangleFilter(); tf.SetInputData(r.GetOutput()); tf.Update()
    c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(tf.GetOutput()); c2p.Update()
    pd = c2p.GetOutput()
    pts = vtk_to_numpy(pd.GetPoints().GetData())
    tri = vtk_to_numpy(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:4]
    U = vtk_to_numpy(pd.GetPointData().GetArray("U"))
    p = vtk_to_numpy(pd.GetPointData().GetArray("p"))
    Q = vtk_to_numpy(pd.GetPointData().GetArray("Q"))
    return pts, tri, U, p, Q

def masked_tri(a, b, tri, lim=12):
    T = mtri.Triangulation(a, b, tri)
    at = a[tri]; bt = b[tri]
    T.set_mask(np.maximum(at.max(1)-at.min(1), bt.max(1)-bt.min(1)) > lim)
    return T

def plot(T, val, title, fname, cmap, label, sym=False, xl="z (mm)", yl="rayon x (mm)",
         xlim=(-120, 250), ylim=(-110, 110)):
    fig, ax = plt.subplots(figsize=(11, 5.2))
    if sym:
        v = np.percentile(np.abs(val), 97); lv = np.linspace(-v, v, 25)
    else:
        lv = np.linspace(max(0, np.percentile(val, 1)), np.percentile(val, 99), 25)
    c = ax.tricontourf(T, val, levels=lv, cmap=cmap, extend="both")
    plt.colorbar(c, ax=ax, label=label, shrink=0.85)
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title)
    ax.set_aspect("equal"); ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, fname), dpi=130); plt.close()

# ===== MERIDIENNE =====
pts, tri, U, p, Q = read(os.path.join(BASE, "flow_vtk", "meridien.vtp"))
z = pts[:, 2]*1000; x = pts[:, 0]*1000
Tm = masked_tri(z, x, tri)
Umag = np.linalg.norm(U, axis=1)
plot(Tm, Umag, "Coupe meridienne - norme de la vitesse |U|", "coupe_mer_Umag.png", "turbo", "|U| (m/s)")
plot(Tm, U[:, 2], "Coupe meridienne - vitesse AXIALE U_z (le jet)", "coupe_mer_Uaxial.png", "RdBu_r", "U_z (m/s)", sym=True)
plot(Tm, U[:, 0], "Coupe meridienne - vitesse RADIALE U_x", "coupe_mer_Uradial.png", "PuOr_r", "U_x (m/s)", sym=True)
plot(Tm, U[:, 1], "Coupe meridienne - SWIRL (vitesse tangentielle)", "coupe_mer_swirl.png", "RdBu_r", "U_theta (m/s)", sym=True)
plot(Tm, p*1.2, "Coupe meridienne - pression", "coupe_mer_pression.png", "RdBu_r", "p (Pa)", sym=True)
plot(Tm, np.clip(Q, 0, np.percentile(Q, 99.5)), "Coupe meridienne - critere Q (tourbillons)", "coupe_mer_Q.png", "hot_r", "Q (1/s2)")

# ===== AVAL (z=60mm) =====
pts2, tri2, U2, p2, Q2 = read(os.path.join(BASE, "flow_vtk", "avalZ.vtp"))
xx = pts2[:, 0]*1000; yy = pts2[:, 1]*1000
Ta = masked_tri(xx, yy, tri2)
Um2 = np.linalg.norm(U2, axis=1)
plot(Ta, Um2, "Coupe aval (60mm) - |U| (empreinte des pales)", "coupe_aval_Umag.png", "turbo", "|U| (m/s)", xl="x (mm)", yl="y (mm)", xlim=(-160, 160), ylim=(-160, 160))
plot(Ta, U2[:, 2], "Coupe aval (60mm) - vitesse axiale U_z", "coupe_aval_Uaxial.png", "RdBu_r", "U_z (m/s)", sym=True, xl="x (mm)", yl="y (mm)", xlim=(-160, 160), ylim=(-160, 160))
plot(Ta, p2*1.2, "Coupe aval (60mm) - pression", "coupe_aval_pression.png", "RdBu_r", "p (Pa)", sym=True, xl="x (mm)", yl="y (mm)", xlim=(-160, 160), ylim=(-160, 160))

print("OK -", len(os.listdir(OUT)), "coupes dans pack_visuels/")
