# pvpython : iso-Q multi-niveaux + vorticite, depuis le volume local internal.vtu
from paraview.simple import *
import os
BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
VOL = os.path.join(BASE, "wsl_data", "vol", "internal.vtu")
PALES = os.path.join(BASE, "flow_vtk", "pales.vtp")
OUT = os.path.join(BASE, "pack_visuels"); os.makedirs(OUT, exist_ok=True)
paraview.simple._DisableFirstRenderCameraReset()

vol = XMLUnstructuredGridReader(FileName=[VOL])
c2p = CellDatatoPointData(Input=vol)

rv = GetActiveViewOrCreate("RenderView")
rv.ViewSize = [1400, 1000]; rv.OrientationAxesVisibility = 0
for a, v in [("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"), ("EnableRayTracing", 0)]:
    try: setattr(rv, a, v)
    except Exception: pass
rv.Background = [1, 1, 1]

def setcam():
    rv.CameraPosition = [0.32, -0.32, 0.22]; rv.CameraFocalPoint = [0, 0, 0.03]; rv.CameraViewUp = [0, 0, 1]

# pales gris (contexte)
pales = XMLPolyDataReader(FileName=[PALES]); dpz = Show(pales, rv)
ColorBy(dpz, ("CELLS", "p")); dpz.SetScalarBarVisibility(rv, False)
GetColorTransferFunction("p").ApplyPreset("Grayscale", True)

# --- ISO-Q multi-niveaux ---
for val, tag in [(20000, "faible"), (80000, "moyen"), (250000, "fort")]:
    ct = Contour(Input=c2p); ct.ContourBy = ["POINTS", "Q"]; ct.Isosurfaces = [val]
    dc = Show(ct, rv)
    ColorBy(dc, ("POINTS", "U", "Magnitude"))
    dc.SetScalarBarVisibility(rv, True)
    uL = GetColorTransferFunction("U"); uL.ApplyPreset("Turbo", True); uL.RescaleTransferFunction(0, 8)
    sb = GetScalarBar(uL, rv); sb.Title = "|U| (m/s)"; sb.ComponentTitle = ""; sb.TitleColor = [0, 0, 0]; sb.LabelColor = [0, 0, 0]
    setcam(); Render()
    SaveScreenshot(os.path.join(OUT, f"isoQ_{tag}_Q{val}.png"), rv, ImageResolution=[1400, 1000])
    Hide(ct, rv); Delete(dc); Delete(ct)
    print("  isoQ", tag)

# --- VORTICITE : coupe meridienne (y=0) ---
sl = Slice(Input=c2p); sl.SliceType = "Plane"; sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [0, 1, 0]
ds = Show(sl, rv)
ColorBy(ds, ("POINTS", "vorticity", "Magnitude"))
ds.RescaleTransferFunctionToDataRange(True, False); ds.SetScalarBarVisibility(rv, True)
vL = GetColorTransferFunction("vorticity"); vL.ApplyPreset("Cool to Warm", True)
sb2 = GetScalarBar(vL, rv); sb2.Title = "|vorticite| (1/s)"; sb2.ComponentTitle = ""; sb2.TitleColor = [0, 0, 0]; sb2.LabelColor = [0, 0, 0]
Hide(pales, rv)
rv.CameraPosition = [0, -0.6, 0.05]; rv.CameraFocalPoint = [0, 0, 0.05]; rv.CameraViewUp = [1, 0, 0]
rv.CameraParallelProjection = 1; rv.CameraParallelScale = 0.13
Render(); SaveScreenshot(os.path.join(OUT, "vorticite_meridienne.png"), rv, ImageResolution=[1400, 900])
print("OK isoQ multi + vorticite")
