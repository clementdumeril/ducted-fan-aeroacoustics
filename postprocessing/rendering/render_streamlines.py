# pvpython : STREAMLINES 3D (lignes de courant) sur le volume, colorees par vitesse
from paraview.simple import *
import os
BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
VOL = os.path.join(BASE, "wsl_data", "vol", "internal.vtu")
PALES = os.path.join(BASE, "flow_vtk", "pales.vtp")
OUT = os.path.join(BASE, "pack_visuels"); os.makedirs(OUT, exist_ok=True)
paraview.simple._DisableFirstRenderCameraReset()

vol = XMLUnstructuredGridReader(FileName=[VOL])
c2p = CellDatatoPointData(Input=vol)   # U en donnees de points pour le traceur

st = StreamTracer(Input=c2p, SeedType="Line")
st.Vectors = ["POINTS", "U"]
st.SeedType.Point1 = [-0.09, 0.0, -0.04]
st.SeedType.Point2 = [0.09, 0.0, -0.04]
st.SeedType.Resolution = 70
st.IntegrationDirection = "BOTH"
st.MaximumStreamlineLength = 0.8
tube = Tube(Input=st); tube.Radius = 0.0009

rv = GetActiveViewOrCreate("RenderView")
rv.ViewSize = [1400, 1000]; rv.OrientationAxesVisibility = 0
for a, v in [("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"), ("EnableRayTracing", 0)]:
    try: setattr(rv, a, v)
    except Exception: pass
rv.Background = [1, 1, 1]

dt = Show(tube, rv)
ColorBy(dt, ("POINTS", "U", "Magnitude"))
dt.RescaleTransferFunctionToDataRange(True, False)
dt.SetScalarBarVisibility(rv, True)
uLUT = GetColorTransferFunction("U"); uLUT.ApplyPreset("Turbo", True); uLUT.RescaleTransferFunction(0, 8)
sb = GetScalarBar(uLUT, rv); sb.Title = "|U| (m/s)"; sb.ComponentTitle = ""
sb.TitleColor = [0, 0, 0]; sb.LabelColor = [0, 0, 0]
# pales en gris pour le contexte
pales = XMLPolyDataReader(FileName=[PALES]); dp = Show(pales, rv)
ColorBy(dp, ("CELLS", "p")); dp.SetScalarBarVisibility(rv, False)
GetColorTransferFunction("p").ApplyPreset("Grayscale", True)

# 3 vues
views = [("3q", [0.34, -0.34, 0.24], [0, 0, 0.05]),
         ("profil", [0.0, -0.5, 0.06], [0, 0, 0.06]),
         ("aval", [0.45, 0.0, 0.10], [0, 0, 0.05])]
for name, pos, foc in views:
    rv.CameraPosition = pos; rv.CameraFocalPoint = foc; rv.CameraViewUp = [0, 0, 1]
    Render(); SaveScreenshot(os.path.join(OUT, f"streamlines_{name}.png"), rv, ImageResolution=[1400, 1000])
print("OK streamlines 3 vues")
