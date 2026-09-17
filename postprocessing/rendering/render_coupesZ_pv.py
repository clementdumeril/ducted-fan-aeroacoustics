# pvpython : coupes transverses multi-z colorees par |U| (rapide)
from paraview.simple import *
import os
BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
SRC = os.path.join(BASE, "wsl_data", "coupesZ")
OUT = os.path.join(BASE, "pack_visuels"); os.makedirs(OUT, exist_ok=True)
paraview.simple._DisableFirstRenderCameraReset()
planes = [("zm50", -50), ("z0", 0), ("z20", 20), ("z40", 40), ("z60", 60),
          ("z80", 80), ("z100", 100), ("z120", 120), ("z150", 150), ("z200", 200)]

rv = GetActiveViewOrCreate("RenderView")
rv.ViewSize = [1000, 1000]; rv.OrientationAxesVisibility = 0
for a, v in [("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"), ("EnableRayTracing", 0)]:
    try: setattr(rv, a, v)
    except Exception: pass
rv.Background = [1, 1, 1]; rv.CameraParallelProjection = 1

for name, zmm in planes:
    f = os.path.join(SRC, name + ".vtp")
    if not os.path.exists(f): print("absent", name); continue
    s = XMLPolyDataReader(FileName=[f])
    d = Show(s, rv)
    ColorBy(d, ("CELLS", "U", "Magnitude"))
    d.SetScalarBarVisibility(rv, True)
    uLUT = GetColorTransferFunction("U"); uLUT.ApplyPreset("Turbo", True); uLUT.RescaleTransferFunction(0, 8)
    sb = GetScalarBar(uLUT, rv); sb.Title = "|U| (m/s)"; sb.ComponentTitle = ""
    sb.TitleColor = [0, 0, 0]; sb.LabelColor = [0, 0, 0]
    # vue de face (le long de +z)
    rv.CameraPosition = [0, 0, zmm/1000.0 + 0.5]; rv.CameraFocalPoint = [0, 0, zmm/1000.0]; rv.CameraViewUp = [0, 1, 0]
    rv.CameraParallelScale = 0.16
    Render()
    SaveScreenshot(os.path.join(OUT, f"coupeZ_{name}.png"), rv, ImageResolution=[1000, 1000])
    Delete(d); Delete(s)
    print("  coupe", name)
print("OK coupes-z pv")
