# pvpython : animation CAMERA ORBITALE autour des vortex Q + pales (aucun re-calcul)
from paraview.simple import *
import os, math
BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
VTK = os.path.join(BASE, "flow_vtk")
OUT = os.path.join(BASE, "anim_cam")
os.makedirs(OUT, exist_ok=True)
paraview.simple._DisableFirstRenderCameraReset()

q = XMLPolyDataReader(FileName=[os.path.join(VTK, "Qiso.vtp")])
pales = XMLPolyDataReader(FileName=[os.path.join(VTK, "pales.vtp")])
rv = GetActiveViewOrCreate("RenderView")
rv.ViewSize = [1280, 960]; rv.OrientationAxesVisibility = 0
for a, v in [("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"), ("EnableRayTracing", 0)]:
    try: setattr(rv, a, v)
    except Exception: pass
rv.Background = [1, 1, 1]

# vortex colores par |U|
dq = Show(q, rv)
ColorBy(dq, ("CELLS", "U", "Magnitude"))
dq.RescaleTransferFunctionToDataRange(True, False)
dq.SetScalarBarVisibility(rv, True)
uLUT = GetColorTransferFunction("U"); uLUT.ApplyPreset("Turbo", True)
sb = GetScalarBar(uLUT, rv); sb.Title = "|U| (m/s)"; sb.ComponentTitle = ""
sb.TitleColor = [0, 0, 0]; sb.LabelColor = [0, 0, 0]
# pales colorees par p (geometrie fixe -> pas de scintillement), barre masquee
dp = Show(pales, rv)
ColorBy(dp, ("CELLS", "p"))
dp.SetScalarBarVisibility(rv, False)
GetColorTransferFunction("p").ApplyPreset("Grayscale", True)

# orbite camera autour de l'axe z
NFR = 120; R = 0.34; H = 0.15; FOC = [0.0, 0.0, 0.02]
for i in range(NFR):
    th = 2 * math.pi * i / NFR
    rv.CameraPosition = [R * math.cos(th), R * math.sin(th), H]
    rv.CameraFocalPoint = FOC
    rv.CameraViewUp = [0, 0, 1]
    Render()
    SaveScreenshot(os.path.join(OUT, f"frame_{i:04d}.png"), rv, ImageResolution=[1280, 960])
    if i % 20 == 0: print(f"frame {i}/{NFR}")
print("OK frames:", NFR)
