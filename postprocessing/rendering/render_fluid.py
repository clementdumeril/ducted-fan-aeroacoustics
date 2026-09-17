# pvpython : ANIMATION FLUIDE - camera fixe, les vortex Q tournent (lit les frames WSL)
from paraview.simple import *
import os
BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
WSL = r"\\wsl.localhost\Ubuntu\home\qltus\ventilo_sst\postProcessing\animFrames"
OUT = os.path.join(BASE, "anim_fluid"); os.makedirs(OUT, exist_ok=True)
times = [l.strip() for l in open(r"C:\Users\PC\frametimes.txt") if l.strip()]
print(len(times), "frames a rendre")
paraview.simple._DisableFirstRenderCameraReset()

rv = GetActiveViewOrCreate("RenderView")
rv.ViewSize = [1280, 960]; rv.OrientationAxesVisibility = 0
for a, v in [("UseColorPaletteForBackground", 0), ("BackgroundColorMode", "Single Color"), ("EnableRayTracing", 0)]:
    try: setattr(rv, a, v)
    except Exception: pass
rv.Background = [1, 1, 1]

for i, t in enumerate(times):
    q = XMLPolyDataReader(FileName=[os.path.join(WSL, t, "Qiso.vtp")])
    pales = XMLPolyDataReader(FileName=[os.path.join(WSL, t, "pales.vtp")])
    dq = Show(q, rv)
    ColorBy(dq, ("CELLS", "U", "Magnitude"))
    dq.SetScalarBarVisibility(rv, True)
    uLUT = GetColorTransferFunction("U"); uLUT.ApplyPreset("Turbo", True)
    uLUT.RescaleTransferFunction(0.0, 8.0)          # PLAGE FIXE -> pas de scintillement
    sb = GetScalarBar(uLUT, rv); sb.Title = "|U| (m/s)"; sb.ComponentTitle = ""
    sb.TitleColor = [0, 0, 0]; sb.LabelColor = [0, 0, 0]
    dp = Show(pales, rv)
    ColorBy(dp, ("CELLS", "p")); dp.SetScalarBarVisibility(rv, False)
    GetColorTransferFunction("p").ApplyPreset("Grayscale", True)
    # CAMERA FIXE (vue 3/4)
    rv.CameraPosition = [0.32, -0.32, 0.22]
    rv.CameraFocalPoint = [0.0, 0.0, 0.02]
    rv.CameraViewUp = [0, 0, 1]
    Render()
    SaveScreenshot(os.path.join(OUT, f"frame_{i:04d}.png"), rv, ImageResolution=[1280, 960])
    Delete(dp); Delete(dq); Delete(pales); Delete(q)
    if i % 20 == 0: print(f"  frame {i}/{len(times)} (t={t})")
print("OK fluide:", len(times), "frames")
