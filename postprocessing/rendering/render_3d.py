# pvpython : vortex (iso-Q) colores par vitesse  +  pales colorees par pression
from paraview.simple import *
import os
BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
VTK = os.path.join(BASE, "flow_vtk")
paraview.simple._DisableFirstRenderCameraReset()

# ---------- VUE 1 : vortex Q colores par |U| ----------
q = XMLPolyDataReader(FileName=[os.path.join(VTK, "Qiso.vtp")])
rv = GetActiveViewOrCreate("RenderView")
rv.ViewSize = [1500, 950]; rv.OrientationAxesVisibility = 0
try:
    rv.UseColorPaletteForBackground = 0
except Exception: pass
try:
    rv.BackgroundColorMode = "Single Color"
except Exception: pass
try:
    rv.EnableRayTracing = 0
except Exception: pass
rv.Background = [1, 1, 1]
dq = Show(q, rv)
ColorBy(dq, ("CELLS", "U", "Magnitude"))
dq.RescaleTransferFunctionToDataRange(True, False)
dq.SetScalarBarVisibility(rv, True)
uLUT = GetColorTransferFunction("U"); uLUT.ApplyPreset("Turbo", True)
sb = GetScalarBar(uLUT, rv); sb.Title = "|U| (m/s)"; sb.ComponentTitle = ""
sb.TitleColor = [0, 0, 0]; sb.LabelColor = [0, 0, 0]
ResetCamera(rv)
rv.CameraPosition = [0.34, -0.34, 0.24]
rv.CameraFocalPoint = [0.0, 0.0, 0.04]
rv.CameraViewUp = [0, 0, 1]
Render()
SaveScreenshot(os.path.join(BASE, "flow_5_vortex_Q.png"), rv, ImageResolution=[1500, 950])

# ---------- VUE 2 : pales colorees par pression ----------
Hide(q, rv)
pales = XMLPolyDataReader(FileName=[os.path.join(VTK, "pales.vtp")])
dp = Show(pales, rv)
ColorBy(dp, ("CELLS", "p"))
dp.RescaleTransferFunctionToDataRange(True, False)
dp.SetScalarBarVisibility(rv, True)
pLUT = GetColorTransferFunction("p"); pLUT.ApplyPreset("Cool to Warm", True)
sb2 = GetScalarBar(pLUT, rv); sb2.Title = "p (m2/s2)"; sb2.ComponentTitle = ""
sb2.TitleColor = [0, 0, 0]; sb2.LabelColor = [0, 0, 0]
ResetCamera(rv)
rv.CameraPosition = [0.02, -0.28, 0.14]
rv.CameraFocalPoint = [0.0, 0.0, 0.0]
rv.CameraViewUp = [0, 0, 1]
Render()
SaveScreenshot(os.path.join(BASE, "flow_6_pales_pression.png"), rv, ImageResolution=[1400, 950])
print("OK 3D")
