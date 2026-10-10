"""Verification of the FW-H output against a compact-dipole estimate built from the blade forces.

At 0 Hz, 25 Hz (wavelength 13.6 m) and the BPF (2.7 m) the rotor (tip radius 69 mm) is acoustically compact,
so its loading noise reduces to that of the net force F(t) on the blades, applied at the hub:
    p'(x, t) = 1 / (4 pi) * [ x_i / (c0 r^2) dF_i/dt + x_i F_i / r^3 ]   at retarded time t - r / c0
(far-field and near-field dipole terms; at 25 Hz k r = 0.46, so the near-field term dominates). For a
harmonic component of amplitude A along the observer direction: |p'| = A / (4 pi r) sqrt((w/c0)^2 + 1/r^2);
for the mean force: p' = F / (4 pi r^2). The FW-H surface integrates the blade pressure only, so the
pressure part of the force output is used (columns pressure_x/y/z of data/forces_data/force_0.19.dat,
0.19-0.29 s, the same window as the first FW-H record).
FW-H values: mean of the raw signal per microphone, and tone levels from
postprocessing/acoustics/post_acoustique_combine.py (raw signal, per-segment spectra).
Writes paper/tables/dipole_check.tex.
"""
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
C0, R, PREF = 343.0, 1.0, 20e-6
f = np.loadtxt(REPO / "data" / "forces_data" / "force_0.19.dat")
t = f[:, 0] - f[0, 0]
P = {"x": f[:, 4], "y": f[:, 5], "z": f[:, 6]}          # pressure part of the blade force (N)
fwh = np.loadtxt(REPO / "data" / "FWH_helice-time.dat", skiprows=1)
fwh = fwh[fwh[:, 0] >= 0.015]
FWH_TONES = {"25 Hz, radial": 21.2, "BPF, downstream": -9.8}   # dB, from post_acoustique_combine.py


def amp(x, f0):
    A = np.column_stack([np.ones_like(t), t, np.cos(2 * np.pi * f0 * t), np.sin(2 * np.pi * f0 * t)])
    c, *_ = np.linalg.lstsq(A, x, rcond=None)
    return np.hypot(c[2], c[3])


def db(p_amp):
    return 20 * np.log10(p_amp / np.sqrt(2) / PREF)


rows = []
# mean (static) loading: steady axial force, observer downstream on the axis
fz_mean = abs(P["z"].mean())
est = fz_mean / (4 * np.pi * R ** 2)
meas = fwh[:, 1].mean()
rows.append(("Mean, downstream axis", f"{fz_mean*1e3:.1f}", f"{est*1e3:.2f} mPa", f"{meas*1e3:.2f} mPa",
             20 * np.log10(abs(meas) / est)))
print(f"mean: |Fz_p| = {fz_mean*1e3:.2f} mN -> {est*1e3:.2f} mPa; FW-H mean (mPa) down/up/radial/45: "
      f"{np.round(fwh[:, 1:].mean(0)*1e3, 3)}; ratio {abs(meas)/est:.2f}")
# once-per-revolution: rotating in-plane force, observer in the rotor plane (+x)
ax = amp(P["x"], 25.15)
e = db(ax / (4 * np.pi * R) * np.sqrt((2 * np.pi * 25.15 / C0) ** 2 + 1 / R ** 2))
rows.append(("25 Hz, radial", f"{ax*1e3:.2f}", f"{e:.1f} dB", f"{FWH_TONES['25 Hz, radial']:.1f} dB",
             FWH_TONES["25 Hz, radial"] - e))
print(f"25 Hz: in-plane pressure force Fx {ax*1e3:.3f} mN, Fy {amp(P['y'], 25.15)*1e3:.3f} mN -> {e:.1f} dB")
# BPF: axial force fluctuation, observer downstream on the axis
grid = np.arange(110.0, 140.0, 0.5)
fz = [amp(P["z"], g) for g in grid]; k = int(np.argmax(fz))
e = db(fz[k] / (4 * np.pi * R) * np.sqrt((2 * np.pi * grid[k] / C0) ** 2 + 1 / R ** 2))
rows.append((f"BPF ({grid[k]:.1f} Hz), downstream axis", f"{fz[k]*1e3:.2f}", f"{e:.1f} dB",
             f"{FWH_TONES['BPF, downstream']:.1f} dB", FWH_TONES["BPF, downstream"] - e))
print(f"BPF: axial pressure force {fz[k]*1e3:.3f} mN at {grid[k]} Hz -> {e:.1f} dB")

lines = [r"\begin{tabular}{lrrrr}", r"\toprule",
         r"Component, microphone & force (mN) & compact dipole & FW-H & FW-H $-$ dipole (dB) \\", r"\midrule"]
for name, F, est_s, meas_s, gap in rows:
    lines.append(rf"{name} & {F} & {est_s} & {meas_s} & {gap:+.1f} \\")
    print(f"{name}: gap {gap:+.1f} dB")
lines += [r"\bottomrule", r"\end{tabular}"]
out = REPO / "paper" / "tables"; out.mkdir(parents=True, exist_ok=True)
(out / "dipole_check.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("written", out / "dipole_check.tex")
