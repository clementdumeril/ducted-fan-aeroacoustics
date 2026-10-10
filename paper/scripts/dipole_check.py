"""Cross-check of the FW-H tone levels against a compact-dipole estimate built from the blade forces.

At 25 Hz (wavelength 13.6 m) and at the BPF (2.7 m) the rotor (tip radius 69 mm) is acoustically compact,
so the loading noise of the whole rotor reduces to that of the net force F(t) applied at the hub:
    p'(x, t) = 1 / (4 pi) * [ x_i / (c0 r^2) dF_i/dt + x_i F_i / r^3 ]   evaluated at t - r / c0,
(far-field and near-field dipole terms; the near-field term matters at 25 Hz, where k r = 0.46).
For a harmonic force component of amplitude A at angular frequency w, along the observer direction:
    |p'| = A / (4 pi r) * sqrt((w / c0)^2 + 1 / r^2).
Force amplitudes come from a least-squares sine fit on data/forces_data/force_0.19.dat (0.19-0.29 s,
the same window as the first FW-H record); FW-H levels from postprocessing/acoustics (raw signal).
Writes paper/tables/dipole_check.tex.
"""
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
C0, R, PREF = 343.0, 1.0, 20e-6
f = np.loadtxt(REPO / "data" / "forces_data" / "force_0.19.dat")
t = f[:, 0] - f[0, 0]
F = {"x": f[:, 1], "y": f[:, 2], "z": f[:, 3]}


def amp(x, f0):
    A = np.column_stack([np.ones_like(t), t, np.cos(2 * np.pi * f0 * t), np.sin(2 * np.pi * f0 * t)])
    c, *_ = np.linalg.lstsq(A, x, rcond=None)
    return np.hypot(c[2], c[3])


def spl(a, f0):
    p = a / (4 * np.pi * R) * np.sqrt((2 * np.pi * f0 / C0) ** 2 + 1 / R ** 2)
    return 20 * np.log10(p / np.sqrt(2) / PREF)


grid = np.arange(110.0, 140.0, 0.5)
fz = [amp(F["z"], g) for g in grid]; kz = int(np.argmax(fz))
rows = [
    ("25 Hz, radial (+x)", 25.15, amp(F["x"], 25.15), 21.2),
    ("25 Hz, axis", 25.15, None, None),
    (f"BPF ({grid[kz]:.1f} Hz), downstream axis", grid[kz], fz[kz], -9.8),
]
lines = [r"\begin{tabular}{lrrr}", r"\toprule",
         r"Tone, microphone & force amplitude (mN) & compact dipole (dB) & FW-H (dB) \\", r"\midrule"]
for name, f0, a, fwh in rows:
    if a is None:
        lines.append(rf"{name} & -- & no radiation & at floor \\")
        print(f"{name}: compact dipole gives no radiation (force in the rotor plane)")
        continue
    est = spl(a, f0)
    lines.append(rf"{name} & {a*1e3:.2f} & {est:.1f} & {fwh:.1f} \\")
    print(f"{name}: |F| = {a*1e3:.3f} mN -> compact dipole {est:.1f} dB, FW-H {fwh:.1f} dB, difference {fwh-est:+.1f} dB")
print(f"in-plane force at 25 Hz: Fx {amp(F['x'],25.15)*1e3:.3f} mN, Fy {amp(F['y'],25.15)*1e3:.3f} mN")
lines += [r"\bottomrule", r"\end{tabular}"]
out = REPO / "paper" / "tables"; out.mkdir(parents=True, exist_ok=True)
(out / "dipole_check.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("written", out / "dipole_check.tex")
