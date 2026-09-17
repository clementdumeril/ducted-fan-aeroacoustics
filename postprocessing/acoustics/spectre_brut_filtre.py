# -*- coding: utf-8 -*-
"""Comparaison spectre BRUT vs FILTRE (transparence sur le post-traitement)."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import medfilt, detrend

BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
PREF = 20e-6; BPF = 125.7; T_CUT = 0.015
d1 = np.loadtxt(BASE + r"\FWH_helice-time.dat", skiprows=1)
d2 = np.loadtxt(BASE + r"\FWH_helice-time_2.dat", skiprows=1)
dt = np.median(np.diff(d1[:, 0])); fs = 1/dt
MIC = 3  # radial (col 3) = le plus contamine

def spectrum(sig):
    N = len(sig); win = np.hanning(N); cg = win.mean()
    X = np.fft.rfft(sig*win); amp = 2*np.abs(X)/N/cg; amp[0] /= 2
    f = np.fft.rfftfreq(N, d=dt)
    return f, 20*np.log10(amp/np.sqrt(2)/PREF + 1e-30)

# concatenation (brut : moyenne retiree par segment ; filtre : + median + detrend)
raw, clean = [], []
for d in (d1, d2):
    m = d[:, 0] >= T_CUT
    s = d[m, MIC]
    raw.append(s - s.mean())
    clean.append(detrend(medfilt(s, 5), type="linear"))
raw = np.concatenate(raw); clean = np.concatenate(clean)

fr, sr = spectrum(raw); fc, sc = spectrum(clean)
fig, ax = plt.subplots(figsize=(11, 5))
ax.plot(fr, sr, lw=0.9, color="#B0B0B0", label="spectre BRUT (moyenne retiree)")
ax.plot(fc, sc, lw=1.3, color="#0F6E56", label="spectre FILTRE (median 5pts + detrend)")
for k in range(1, 4):
    ax.axvline(BPF*k, color="k", ls="--", lw=0.6, alpha=.5)
    ax.text(BPF*k, 24, f"{k}xBPF" if k > 1 else "BPF", rotation=90, fontsize=7, va="top", ha="right", alpha=.7)
ax.axvline(25.15, color="purple", ls=":", lw=0.9, alpha=.7); ax.text(25, 24, "25Hz", rotation=90, fontsize=7, va="top", ha="right", color="purple")
ax.set_xlim(0, 700); ax.set_ylim(-45, 28); ax.set_xlabel("frequence (Hz)"); ax.set_ylabel("SPL (dB ref 20uPa)")
ax.set_title("Micro radial : spectre brut vs filtre (le filtrage n'invente pas les raies, il reduit le bruit)")
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(BASE + r"\acou_7_brut_vs_filtre.png", dpi=130); plt.close()
print("OK. Figure: acou_7_brut_vs_filtre.png")
print(f"protocole FFT : fenetre Hann, N={len(clean)} pts, fs={fs:.0f}Hz, 1 bloc (pas de Welch), single-sided, sans ponderation A")
