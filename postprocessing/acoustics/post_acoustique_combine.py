# -*- coding: utf-8 -*-
"""Spectre intermediaire : concatenation des 2 enregistrements FW-H -> Delta f plus fin."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import medfilt, detrend

BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
F1 = BASE + r"\FWH_helice-time.dat"      # 0.19 -> 0.29
F2 = BASE + r"\FWH_helice-time_2.dat"    # 0.29 -> 0.374 (partiel)
PREF = 20e-6; BPF = 125.7; SHAFT = 25.15; T_CUT = 0.015
mics = ["aval_1m","amont_1m","radial_1m","diag45_1m"]
coul = {"aval_1m":"#185FA5","amont_1m":"#D85A30","radial_1m":"#0F6E56","diag45_1m":"#993556"}

d1 = np.loadtxt(F1, skiprows=1); d2 = np.loadtxt(F2, skiprows=1)
dt = np.median(np.diff(d1[:,0])); fs = 1/dt

# garde la partie valide de chaque segment (apres transitoire), retire la moyenne, concatene
sig = {m: [] for m in mics}
for d in (d1, d2):
    m1 = d[:,0] >= T_CUT
    for i,m in enumerate(mics):
        s = d[m1, i+1]; sig[m].append(s - s.mean())
for m in mics:
    sig[m] = np.concatenate(sig[m])
N = len(sig["aval_1m"]); T = N*dt
print(f"signal combine : {N} pts, {T*1000:.0f} ms, fs={fs:.0f}Hz, df={1/T:.1f}Hz")

clean = {m: detrend(medfilt(sig[m],5), type="linear") for m in mics}

# === spectre ===
fig, ax = plt.subplots(figsize=(11,5.2))
win = np.hanning(N); cg = win.mean(); f = np.fft.rfftfreq(N, d=dt)
res = {}
for m in mics:
    X = np.fft.rfft(clean[m]*win); amp = 2*np.abs(X)/N/cg; amp[0]/=2
    spl = 20*np.log10(amp/np.sqrt(2)/PREF + 1e-30)
    ax.plot(f, spl, lw=1.1, color=coul[m], label=m)
    def lvl(fc, bw=6):
        w = np.abs(f-fc)<bw; return spl[w].max() if w.any() else np.nan
    res[m] = {"shaft":lvl(SHAFT),"BPF":lvl(BPF),"2BPF":lvl(2*BPF),"3BPF":lvl(3*BPF),
              "OASPL":20*np.log10(clean[m].std()/PREF)}
for k in range(1,5):
    ax.axvline(BPF*k, color="k", ls="--", lw=0.6, alpha=.5)
    ax.text(BPF*k, 22, f"{k}xBPF" if k>1 else "BPF 126Hz", rotation=90, fontsize=7, va="top", ha="right", alpha=.7)
ax.axvline(SHAFT, color="purple", ls=":", lw=0.9, alpha=.7)
ax.text(SHAFT, 22, "arbre 25Hz", rotation=90, fontsize=7, va="top", ha="right", color="purple", alpha=.8)
ax.set_xlim(0,700); ax.set_ylim(-45,28); ax.set_xlabel("frequence (Hz)"); ax.set_ylabel("SPL (dB ref 20uPa)")
ax.set_title(f"Spectre SPL intermediaire - signal combine {T*1000:.0f} ms (df={1/T:.1f}Hz)")
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(BASE+r"\acou_6_spectre_intermediaire.png", dpi=125); plt.close()

print(f"\n{'micro':12}{'arbre25':>9}{'BPF126':>9}{'2xBPF':>9}{'3xBPF':>9}{'OASPL':>9}")
for m in mics:
    r=res[m]; print(f"{m:12}{r['shaft']:9.1f}{r['BPF']:9.1f}{r['2BPF']:9.1f}{r['3BPF']:9.1f}{r['OASPL']:9.1f}")
print("\nFigure: acou_6_spectre_intermediaire.png")
