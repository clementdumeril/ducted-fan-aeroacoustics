# -*- coding: utf-8 -*-
"""Post-traitement acoustique FW-H : signal p'(t) aux 4 micros -> spectre SPL, BPF, OASPL, directivite."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
DAT  = BASE + r"\FWH_helice-time.dat"
PREF = 20e-6          # reference SPL (Pa)
BPF  = 125.7          # frequence de passage des pales (5 pales x 25.15 tr/s)
T_CUT = 0.02          # on jette le transitoire de demarrage (s)

# --- lecture ---
mics = ["aval_1m", "amont_1m", "radial_1m", "diag45_1m"]
data = np.loadtxt(DAT, skiprows=1)
t = data[:, 0]
P = {m: data[:, i+1] for i, m in enumerate(mics)}
dt = np.median(np.diff(t)); fs = 1.0/dt
print(f"echantillons={len(t)}  dt={dt:.3e}s  fs={fs:.1f}Hz  Nyquist={fs/2:.0f}Hz  duree={t[-1]:.3f}s")

# --- partie stationnaire (apres transitoire) ---
mask = t >= T_CUT
ts = t[mask]
N = mask.sum()
win = np.hanning(N); cg = win.mean()
Tstat = ts[-1]-ts[0]
print(f"partie analysee: {ts[0]:.3f} -> {ts[-1]:.3f}s ({N} pts)  -> resolution df={1/Tstat:.1f}Hz")

couleurs = {"aval_1m":"#185FA5", "amont_1m":"#D85A30", "radial_1m":"#0F6E56", "diag45_1m":"#993556"}

# === FIG 1 : signaux temporels ===
fig, ax = plt.subplots(figsize=(11,4.5))
for m in mics:
    ax.plot(t*1000, P[m]*1000, lw=0.8, color=couleurs[m], label=m)
ax.axvspan(0, T_CUT*1000, color="gray", alpha=0.15, label="transitoire (ignore)")
ax.set_xlabel("temps (ms)"); ax.set_ylabel("pression acoustique p' (mPa)")
ax.set_title("Signal acoustique FW-H aux 4 micros (1 m)")
ax.grid(alpha=.3); ax.legend(fontsize=8, ncol=5)
plt.tight_layout(); plt.savefig(BASE+r"\acou_1_signal.png", dpi=120); plt.close()

# === FIG 2 : spectre SPL ===
fig, ax = plt.subplots(figsize=(11,5))
f = np.fft.rfftfreq(N, d=dt)
spl_bpf = {}; oaspl = {}
for m in mics:
    sig = P[m][mask] - P[m][mask].mean()
    X = np.fft.rfft(sig*win)
    amp = 2*np.abs(X)/N/cg; amp[0] /= 2
    prms = amp/np.sqrt(2)
    spl = 20*np.log10(prms/PREF + 1e-30)
    ax.plot(f, spl, lw=1.0, color=couleurs[m], label=m)
    # SPL au pic BPF (max dans une fenetre autour de 126 Hz)
    ib = np.argmin(np.abs(f-BPF)); w = np.abs(f-BPF) < 15
    spl_bpf[m] = spl[w].max()
    oaspl[m] = 20*np.log10(sig.std()/PREF)
for k in range(1,6):
    ax.axvline(BPF*k, color="k", ls="--", lw=0.6, alpha=0.5)
    ax.text(BPF*k, ax.get_ylim()[1]-3, f"{k}xBPF" if k>1 else "BPF", rotation=90, fontsize=7, va="top", ha="right", alpha=.7)
ax.set_xlim(0, 1500); ax.set_xlabel("frequence (Hz)"); ax.set_ylabel("SPL (dB ref 20 uPa)")
ax.set_title(f"Spectre SPL aux 4 micros - BPF = {BPF:.0f} Hz (raies tonales)")
ax.grid(alpha=.3); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(BASE+r"\acou_2_spectre.png", dpi=120); plt.close()

# === FIG 3 : directivite (OASPL + niveau BPF par micro) ===
fig, ax = plt.subplots(figsize=(8,5))
x = np.arange(len(mics)); w=0.38
ax.bar(x-w/2, [oaspl[m] for m in mics], w, color="#378ADD", label="OASPL (global)")
ax.bar(x+w/2, [spl_bpf[m] for m in mics], w, color="#D85A30", label="niveau BPF (126 Hz)")
for i,m in enumerate(mics):
    ax.text(i-w/2, oaspl[m]+0.3, f"{oaspl[m]:.1f}", ha="center", fontsize=8)
    ax.text(i+w/2, spl_bpf[m]+0.3, f"{spl_bpf[m]:.1f}", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(mics); ax.set_ylabel("SPL (dB)")
ax.set_title("Directivite : niveau sonore selon la position du micro")
ax.grid(alpha=.3, axis="y"); ax.legend()
plt.tight_layout(); plt.savefig(BASE+r"\acou_3_directivite.png", dpi=120); plt.close()

# === RESUME ===
print("\n=== RESULTATS ACOUSTIQUES ===")
print(f"{'micro':12} {'OASPL(dB)':>10} {'BPF 126Hz(dB)':>14}")
for m in mics:
    print(f"{m:12} {oaspl[m]:10.1f} {spl_bpf[m]:14.1f}")
print(f"\nmicro le + bruyant (OASPL) : {max(oaspl,key=oaspl.get)} ({max(oaspl.values()):.1f} dB)")
print("Figures: acou_1_signal.png, acou_2_spectre.png, acou_3_directivite.png")
