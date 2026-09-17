# -*- coding: utf-8 -*-
"""Post-traitement acoustique FW-H v2 : de-spike + retrait derive -> isole le tonal (BPF)."""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import medfilt, detrend, welch

BASE = r"C:\Users\PC\OneDrive - ensam.eu\Bureau\polytechnia\simulation"
DAT  = BASE + r"\FWH_helice-time.dat"
PREF = 20e-6; BPF = 125.7; T_CUT = 0.02
mics = ["aval_1m","amont_1m","radial_1m","diag45_1m"]
coul = {"aval_1m":"#185FA5","amont_1m":"#D85A30","radial_1m":"#0F6E56","diag45_1m":"#993556"}

data = np.loadtxt(DAT, skiprows=1); t = data[:,0]
dt = np.median(np.diff(t)); fs = 1/dt
mask = t >= T_CUT; ts = t[mask]; N = mask.sum()
print(f"fs={fs:.0f}Hz Nyq={fs/2:.0f}Hz  analyse {ts[0]:.3f}->{ts[-1]:.3f}s ({N}pts) df={1/(ts[-1]-ts[0]):.1f}Hz")

# nettoyage : de-spike (median 5pts) puis retrait derive (detrend lineaire)
clean = {}
for m in mics:
    sig = data[mask, mics.index(m)+1]
    sig_ds = medfilt(sig, 5)             # enleve les spikes isoles
    clean[m] = detrend(sig_ds, type="linear")  # enleve la derive

# === FIG : signal brut vs nettoye (radial, le + contamine) ===
fig, ax = plt.subplots(figsize=(11,4))
ax.plot(ts*1000, (data[mask, 3]-data[mask,3].mean())*1000, lw=0.6, color="gray", label="radial brut", alpha=.6)
ax.plot(ts*1000, clean["radial_1m"]*1000, lw=0.9, color="#0F6E56", label="radial nettoye (de-spike + detrend)")
ax.set_xlabel("temps (ms)"); ax.set_ylabel("p' (mPa)"); ax.set_title("Nettoyage du signal")
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(BASE+r"\acou_4_nettoyage.png", dpi=120); plt.close()

# === FIG : spectre nettoye (Welch + FFT) ===
fig, ax = plt.subplots(figsize=(11,5))
win = np.hanning(N); cg = win.mean(); f = np.fft.rfftfreq(N, d=dt)
res = {}
for m in mics:
    X = np.fft.rfft(clean[m]*win); amp = 2*np.abs(X)/N/cg; amp[0]/=2
    spl = 20*np.log10(amp/np.sqrt(2)/PREF + 1e-30)
    ax.plot(f, spl, lw=1.0, color=coul[m], label=m)
    # niveau aux raies
    def lvl(fc):
        w = np.abs(f-fc)<8; return spl[w].max() if w.any() else np.nan
    res[m] = {"BPF":lvl(BPF), "2BPF":lvl(2*BPF), "3BPF":lvl(3*BPF), "shaft":lvl(25.15),
              "OASPL":20*np.log10(clean[m].std()/PREF)}
for k in range(1,5):
    ax.axvline(BPF*k, color="k", ls="--", lw=0.6, alpha=.5)
    ax.text(BPF*k, 18, f"{k}xBPF" if k>1 else "BPF 126Hz", rotation=90, fontsize=7, va="top", ha="right", alpha=.7)
ax.axvline(25.15, color="purple", ls=":", lw=0.8, alpha=.6); ax.text(25, 18, "arbre 25Hz", rotation=90, fontsize=7, va="top", ha="right", color="purple", alpha=.7)
ax.set_xlim(0,700); ax.set_ylim(-50,25); ax.set_xlabel("frequence (Hz)"); ax.set_ylabel("SPL (dB ref 20uPa)")
ax.set_title("Spectre SPL nettoye (de-spike + detrend) - raies tonales")
ax.grid(alpha=.3); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(BASE+r"\acou_5_spectre_propre.png", dpi=120); plt.close()

print("\n=== NIVEAUX TONAUX (dB, signal nettoye) ===")
print(f"{'micro':12}{'arbre25':>9}{'BPF126':>9}{'2xBPF':>9}{'3xBPF':>9}{'OASPL':>9}")
for m in mics:
    r=res[m]; print(f"{m:12}{r['shaft']:9.1f}{r['BPF']:9.1f}{r['2BPF']:9.1f}{r['3BPF']:9.1f}{r['OASPL']:9.1f}")
print("\nFigures: acou_4_nettoyage.png, acou_5_spectre_propre.png")
