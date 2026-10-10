"""High-frequency line comb in the FW-H signals and periodic bursts in the blade force.

1. Lines above 600 Hz in the per-segment averaged spectra (same processing as
   postprocessing/acoustics/post_acoustique_combine.py): frequency and height above the local floor.
2. Envelope of the high-passed (> 600 Hz) signals folded on a 10.00 ms period (the field write interval,
   writeInterval 0.01 s) and on 9.94 ms (a quarter revolution at 25.15 Hz): peak-to-mean ratio of the folded
   envelope, and phase of its maximum for the axial force F_z and the downstream FW-H signal.
Writes paper/figures/comb_fold.pdf.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, detrend, find_peaks, hilbert, sosfiltfilt

REPO = Path(__file__).resolve().parents[2]
PREF, DT, NFFT, T_CUT = 20e-6, 1.75e-4, 8192, 0.015
recs = [(np.loadtxt(REPO / "data" / n, skiprows=1), t0) for n, t0 in
        (("FWH_helice-time.dat", 0.19), ("FWH_helice-time_2.dat", 0.29))]
mics = ["downstream", "upstream", "radial", "45 deg"]

# 1. lines above 600 Hz
n = int((recs[0][0][:, 0] >= T_CUT).sum()); win = np.hanning(n); cg = win.mean(); f = np.fft.rfftfreq(NFFT, DT)
for i, m in enumerate(mics):
    pw = []
    for d, _ in recs:
        x = d[d[:, 0] >= T_CUT, i + 1]; x = detrend(x - x.mean())
        a = 2 * np.abs(np.fft.rfft(x * win, NFFT)) / n / cg; pw.append((a / np.sqrt(2)) ** 2)
    s = 10 * np.log10(np.mean(pw, 0) / PREF ** 2); fl = np.median(s[(f > 140) & (f < 240)])
    k, _ = find_peaks(s, prominence=8); k = k[(f[k] > 600) & (f[k] < 1500)]
    print(f"{m:10s} lines (Hz, dB above floor):", [(int(round(f[j])), round(float(s[j] - fl), 1)) for j in k])

# 2. folded envelopes
sos_f = butter(4, 600, btype="high", fs=1 / DT, output="sos")
F = np.loadtxt(REPO / "data" / "forces_data" / "force_0.19.dat"); tf = F[:, 0]; fz = F[:, 6] - F[:, 6].mean()
sos_z = butter(4, 600, btype="high", fs=1 / 3.5e-4, output="sos")
env_z = np.abs(hilbert(sosfiltfilt(sos_z, fz))) ** 2
NB = 25


def fold(t, e, T):
    b = ((t % T) / T * NB).astype(int) % NB
    p = np.array([e[b == j].mean() for j in range(NB)])
    return p / p.mean()


fig, axs = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
for ax, T in zip(axs, (0.010, 0.25 / 25.15)):
    pz = fold(tf, env_z, T)
    pw = np.mean([fold(t0 + d[d[:, 0] >= T_CUT, 0], np.abs(hilbert(sosfiltfilt(sos_f, d[d[:, 0] >= T_CUT, 1]))) ** 2, T)
                  for d, t0 in recs], axis=0)
    ph = (np.arange(NB) + 0.5) / NB * T * 1e3
    ax.plot(ph, pz, "o-", ms=3, label="blade force $F_z$ (pressure)")
    ax.plot(ph, pw, "s-", ms=3, label="FW-H, downstream")
    ax.set_xlabel(f"phase within a {T*1e3:.2f} ms period (ms)")
    ax.grid(alpha=0.3)
    lag = (ph[np.argmax(pw)] - ph[np.argmax(pz)]) % (T * 1e3)
    print(f"fold T = {T*1e3:.2f} ms: F_z peak/mean {pz.max():.1f} at {ph[np.argmax(pz)]:.1f} ms, "
          f"FW-H peak/mean {pw.max():.1f} at {ph[np.argmax(pw)]:.1f} ms, lag {lag:.1f} ms (r/c0 = 2.9 ms)")
axs[0].set_ylabel("folded high-pass envelope / mean"); axs[0].legend(fontsize=8)
fig.tight_layout()
out = REPO / "paper" / "figures"; out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "comb_fold.pdf"); print("written", out / "comb_fold.pdf")
