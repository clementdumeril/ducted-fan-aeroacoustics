# -*- coding: utf-8 -*-
"""Far-field SPL spectrum, tone levels, directivity and raw-vs-filtered check from the FW-H signals.

The two FW-H records in data/ (sim time 0.19-0.29 s and 0.29-0.39 s) each lose their first 15 ms
(FW-H start-up), which leaves two 85 ms segments separated by a 15 ms gap. They are NOT concatenated:
each segment gets its own Hann-windowed, zero-padded spectrum and the two power spectra are averaged.
Tone levels are the spectrum peaks within +-6 Hz of the shaft frequency and of the BPF; the local
floor is the median level between 140 and 240 Hz. Results use the raw signal (mean and linear trend
removed per segment). A 5-point median filter, used in an earlier version, is kept only as a check:
the few outlier samples in the segments do not need it, and it raises the floor by 9 to 24 dB and
shifts the radial BPF level by about 4 dB.
Figures are written to figures/acoustics/.
"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import medfilt, detrend

REPO = Path(__file__).resolve().parents[2]
RECORDS = [(REPO / "data" / "FWH_helice-time.dat", 0.19),     # (file, sim time at the start of the record)
           (REPO / "data" / "FWH_helice-time_2.dat", 0.29)]
OUT = REPO / "figures" / "acoustics"
PREF = 20e-6; BPF = 125.7; SHAFT = 25.15; T_CUT = 0.015; NFFT = 8192
mics = ["aval_1m", "amont_1m", "radial_1m", "diag45_1m"]
names = {"aval_1m": "downstream", "amont_1m": "upstream", "radial_1m": "radial", "diag45_1m": "45 deg"}
colors = {"aval_1m": "#185FA5", "amont_1m": "#D85A30", "radial_1m": "#0F6E56", "diag45_1m": "#993556"}

segments = []
for path, t0 in RECORDS:
    d = np.loadtxt(path, skiprows=1)
    keep = d[:, 0] >= T_CUT
    segments.append(d[keep, 1:])
dt = np.median(np.diff(np.loadtxt(RECORDS[0][0], skiprows=1)[:, 0]))
n = len(segments[0])
print(f"{len(segments)} segments of {n} samples ({n*dt*1000:.0f} ms each), fs = {1/dt:.0f} Hz, "
      f"resolution {1/(n*dt):.1f} Hz (zero-padded to {NFFT})")

win = np.hanning(n); cg = win.mean(); f = np.fft.rfftfreq(NFFT, d=dt)


def spl_avg(i, filtered):
    """Average of the per-segment power spectra for microphone column i, in dB re 20 uPa."""
    power = []
    for seg in segments:
        x = seg[:, i] - seg[:, i].mean()
        if filtered:
            x = medfilt(x, 5)          # removes isolated numerical spikes
        x = detrend(x, type="linear")
        amp = 2 * np.abs(np.fft.rfft(x * win, NFFT)) / n / cg
        power.append((amp / np.sqrt(2)) ** 2)
    return 10 * np.log10(np.mean(power, axis=0) / PREF**2 + 1e-30)


def peak(spl, fc, bw=6):
    w = np.abs(f - fc) <= bw
    return spl[w].max()


def floor(spl):
    return np.median(spl[(f > 140) & (f < 240)])


res = {}
for i, m in enumerate(mics):
    for kind in ("raw", "filtered"):
        s = spl_avg(i, kind == "filtered")
        res[m, kind] = {"spl": s, "shaft": peak(s, SHAFT), "BPF": peak(s, BPF), "floor": floor(s)}

# --- SPL spectrum (raw, average of the two segments) ---
fig, ax = plt.subplots(figsize=(11, 5.2))
for m in mics:
    ax.plot(f, res[m, "raw"]["spl"], lw=1.1, color=colors[m], label=f"{names[m]} (1 m)")
for k in range(1, 5):
    ax.axvline(BPF * k, color="k", ls="--", lw=0.6, alpha=.5)
    ax.text(BPF * k, 27, f"{k}xBPF" if k > 1 else "BPF 126 Hz", rotation=90, fontsize=8, va="top", ha="right", alpha=.7)
ax.axvline(SHAFT, color="purple", ls=":", lw=0.9, alpha=.7)
ax.text(SHAFT, 27, "shaft 25 Hz", rotation=90, fontsize=8, va="top", ha="right", color="purple", alpha=.8)
ax.set_xlim(0, 700); ax.set_ylim(-50, 30)
ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("SPL (dB re 20 µPa)")
ax.set_title("Far-field SPL from FW-H at 1 m: average of two 85 ms segments (Hann, zero-padded)")
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(OUT / "spl_spectrum.png", dpi=125); plt.close()

# --- directivity: shaft-frequency and BPF levels and local floor per microphone (raw) ---
fig, ax = plt.subplots(figsize=(9, 4.8))
xk = np.arange(len(mics)); w = 0.26
series = [("level at 25 Hz", "raw", "shaft", "#185FA5"), ("BPF 126 Hz", "raw", "BPF", "#D85A30"),
          ("local floor 140-240 Hz", "raw", "floor", "#9E9E9E")]
for j, (lab, kind, key, col) in enumerate(series):
    vals = [res[m, kind][key] for m in mics]
    bars = ax.bar(xk + (j - 1) * w, vals, w, color=col, label=lab)
    for b in bars:
        h = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, h + (0.6 if h >= 0 else -2.6), f"{h:.0f}", ha="center", fontsize=7)
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(xk, [names[m] for m in mics]); ax.set_ylabel("Level (dB re 20 µPa)")
ax.set_ylim(-48, 28)
ax.set_title("Directivity at 1 m: 25 Hz tone in the rotor plane only (on the axis it is at the floor); BPF everywhere")
ax.legend(fontsize=8, ncol=2, loc="upper left"); ax.grid(axis="y", alpha=.3)
plt.tight_layout(); plt.savefig(OUT / "directivity.png", dpi=125); plt.close()

# --- check: raw vs filtered spectrum on the radial microphone ---
m = "radial_1m"
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.plot(f, res[m, "raw"]["spl"], color=colors[m], lw=1.3, label="raw (mean and trend removed), used for the results")
ax.plot(f, res[m, "filtered"]["spl"], color="#C0392B", lw=1.0, label="with a 5-point median filter (not used)")
for k in range(1, 4):
    ax.axvline(BPF * k, color="k", ls="--", lw=0.6, alpha=.5)
ax.set_xlim(0, 700); ax.set_ylim(-50, 30)
ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("SPL (dB re 20 µPa)")
ax.set_title("Median-filter check, radial microphone: same 25 Hz tone, but the filter raises the floor by about 24 dB")
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(OUT / "raw_vs_filtered.png", dpi=125); plt.close()

print(f"\n{'mic':12}{'':10}{'25 Hz':>8}{'BPF':>8}{'floor':>8}{'BPF-floor':>11}")
for m in mics:
    for kind in ("raw", "filtered"):
        r = res[m, kind]
        print(f"{m:12}{kind:10}{r['shaft']:8.1f}{r['BPF']:8.1f}{r['floor']:8.1f}{r['BPF']-r['floor']:11.1f}")
print(f"\nFigures written to {OUT}")
