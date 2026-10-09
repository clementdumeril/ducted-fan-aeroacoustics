# -*- coding: utf-8 -*-
"""Far-field SPL spectrum, directivity and raw-vs-filtered check from the FW-H microphone signals.

Concatenates the two FW-H records in data/ (after the start-up transient of each segment),
applies a Hann window and prints the levels at the shaft frequency, BPF and harmonics, plus OASPL.
Figures are written to figures/acoustics/.
"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import medfilt, detrend

REPO = Path(__file__).resolve().parents[2]
F1 = REPO / "data" / "FWH_helice-time.dat"      # sim time 0.19 -> 0.29 s
F2 = REPO / "data" / "FWH_helice-time_2.dat"    # sim time 0.29 -> 0.374 s (partial)
OUT = REPO / "figures" / "acoustics"
PREF = 20e-6; BPF = 125.7; SHAFT = 25.15; T_CUT = 0.015
mics = ["aval_1m", "amont_1m", "radial_1m", "diag45_1m"]
names = {"aval_1m": "downstream", "amont_1m": "upstream", "radial_1m": "radial", "diag45_1m": "45 deg"}
colors = {"aval_1m": "#185FA5", "amont_1m": "#D85A30", "radial_1m": "#0F6E56", "diag45_1m": "#993556"}

d1 = np.loadtxt(F1, skiprows=1); d2 = np.loadtxt(F2, skiprows=1)
dt = np.median(np.diff(d1[:, 0])); fs = 1 / dt

# keep the valid part of each segment, remove its mean, concatenate
sig = {m: [] for m in mics}
for d in (d1, d2):
    keep = d[:, 0] >= T_CUT
    for i, m in enumerate(mics):
        s = d[keep, i + 1]; sig[m].append(s - s.mean())
for m in mics:
    sig[m] = np.concatenate(sig[m])
N = len(sig["aval_1m"]); T = N * dt
print(f"combined signal: {N} samples, {T*1000:.0f} ms, fs = {fs:.0f} Hz, df = {1/T:.1f} Hz")

# 5-point median filter removes isolated numerical spikes, then linear detrend
clean = {m: detrend(medfilt(sig[m], 5), type="linear") for m in mics}
win = np.hanning(N); cg = win.mean(); f = np.fft.rfftfreq(N, d=dt)


def spl_of(x):
    X = np.fft.rfft(x * win); amp = 2 * np.abs(X) / N / cg; amp[0] /= 2
    return 20 * np.log10(amp / np.sqrt(2) / PREF + 1e-30)


def level(spl, fc, bw=6):
    w = np.abs(f - fc) < bw
    return spl[w].max() if w.any() else np.nan


# --- SPL spectrum ---
fig, ax = plt.subplots(figsize=(11, 5.2))
res = {}
for m in mics:
    spl = spl_of(clean[m])
    ax.plot(f, spl, lw=1.1, color=colors[m], label=f"{names[m]} (1 m)")
    res[m] = {"shaft": level(spl, SHAFT), "BPF": level(spl, BPF), "2BPF": level(spl, 2 * BPF),
              "3BPF": level(spl, 3 * BPF), "OASPL": 20 * np.log10(clean[m].std() / PREF)}
for k in range(1, 5):
    ax.axvline(BPF * k, color="k", ls="--", lw=0.6, alpha=.5)
    ax.text(BPF * k, 22, f"{k}xBPF" if k > 1 else "BPF 126 Hz", rotation=90, fontsize=8, va="top", ha="right", alpha=.7)
ax.axvline(SHAFT, color="purple", ls=":", lw=0.9, alpha=.7)
ax.text(SHAFT, 22, "shaft 25 Hz", rotation=90, fontsize=8, va="top", ha="right", color="purple", alpha=.8)
ax.set_xlim(0, 700); ax.set_ylim(-45, 28)
ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("SPL (dB re 20 µPa)")
ax.set_title(f"Far-field SPL from FW-H, 4 microphones at 1 m ({T*1000:.0f} ms record, Δf = {1/T:.1f} Hz)")
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(OUT / "spl_spectrum.png", dpi=125); plt.close()

# levels from the raw signal (mean removed, detrended), to check what the median filter changes
raw = {}
for m in mics:
    s = spl_of(detrend(sig[m], type="linear"))
    raw[m] = {"shaft": level(s, SHAFT), "BPF": level(s, BPF), "2BPF": level(s, 2 * BPF),
              "OASPL": 20 * np.log10(detrend(sig[m], type="linear").std() / PREF)}

# --- directivity: once-per-revolution and BPF tones per microphone, raw and filtered ---
fig, ax = plt.subplots(figsize=(9, 4.8))
xk = np.arange(len(mics)); w = 0.2
series = [("1/rev 25 Hz, raw", [raw[m]["shaft"] for m in mics], "#9BB8D8"),
          ("1/rev 25 Hz, filtered", [res[m]["shaft"] for m in mics], "#185FA5"),
          ("BPF 126 Hz, raw", [raw[m]["BPF"] for m in mics], "#F0A88C"),
          ("BPF 126 Hz, filtered", [res[m]["BPF"] for m in mics], "#D85A30")]
for j, (lab, vals, col) in enumerate(series):
    bars = ax.bar(xk + (j - 1.5) * w, vals, w, color=col, label=lab)
    for b in bars:
        h = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, h + (0.6 if h >= 0 else -2.4), f"{h:.0f}", ha="center", fontsize=7)
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(xk, [names[m] for m in mics]); ax.set_ylabel("Tone level (dB re 20 µPa)")
ax.set_ylim(-30, 25)
ax.set_title("Directivity at 1 m: the 1/rev tone radiates in the rotor plane, the BPF tone is filter-independent")
ax.legend(fontsize=8, ncol=2, loc="upper left"); ax.grid(axis="y", alpha=.3)
plt.tight_layout(); plt.savefig(OUT / "directivity.png", dpi=125); plt.close()

# --- check: raw vs filtered spectrum on the radial microphone ---
m = "radial_1m"
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.plot(f, spl_of(detrend(sig[m], type="linear")), color="0.6", lw=1.0, label="raw (mean removed)")
ax.plot(f, spl_of(clean[m]), color=colors[m], lw=1.2, label="5-point median filter + detrend")
for k in range(1, 4):
    ax.axvline(BPF * k, color="k", ls="--", lw=0.6, alpha=.5)
ax.set_xlim(0, 700); ax.set_ylim(-45, 28)
ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("SPL (dB re 20 µPa)")
ax.set_title("Filter check, radial microphone: 1/rev and BPF agree within 2 dB; above 130 Hz the filter adds artifacts")
ax.grid(alpha=.3); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(OUT / "raw_vs_filtered.png", dpi=125); plt.close()

print(f"\n{'mic':12}{'shaft25':>9}{'BPF126':>9}{'2xBPF':>9}{'3xBPF':>9}{'OASPL':>9}   (filtered)")
for m in mics:
    r = res[m]; print(f"{m:12}{r['shaft']:9.1f}{r['BPF']:9.1f}{r['2BPF']:9.1f}{r['3BPF']:9.1f}{r['OASPL']:9.1f}")
print(f"\n{'mic':12}{'shaft25':>9}{'BPF126':>9}{'2xBPF':>9}{'OASPL':>9}   (raw)")
for m in mics:
    r = raw[m]; print(f"{m:12}{r['shaft']:9.1f}{r['BPF']:9.1f}{r['2BPF']:9.1f}{r['OASPL']:9.1f}")
print(f"\nFigures written to {OUT}")
