# Methodology

Notes that go with the report (`paper/main.pdf`): the numerical choices and the practical problems I ran into.

## 1. Source of the noise

Tonal fan noise comes from the unsteady blade loading, a rotating dipole that radiates at BPF = *Z·n* (here 5 × 25.15 = 126 Hz) and its harmonics. To capture it I used URANS with a sliding mesh, not a frozen rotor.

## 2. The Ffowcs Williams–Hawkings analogy

The microphones at 1 m are outside the 0.4 m domain, and the incompressible solver carries no sound waves. FW-H records the flow on a surface *S* and integrates the wave equation analytically to each observer.

I used Farassat formulation 1A (`libAcoustics`) on the blade skin (`helice` patch, rotating surface), with c₀ = 343 m/s, ρ∞ = 1.2 kg/m³ and p_ref = 20 µPa. Quadrupole terms are not integrated, and *S* covers the blades only: no casing, arms or hub.

## 3. Mesh

Hex-dominant `snappyHexMesh` from a 10 mm base cell: blades at level 5 (0.31 mm), casing and arms level 3, hub level 2, plus refined rotor and wake regions. For the sliding interface, `snappyHexMesh` creates a `faceZone` on a cylinder (R = 72 mm) and `createBaffles` splits it into a cyclicAMI pair (about 35 k faces).

`checkMesh` on the retained mesh gives:

| Quantity | Value |
|---|---:|
| Cells / points | 1,949,584 / 2,205,401 |
| Max aspect ratio | 6.61 |
| Max / average non-orthogonality | 61.4 / 8.2 |
| Max skewness | 7.25 (12 faces above 4) |
| Concave cells | 49,827 (2.6 %) |
| Cells with determinant < 0.001 | 1851 |
| AMI weights (min / avg / max) | 0.90 / 1.002 / 1.27 |

The solver runs with `FOAM_SIGFPE=false` to tolerate the bad cells. That does not fix them, and I have not measured their effect on the spectrum.

There are no prism layers (`addLayers false`; an earlier attempt collapsed on the thin serrated blades). Blade y⁺ averages about 4 and reaches 445 locally, so the surface loading is preliminary.

## 4. Numerical model

- `pimpleFoam`, air ρ = 1.2 kg/m³, ν = 1.5·10⁻⁵ m²/s, k-ω SST with automatic wall treatment.
- Solid-body rotation ω = 158 rad/s, ramped from 0 over the first 5 ms.
- Phase 1 (`system/controlDict`): adaptive Δt (`maxCo 3.0`, `maxDeltaT 1e-4`), FW-H off, to 0.19 s.
- Phase 2 (`system/controlDict.phase2`): fixed Δt = 3.5·10⁻⁵ s for uniform FFT sampling, FW-H on.

## 5. Acoustic post-processing (FFT protocol)

Two 100 ms FW-H records, with the first 15 ms of each discarded, leave two 85 ms segments (486 points, fs = 5714 Hz) separated by a 15 ms gap. Each segment is detrended, Hann-windowed and zero-padded to 8192 points (true resolution 11.8 Hz), and the two power spectra are averaged. Tone levels are the peaks within ±6 Hz of 25.15 Hz and 125.7 Hz; the floor is the median level between 140 and 240 Hz.

I dropped two earlier choices. Concatenating the segments put a phase jump of about 136° into the 25 Hz component. A 5-point median filter raised the floor by 9 to 24 dB and hid the BPF harmonics (`figures/acoustics/raw_vs_filtered.png`). The June 2026 reports (`docs/archive/`) still use the concatenated, filtered signal.

OASPL is not reported. Absolute FW-H levels are inconsistent with a compact-dipole estimate from the blade forces (`paper/scripts/dipole_check.py`), so only tone frequencies and directivity are used as results.

## 6. Practical problems and fixes

| Symptom | Root cause | Fix |
|---|---|---|
| Solver crashes at step 1 (`sigFpe`) | about 50 k degenerate cells give 1/Vol = Inf/NaN | `export FOAM_SIGFPE=false` (tolerates them; does not fix quality) |
| Parallel `snappyHexMesh` aborts | load-balancing of the AMI faceZone across ranks | mesh serially (slower, reliable) |
| Time step collapses or diverges | one sliver cell gives a high local Courant number and an unbounded ω | adaptive Δt with `maxCo 3.0`, so the worst cell sets the step |
| Pressure solver gives NaN (coarse mesh) | a zero-determinant cell makes GAMG ill-conditioned | pressure solver switched from GAMG to PCG/DIC |
| Distributed run stalls | AMI + GAMG communication over Ethernet | run on a single node; use the cluster for throughput (one case per machine) |

## 7. HPC note

On a 0.78 M-cell test case, one node saturates at 4 cores (1.51× over 1 core), and 12 cores across both machines reach 1.93×. For the fan case, splitting one run across the two machines gave no gain, so I used the cluster to run independent cases side by side.
