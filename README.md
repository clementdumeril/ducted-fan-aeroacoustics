# Aeroacoustic CFD of a ducted fan with owl-inspired serrated blades

I simulated a 5-blade ducted desk fan at 1509 rpm with URANS in OpenFOAM (k-ω SST, sliding mesh) and coupled it to a Ffowcs Williams–Hawkings analogy to predict the tonal noise at 1 m. The run gives 43 mN of thrust, 108 m³/h, and a blade-passing tone at 126 Hz.

<p align="center">
  <img src="media/anim_pales_orbit.gif" width="640"><br>
  <em>Five serrated blades, colored by surface pressure (kinematic, m²/s²).</em>
</p>

This was internal work at Polytechnia, the engineering consultancy I co-founded in 2026, on a quiet desk fan with trailing-edge serrations inspired by owl feathers (Wei et al., 2024). No client has received it. Between late May and June 2026 I built the chain from the FreeCAD model: STL export, mesh, rotating URANS, FW-H, spectra. The structural side of the same fan is in [ventilator-modal-analysis-petgcf](https://github.com/clementdumeril/ventilator-modal-analysis-petgcf). The full report is [`paper/main.pdf`](paper/main.pdf).

## Results

| Quantity (ω = 158 rad/s) | Value |
|---|---:|
| Axial thrust / torque | 43.3 mN / 1.02 mN·m |
| Volume flow rate, 60 mm downstream (net flux through the cut; needs ParaView slices not in the repo) | 108 m³/h |
| Blade-passing tone (BPF = 5 × 25.15 Hz) | 126 Hz, 20 to 28 dB above the local floor on all four microphones |
| Strongest off-axis tone | 25 Hz (once per revolution): 61 dB above the floor radially, 48 dB at 45°, at the floor on the axis |

<p align="center">
  <img src="figures/acoustics/spl_spectrum.png" width="85%"><br>
  <em>FW-H spectra relative to each microphone's floor. Shaded: a comb of lines spaced 100.5 Hz, probably numerical.</em>
</p>

The 25 Hz tone comes from the in-plane blade force, which rotates with the rotor. The comb above 600 Hz is locked to a quarter revolution, so I treat it as numerical ([`paper/scripts/comb_check.py`](paper/scripts/comb_check.py)).

## Setup

```
FreeCAD → STL → snappyHexMesh + cyclicAMI → pimpleFoam k-ω SST → FW-H (blade surface) → p'(t) at 4 mics → FFT
```

- Mesh: 1.95 M cells (checkMesh log not versioned), blades refined to 0.31 mm so the serrations are resolved, cyclicAMI sliding interface at R = 72 mm.
- Solver: `pimpleFoam` at 158 rad/s, first with an adaptive Δt to 0.19 s, then a fixed Δt so the FW-H signal is uniformly sampled.
- Acoustics: Farassat 1A (libAcoustics) on the blade surface, four microphones at 1 m.
- Runs: a 2-node Linux cluster I set up (8 + 4 cores, MPI over Ethernet).

The problems I hit and how I fixed them are in [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md#6-practical-problems-and-fixes).

## Limitations

- No mesh or time-step independence study, and no prism layers on the blades (y⁺ peaks at 445).
- URANS resolves tonal noise only, and there is no smooth-blade reference run, so this does not measure the effect of the serrations.
- Absolute levels disagree with a dipole estimate from the blade forces, so I only use tone frequencies and directivity. Two 85 ms records give an 11.8 Hz resolution, and the cluster is currently offline.

## Running it

The STL geometry is not published, so only post-processing runs from this repo (`pip install -r requirements.txt`, Python >= 3.10):

```bash
python postprocessing/acoustics/post_acoustique_combine.py   # spectra, tone levels, directivity
python postprocessing/aero/perf_aero.py                      # thrust, torque, power
```

With OpenFOAM ESI v2506 (used for the published records; the case dictionaries still carry v2406 headers from an earlier run), [libAcoustics](https://github.com/unicfdlab/libAcoustics) and the STL files in `constant/triSurface/`:

```bash
cd openfoam_case
cp -r 0.orig 0
blockMesh && surfaceFeatureExtract
snappyHexMesh -overwrite          # serial only
createBaffles -overwrite
decomposePar
export FOAM_SIGFPE=false
mpirun -np 8 pimpleFoam -parallel     # phase 1
cp system/controlDict.phase2 system/controlDict
mpirun -np 8 pimpleFoam -parallel     # phase 2, FW-H recording
```

The June 2026 reports, which predate the re-analysis, are in [`docs/archive/`](docs/archive).

## Credits

- Serration concept and acoustic method: Z. Wei et al., *Nature Communications*, 2024 (G. X. Gu group, UC Berkeley).
- FW-H: [libAcoustics](https://github.com/unicfdlab/libAcoustics). Solver: OpenFOAM (ESI).
- Fan geometry: Polytechnia FreeCAD model. Simulation setup, runs, post-processing and reports: Clément Dumeril.

[clement.dumeril.net](https://clement.dumeril.net) · [github.com/clementdumeril](https://github.com/clementdumeril) · MIT license
