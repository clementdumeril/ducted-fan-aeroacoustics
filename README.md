# Aeroacoustic CFD of a ducted fan with owl-inspired serrated blades

URANS simulation (OpenFOAM, k-ω SST, sliding mesh) of a 5-blade ducted desk fan at 1509 rpm, coupled to a Ffowcs Williams–Hawkings analogy to predict the tonal noise at 1 m: 43 mN of thrust, 108 m³/h, and a blade-passing tone at 126 Hz.

<p align="center">
  <img src="media/anim_pales_orbit.gif" width="640"><br>
  <em>Five serrated blades, coloured by surface pressure (kinematic, m²/s²).</em>
</p>

Polytechnia is the engineering consultancy I co-founded in 2026. One of its first projects was a quiet desk fan whose blades carry trailing-edge serrations inspired by owl feathers (Wei et al., *Nature Communications*, 2024). Between late May and June 2026, starting from the Polytechnia FreeCAD model, I built the simulation chain to see what such a fan sounds like before printing it: STL export and checks, mesh, rotating URANS, acoustic analogy, spectra. This was internal Polytechnia work; no client has received it. The structural side of the same fan (modal and drop-shock analysis of the 3D-printed housing) is in [ventilator-modal-analysis-petgcf](https://github.com/clementdumeril/ventilator-modal-analysis-petgcf).

Skills: OpenFOAM case set-up and debugging (sliding mesh, solver stability), signal processing (FFT, windowing, spike filtering checked against raw), Linux/MPI cluster built and benchmarked, Python/pvpython pipelines.

<!-- TODO Clément: origin of the 1509 rpm (158 rad/s) operating point; the modal study uses 2000 rpm -->

## Results

| Quantity (established regime, ω = 158 rad/s) | Value | Source |
|---|---:|---|
| Axial thrust / torque | 43.3 mN / 1.02 mN·m | [`perf_aero.py`](postprocessing/aero/perf_aero.py), [`data/forces_data`](data/forces_data) |
| Aerodynamic power (torque × ω) | 0.162 W | same |
| Volume flow rate, 60 mm downstream | 108 m³/h (30 L/s) | same, integrated on the downstream cut |
| Blade-passing tone (BPF = 5 × 25.15 Hz) | 126 Hz, 12 to 16 dB above the local floor on the axial microphones | [`post_acoustique_combine.py`](postprocessing/acoustics/post_acoustique_combine.py), [`data/`](data) |
| Strongest tone | 25 Hz (once per revolution), about 38 dB stronger in the rotor plane than on the axis | same |

<p align="center">
  <img src="figures/flow/meridional_velocity.png" width="85%"><br>
  <em>Meridional cut: the blades drive an annular jet along the duct wall, with a slow recirculating wake behind the motor hub.</em>
</p>

<p align="center">
  <img src="figures/acoustics/spl_spectrum.png" width="85%"><br>
  <em>Far-field SPL at the four microphones (170 ms record, Hann window, Δf = 5.9 Hz). The BPF tone is visible on every microphone; the 25 Hz tone dominates in the rotor plane.</em>
</p>

The 25 Hz tone was not expected for a 5-blade rotor. The force history explains it: the in-plane blade force (about 0.95 mN) rotates with the rotor and oscillates at the shaft frequency, while the thrust fluctuates near the BPF. A rotating in-plane force radiates mostly in the rotor plane, which matches the radial and 45° microphones being 35 to 38 dB louder at 25 Hz than the axial ones. I have not established whether that imbalance is physical (blade-to-blade geometry) or comes from the mesh.

<p align="center">
  <img src="figures/acoustics/directivity.png" width="70%"><br>
  <em>Absolute levels are underpredicted (incompressible URANS, blade-only FW-H surface); compare levels between microphones only.</em>
</p>

## What I built

- Geometry: exported the casing, rotor, hub and motor from the Polytechnia FreeCAD assembly to STL with Python scripts, and checked the tip clearance and duct bore on the exported surfaces.
- Mesh ([`snappyHexMeshDict`](openfoam_case/system/snappyHexMeshDict)): hex-dominant snappyHexMesh, 10 mm background cells refined to 0.31 mm on the blades (level 5) so that the serrations are resolved. Final mesh: 1.95 M cells, max non-orthogonality 61.
- Sliding interface ([`createBafflesDict`](openfoam_case/system/createBafflesDict)): a cylindrical faceZone (R = 72 mm) around the rotor, split into a cyclicAMI pair with `createBaffles`, about 35 k faces interpolated at each step.
- Flow model ([`openfoam_case/`](openfoam_case)): `pimpleFoam`, k-ω SST, solid-body rotation at 158 rad/s with a 5 ms start-up ramp. Phase 1 ([`controlDict`](openfoam_case/system/controlDict)) runs to 0.19 s with an adaptive Δt (`maxCo 3.0`, `maxDeltaT 1e-4`) until the loads settle; phase 2 ([`controlDict.phase2`](openfoam_case/system/controlDict.phase2)) restarts with a fixed Δt = 3.5·10⁻⁵ s so the FW-H signal is uniformly sampled.
- Acoustics ([`fwhControl`](openfoam_case/system/fwhControl)): Farassat formulation 1A (libAcoustics) on the rotating blade surface, four observers at 1 m: downstream, upstream, radial and 45°.
- Post-processing ([`postprocessing/`](postprocessing)): Python for the spectra, loads and flow rate; pvpython for the 3-D renders and animations.
- Compute: a 2-node Linux cluster I set up for this project and for the Code_Aster runs (8 + 4 cores, MPI over Ethernet, shared NFS case directory).

I started at the end of May with steady MRF runs (3000, then 4320 rpm) to get flow rates quickly. For the transient run I first tried the non-conformal coupling (NCC) of OpenFOAM.org. Over six mesh variants (re-ordered interfaces, Delaunay caps, rounded caps) the runs kept stopping in `fvMeshStitcher` around 200 µs of simulated time, so on 12 June I moved to the ESI OpenFOAM and cyclicAMI. The other practical fixes are listed in [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md#6-practical-pitfalls-solved-lessons-learned): serial snappyHexMesh (the parallel run aborted on the AMI zone) and PCG instead of GAMG for the pressure on the coarse mesh.

## How it works

Tonal fan noise comes from the unsteady pressure the blades exert on the air, a rotating dipole. The domain is only 0.3 m wide and the solver is incompressible, so it carries no sound waves. FW-H solves this: the pressure on the blade surface is recorded during the run, and the wave equation is integrated analytically to each microphone at its retarded time.

```
FreeCAD → STL → snappyHexMesh + cyclicAMI → pimpleFoam k-ω SST → FW-H (blade surface) → p'(t) at 4 mics → FFT
```

<p align="center">
  <img src="media/anim_ecoulement_fluide.gif" width="560"><br>
  <em>Q-criterion vortices shed by the rotating blades over one revolution, coloured by velocity: tip vortices and the helical wake.</em>
</p>

<p align="center">
  <img src="figures/mesh/flow_7_maillage_coupe.png" width="48%">
  <img src="figures/flow/extra_dentelures_zoom.png" width="40%"><br>
  <em>Left: mesh refinement from the far field to the blade band and AMI cylinder. Right: the serrated trailing edges as resolved by the mesh.</em>
</p>

## How I checked it

- Established regime. Thrust varies by less than 1 % per revolution after about 4.75 revolutions, before FW-H recording starts.

  <img src="figures/acoustics/blade_loads_time.png" width="80%">

- Results recomputed from `data/`: `perf_aero.py` and `post_acoustique_combine.py` recompute thrust, torque, flow rate and tone levels from the files in [`data/`](data). The values in the Results table come from those runs.
- Raw vs filtered spectra: the pressure signals carry isolated numerical spikes, which the processing removes with a 5-point median filter. The 25 Hz and BPF levels agree within about 3 dB between raw and filtered signals on every microphone (2 dB or less for the BPF). Above about 130 Hz the filter adds artifacts and removes the 2×BPF peak, so I only claim the two low-frequency tones. The overall level (OASPL) moves by up to 7 dB on the axial microphones depending on filtering, so I do not use it as a result.

  <img src="figures/acoustics/raw_vs_filtered.png" width="80%">

- Mesh quality: `checkMesh` figures and the cells the solver has to tolerate are listed in [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md#3-mesh); their effect on the spectrum has not been measured.
- Cluster benchmark. On a 0.78 M-cell test case, the gain on one node saturates at 4 cores (1.51× over 1 core, 1.44× on 8). 12 cores across both machines reached 1.93×, at 16 % parallel efficiency. For this AMI case, splitting one run across the two machines gave no gain, so I used the cluster to run independent cases side by side.

## Limitations

- No mesh-independence study: one 1.95 M-cell mesh was retained, not validated against coarser and finer meshes.
- No prism layers on the blades (they could not be inserted on the thin serrated edges). Blade y⁺ averages about 4 and peaks at 445, so the surface loading is preliminary.
- Time-step independence was not checked with a second Δt.
- URANS resolves tonal noise only. Serrations act mainly on broadband noise, and there is no smooth-blade reference run, so this study does not measure the effect of the serrations.
- The FW-H surface covers the blades only (no casing, arms, hub or quadrupole terms), and absolute levels are not compared with any measurement.

## Reproduce

The mesh and solver steps need the STL geometry (`constant/triSurface/{helice,carter,moteur,bras}.stl`), which is not published. Post-processing runs from the files in `data/`.

Post-processing, from the repository root (Python 3 with numpy, scipy, matplotlib; vtk for the flow figures; a few seconds each):

```bash
python postprocessing/acoustics/post_acoustique_combine.py      # spectra, tone levels, directivity, filter check
python postprocessing/aero/perf_aero.py                         # thrust, torque, power (+ flow rate if the slices are given)
python postprocessing/rendering/render_flow_2d.py <flow_vtk>    # needs meridien.vtp and avalZ.vtp exported from ParaView
```

Solver, with OpenFOAM v2406/v2506 (ESI) and [libAcoustics](https://github.com/unicfdlab/libAcoustics), once the STL files are in place:

```bash
cd openfoam_case
cp -r 0.orig 0
blockMesh && surfaceFeatureExtract
snappyHexMesh -overwrite          # serial: the parallel run aborts on the AMI faceZone
createBaffles -overwrite          # builds the cyclicAMI pair
decomposePar
mpirun -np 8 pimpleFoam -parallel                    # phase 1, to 0.19 s
cp system/controlDict.phase2 system/controlDict
mpirun -np 8 pimpleFoam -parallel                    # phase 2, FW-H recording to 0.39 s
```

The full report is in [`docs/report_EN.pdf`](docs/report_EN.pdf), with a French version alongside.

## Credits

- Serration concept and acoustic method: Z. Wei et al., *Towards silent and efficient flight by combining bioinspired owl feather serrations with cicada wing geometry*, Nature Communications, 2024 (G. X. Gu group, UC Berkeley).
- FW-H implementation: [libAcoustics](https://github.com/unicfdlab/libAcoustics). Solver: OpenFOAM (ESI).
- Fan geometry: Polytechnia FreeCAD model. Simulation setup, runs, post-processing and reports: Clément Dumeril.

[clement.dumeril.net](https://clement.dumeril.net) · [github.com/clementdumeril](https://github.com/clementdumeril) · MIT licence
