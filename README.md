# Aeroacoustic CFD of a ducted fan with owl-inspired serrated blades

I simulated a 5-blade ducted desk fan at 1509 rpm with URANS in OpenFOAM (k-ω SST, sliding mesh) and coupled it to a Ffowcs Williams–Hawkings analogy to predict the tonal noise at 1 m. The run gives 43 mN of thrust, 108 m³/h, and a blade-passing tone at 126 Hz.

<p align="center">
  <img src="media/anim_pales_orbit.gif" width="640"><br>
  <em>Five serrated blades, colored by surface pressure (kinematic, m²/s²).</em>
</p>

Polytechnia is the engineering consultancy I co-founded in 2026. One of its first projects was a quiet desk fan whose blades carry trailing-edge serrations inspired by owl feathers (Wei et al., *Nature Communications*, 2024). Between late May and June 2026 I started from the Polytechnia FreeCAD model and built the simulation chain to see what such a fan sounds like before printing it: STL export and checks, mesh, rotating URANS, acoustic analogy, spectra. This was internal Polytechnia work; no client has received it. The structural side of the same fan (modal and drop-shock analysis of the 3D-printed housing) is in [ventilator-modal-analysis-petgcf](https://github.com/clementdumeril/ventilator-modal-analysis-petgcf).

The full report, in preprint form, is [`paper/main.pdf`](paper/main.pdf).

## Results

| Quantity (established regime, ω = 158 rad/s) | Value | Source |
|---|---:|---|
| Axial thrust / torque | 43.3 mN / 1.02 mN·m | [`perf_aero.py`](postprocessing/aero/perf_aero.py), [`data/forces_data`](data/forces_data) |
| Aerodynamic power (torque × ω) | 0.162 W | same |
| Volume flow rate, 60 mm downstream | 108 m³/h (30 L/s) | integrated on the downstream cut; needs the unpublished ParaView slices |
| Blade-passing tone (BPF = 5 × 25.15 Hz) | 126 Hz, 20 to 28 dB above the local floor on all four microphones, with harmonics at 252 and 378 Hz | [`post_acoustique_combine.py`](postprocessing/acoustics/post_acoustique_combine.py), [`data/`](data) |
| Strongest tone | 25 Hz (once per revolution), in the rotor plane only: more than 45 dB above its level on the axis, where it stays at the floor | same |

<p align="center">
  <img src="figures/flow/meridional_velocity.png" width="85%"><br>
  <em>Meridional cut: the blades drive an annular jet along the duct wall, with a slow recirculating wake behind the motor hub.</em>
</p>

<p align="center">
  <img src="figures/acoustics/spl_spectrum.png" width="85%"><br>
  <em>FW-H spectra relative to each microphone's floor: two 85 ms records (15 ms start-up removed from each), one Hann-windowed spectrum per record, averaged. The BPF and its second harmonic appear on every microphone, the third on all but the radial one; the 25 Hz tone only in the rotor plane. Shaded: a comb of lines spaced 100.5 Hz, probably numerical (see below).</em>
</p>

I did not expect a 25 Hz tone from a 5-blade rotor. The force history explains it: the in-plane blade force (0.95 mN, of which 0.77 mN from pressure) rotates with the rotor and oscillates at the shaft frequency, while the thrust fluctuates near the BPF. A rotating in-plane force radiates mostly in the rotor plane. Consistent with that, the 25 Hz tone is 3 dB weaker at 45° than at the radial microphone, with nothing above the floor on the axis. It has the same level and phase in both records, so it is a steady feature of the run.

At 25 Hz the rotor is acoustically compact, and a compact dipole rotating in the rotor plane radiates nothing on the axis and 3.0 dB less at 45° than in the plane. The FW-H data give 3.1 dB. I have not established whether that imbalance comes from the STL and mesh (the surfaces are faceted) or from the flow. The on-axis BPF also points to a loss of five-fold symmetry, since a perfectly symmetric rotor with three arms would not radiate the BPF on the axis either (the axisymmetric mode first appears at 3×BPF).

Above 600 Hz the spectra contain lines every 100.5 Hz (642 to 1446 Hz), up to 44 dB above the floor. They come from bursts in the blade force whose period, fitted on ten bursts, is 9.943 ± 0.011 ms, which is a quarter revolution and not the 10 ms field write interval. So the comb is locked to the rotor four times per revolution, probably through the four-fold symmetry of the box and Cartesian mesh seen across the AMI interface. I treat it as numerical and leave it out of the results ([`paper/scripts/comb_check.py`](paper/scripts/comb_check.py)).

<p align="center">
  <img src="figures/acoustics/directivity.png" width="70%"><br>
  <em>Absolute levels are unreliable (see the force-based check below); compare levels between microphones only.</em>
</p>

## How it works

Tonal fan noise comes from the unsteady pressure the blades exert on the air, which acts as a rotating dipole. I could not read it straight from the CFD: the domain is only 0.4 m wide and the solver is incompressible, so it carries no sound waves. FW-H gets around both problems. The pressure on the blade surface is recorded during the run, and the wave equation is integrated analytically to each microphone at its retarded time.

```
FreeCAD → STL → snappyHexMesh + cyclicAMI → pimpleFoam k-ω SST → FW-H (blade surface) → p'(t) at 4 mics → FFT
```

<p align="center">
  <img src="media/anim_ecoulement_fluide.gif" width="560"><br>
  <em>Q-criterion vortices shed by the rotating blades, colored by velocity: tip vortices and the helical wake.</em>
</p>

## What I built

- I exported the casing, rotor, hub and motor from the Polytechnia FreeCAD assembly to STL with Python scripts, and checked the tip clearance and duct bore on the exported surfaces.
- The mesh ([`snappyHexMeshDict`](openfoam_case/system/snappyHexMeshDict)) is hex-dominant snappyHexMesh: 10 mm background cells refined to 0.31 mm on the blades (level 5) so that the serrations are resolved, plus refined rotor and wake cylinders. The final mesh has 1.95 M cells, with max non-orthogonality 61. An earlier 0.91 M-cell mesh (blades at level 4) ran locally under WSL and is not used; its dictionary is kept as `snappyHexMeshDict.baseline_0p91M`.
- For the sliding interface ([`createBafflesDict`](openfoam_case/system/createBafflesDict)), a cylindrical faceZone (R = 72 mm) around the rotor is split into a cyclicAMI pair with `createBaffles`. About 35 k faces are interpolated at each step.
- The flow model ([`openfoam_case/`](openfoam_case)) is `pimpleFoam` with k-ω SST and solid-body rotation at 158 rad/s, after a 5 ms start-up ramp. Phase 1 ([`controlDict`](openfoam_case/system/controlDict)) runs to 0.19 s with an adaptive Δt (`maxCo 3.0`, `maxDeltaT 1e-4`) until the loads settle. Phase 2 ([`controlDict.phase2`](openfoam_case/system/controlDict.phase2)) restarts with a fixed Δt = 3.5·10⁻⁵ s so the FW-H signal is uniformly sampled.
- The acoustics ([`fwhControl`](openfoam_case/system/fwhControl)) use Farassat formulation 1A (libAcoustics) on the rotating blade surface, with four observers at 1 m: downstream, upstream, radial and 45°.
- Post-processing ([`postprocessing/`](postprocessing)) is Python for the spectra, loads and flow rate, and pvpython for the 3-D renders and animations.
- Everything ran on a 2-node Linux cluster I set up for this project and for the Code_Aster runs (8 + 4 cores, MPI over Ethernet, shared NFS case directory).

<p align="center">
  <img src="figures/mesh/flow_7_maillage_coupe.png" width="48%">
  <img src="figures/flow/extra_dentelures_zoom.png" width="40%"><br>
  <em>Left: mesh refinement from the far field to the blade band and AMI cylinder. Right: the serrated trailing edges as resolved by the mesh.</em>
</p>

<details>
<summary>How I got to this setup</summary>

I started at the end of May with steady MRF runs (3000, then 4320 rpm) to get flow rates quickly. For the transient run I first tried the non-conformal coupling (NCC) of OpenFOAM.org. Over six mesh variants (re-ordered interfaces, Delaunay caps, rounded caps) the runs kept stopping in `fvMeshStitcher` around 200 µs of simulated time, so on 12 June I moved to the ESI OpenFOAM and cyclicAMI. The other practical fixes are listed in [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md#6-practical-problems-and-fixes): serial snappyHexMesh (the parallel run aborted on the AMI zone) and PCG instead of GAMG for the pressure on the coarse mesh.

</details>

## How I checked it

The loads settle before FW-H recording starts: thrust varies by less than 1 % per revolution after about 4.75 revolutions.

<img src="figures/acoustics/blade_loads_time.png" width="80%">

`perf_aero.py` and `post_acoustique_combine.py` recompute thrust, torque and tone levels from the files in [`data/`](data), and the values in the Results table come from those runs. The flow rate needs the ParaView slices, which are not published.

The two FW-H records are 85 ms each once their start-up is removed, with a 15 ms gap between them. My first version concatenated them, which puts a 136° phase jump into the 25 Hz component, and applied a 5-point median filter against numerical spikes. Neither held up. The few outlier samples in the segments do not need filtering (the large values on the axial microphones are in the discarded start-up), and the filter raised the floor by 9 to 24 dB and hid the BPF harmonics. The spectra are now computed per record on the raw signal and then averaged. On the radial microphones the 25 Hz tone has the same level and phase in both records.

<img src="figures/acoustics/raw_vs_filtered.png" width="80%">

I also compared the FW-H output with the blade forces ([`paper/scripts/dipole_check.py`](paper/scripts/dipole_check.py)). At 0 Hz, 25 Hz and the BPF the rotor is acoustically compact, so the FW-H signal should match a dipole built from the integrated blade pressure force. The directivities match: the mean signal follows cos θ, and the 25 Hz tone is empty on the axis and 3 dB lower at 45°. The amplitudes do not (−9 dB for the mean load, +14 dB at 25 Hz, −18 dB for the BPF on the axis). The mean-load case has no time derivative and no retardation, so my first suspect is the input to the acoustic step: pressure scaling, normals or the library output convention. Until that is sorted out, I count the tone frequencies and directivities as results, but not the levels.

The `checkMesh` figures and the cells the solver has to tolerate are listed in [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md#3-mesh). I have not measured their effect on the spectrum.

<details>
<summary>Cluster benchmark</summary>

On a 0.78 M-cell test case, the gain on one node saturates at 4 cores (1.51× over 1 core, 1.44× on 8). 12 cores across both machines reached 1.93×, at 16 % parallel efficiency. For this AMI case, splitting one run across the two machines gave no gain, so I used the cluster to run independent cases side by side.

</details>

## Limitations

- There is no mesh-independence study. I kept one 1.95 M-cell mesh and did not validate it against coarser and finer meshes.
- The retained mesh has no prism layers on the blades (an earlier attempt collapsed on the thin serrated edges). Blade y⁺ averages about 4 and peaks at 445, so the surface loading is preliminary.
- I did not check time-step independence with a second Δt.
- URANS resolves tonal noise only. Serrations act mainly on broadband noise, and there is no smooth-blade reference run, so this study does not measure the effect of the serrations.
- The FW-H surface covers the blades only (no casing, arms, hub or quadrupole terms), and absolute levels are not compared with any measurement. They are also inconsistent with the blade forces (see above), a discrepancy I have not resolved.
- Two 85 ms records give an 11.8 Hz resolution. A contiguous 0.4 s record would take more than 100 h on the 1.95 M-cell mesh (16.5 s per step was measured on the 0.91 M-cell mesh), and the cluster holding the restart data was offline in October 2026.
- The solver only runs with floating-point trapping off (`FOAM_SIGFPE=false`) because of degenerate cells, and I did not measure their effect. Time integration is first order (227 steps per BPF period), and the lateral boundaries are 2.6 duct radii from the axis.

## Running it

The mesh and solver steps need the STL geometry (`constant/triSurface/{helice,carter,moteur,bras}.stl`), which is not published. Post-processing runs from the files in `data/`.

Post-processing, from the repository root (Python 3 with numpy, scipy, matplotlib; vtk for the flow figures; a few seconds each):

```bash
python postprocessing/acoustics/post_acoustique_combine.py      # spectra, tone levels, directivity, filter check
python postprocessing/aero/perf_aero.py                         # thrust, torque, power (+ flow rate if the slices are given)
python postprocessing/rendering/render_flow_2d.py <flow_vtk>    # needs meridien.vtp and avalZ.vtp exported from ParaView (not published)
```

Solver, with OpenFOAM (ESI; the published records come from the cluster runs, v2506) and [libAcoustics](https://github.com/unicfdlab/libAcoustics), once the STL files are in place:

```bash
cd openfoam_case
cp -r 0.orig 0
blockMesh && surfaceFeatureExtract
snappyHexMesh -overwrite          # serial: the parallel run aborts on the AMI faceZone
createBaffles -overwrite          # builds the cyclicAMI pair
decomposePar
export FOAM_SIGFPE=false          # the solver aborts at step 1 on degenerate cells otherwise (docs/METHODOLOGY.md, §6)
mpirun -np 8 pimpleFoam -parallel                    # phase 1, to 0.19 s
cp system/controlDict.phase2 system/controlDict
mpirun -np 8 pimpleFoam -parallel                    # phase 2, FW-H recording to 0.39 s
```

`paper/build_paper.sh` rebuilds the report ([`paper/main.pdf`](paper/main.pdf)). The June 2026 reports, which predate the re-analysis, are kept in [`docs/archive/`](docs/archive).

## Credits

- Serration concept and acoustic method: Z. Wei et al., *Towards silent and efficient flight by combining bioinspired owl feather serrations with cicada wing geometry*, Nature Communications, 2024 (G. X. Gu group, UC Berkeley).
- FW-H implementation: [libAcoustics](https://github.com/unicfdlab/libAcoustics). Solver: OpenFOAM (ESI). The published records come from the cluster runs, whose scripts source v2506; the earlier 0.91 M-cell run under WSL used v2406 and is not used. The cluster has been offline since June, so this could not be re-checked from the case directories.
- Fan geometry: Polytechnia FreeCAD model. Simulation setup, runs, post-processing and reports: Clément Dumeril.

[clement.dumeril.net](https://clement.dumeril.net) · [github.com/clementdumeril](https://github.com/clementdumeril) · MIT license
