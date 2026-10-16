# Methodology

Notes that go with the report (`paper/main.pdf`): the physics, the numerical choices, and the
practical problems I ran into along the way.

## 1. Source of the noise

The tonal noise of a fan is dominated by the unsteady aerodynamic force the blades exert on the fluid,
a compact rotating dipole source. In the stationary frame, a loaded blade rotating at frequency *n*
radiates tones at the blade-passing frequency BPF = *Z·n* (here 5 × 25.15 = 126 Hz) and its harmonics.
To capture it, the unsteady blade loading has to be resolved, so I used a URANS simulation with true
rotation (sliding mesh) rather than a frozen-rotor approximation.

## 2. The Ffowcs Williams–Hawkings analogy

Two obstacles prevent reading the pressure at the microphones directly:

1. The domain is about 0.4 m across; the observers (1 m) are outside it, and meshing out to 1 m finely
   enough to propagate sound would cost millions of extra cells.
2. The solver is incompressible, so the solution contains no acoustic waves (infinite sound speed).

The FW-H analogy gets around both. The unsteady flow is recorded on a surface *S* enclosing the source,
and the acoustic wave equation is integrated analytically to any observer **x** at time *t*:

- thickness (monopole) term: blade volume displacement;
- loading (dipole) term: the unsteady surface pressure, dominant for tonal fan noise;
- quadrupole term: volume turbulence, secondary at low Mach and not integrated here.

The implementation is Farassat formulation 1A, the standard retarded-time formulation for moving
surfaces (`libAcoustics`). The source surface *S* is the blade skin (`helice` patch) with
`nonUniformSurfaceMotion` (rotating surface). Reference values: c₀ = 343 m/s, ρ∞ = 1.2 kg/m³,
p_ref = 20 µPa.

*S* covers the blades only, which is a known limitation. The casing, the three support arms, the hub and the
volumetric quadrupoles are not included.

## 3. Mesh

The mesh is hex-dominant `snappyHexMesh`. A Cartesian background box is refined by octree around the STL
surfaces, then snapped onto the geometry. This keeps the cell count low with targeted local refinement,
and it handles the sliding interface natively.

Refinement, from a 10 mm base cell: blades level 5 (0.31 mm), casing and arms level 3, hub level 2, plus a
refined rotor cylinder and wake region.

For the sliding interface (AMI), `snappyHexMesh` creates a `faceZone` on a virtual cylinder (R = 72 mm)
separating the rotating rotor zone from the static zone. `createBaffles` splits each internal face of
this zone into a coupled cyclicAMI pair (AMI1 / AMI2); `createPatch` cannot duplicate internal faces
this way. At each step the AMI interpolates fluxes across the sliding interface (about 35 k faces).

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

The concave and low-determinant cells sit in geometric transitions, which is common when a hex mesh
follows a curved geometry. The solver runs with `FOAM_SIGFPE=false` to tolerate them; this does not fix
their quality, and their effect on the spectrum has not been measured.

On y⁺: prism layers were not used in the retained mesh (`addLayers false`; an earlier attempt collapsed on
the thin serrated blades), so there are no wall layers. Blade y⁺ averages about 4 (buffer layer) and reaches 445 locally, so the surface loading
is preliminary.

## 4. Numerical model

- `pimpleFoam` (incompressible, transient); air ρ = 1.2 kg/m³, ν = 1.5·10⁻⁵ m²/s.
- k-ω SST (URANS), with automatic wall treatment.
- Solid-body rotation ω = 158 rad/s, ramped from 0 over the first 5 ms to avoid an impulsive start.
- Schemes: bounded first-order convection, PIMPLE pressure-velocity coupling with `correctPhi`
  (required for a moving mesh).

### Two time phases

- Phase 1 (establish the flow, `system/controlDict`): adaptive Δt with `maxCo 3.0` and
  `maxDeltaT 1e-4`, FW-H off, run to 0.19 s.
- Phase 2 (acoustics, `system/controlDict.phase2`): restart with a fixed Δt = 3.5·10⁻⁵ s, because an
  FFT needs uniform sampling, and FW-H on. The forces in `data/forces_data/force_0.19.dat` are sampled
  every 10 steps at exactly 3.5·10⁻⁴ s, which confirms the fixed step.

## 5. Acoustic post-processing (FFT protocol)

Two 100 ms FW-H records (sim time 0.19-0.29 s and 0.29-0.39 s). The first 15 ms of each are discarded
(FW-H start-up), which leaves two 85 ms segments (486 points, fs = 5714 Hz) separated by a 15 ms gap.
Each segment gets its mean and linear trend removed, a Hann window, and a zero-padded FFT (8192 points;
true resolution 11.8 Hz). The two power spectra are averaged. Single-sided, no A-weighting. Tone levels
are the peaks within ±6 Hz of 25.15 Hz and 125.7 Hz; the local floor is the median level between 140
and 240 Hz.

I dropped two earlier choices after checking them. Concatenating the segments was the first: the 15 ms
gap introduces a phase jump (about 136° at 25 Hz), which splits and lowers the 25 Hz peak. The second was
a 5-point median filter against numerical spikes. The few outlier samples in the segments do not need
filtering (up to 4 per segment beyond 5 robust standard deviations; the large values on the axial
microphones are in the discarded start-up), and the filter raises the floor by 9 to 24 dB and hides the BPF
harmonics (`figures/acoustics/raw_vs_filtered.png`).

The overall level (OASPL) is not reported: absolute levels are unreliable (incompressible URANS,
blade-only FW-H surface), so only level differences between microphones and tone-to-floor ratios are used.

The June 2026 reports (`docs/archive/`) predate this re-analysis: their acoustic levels come from the
concatenated, filtered signal. The absolute FW-H levels are also inconsistent with a compact-dipole estimate
built from the blade forces (`paper/scripts/dipole_check.py`), so only tone frequencies and directivity are
used as results.

Units: in `pimpleFoam`, `p` is a kinematic pressure (m²/s²); conversion to Pa uses p_Pa = ρ·p_kin.

## 6. Practical problems and fixes

| Symptom | Root cause | Fix |
|---|---|---|
| Solver crashes at step 1 (`sigFpe`) | about 50 k degenerate cells give 1/Vol = Inf/NaN | `export FOAM_SIGFPE=false` (tolerates them; does not fix quality) |
| Parallel `snappyHexMesh` aborts | load-balancing of the AMI faceZone across ranks | mesh serially (slower, reliable) |
| Time step collapses or diverges | one sliver cell gives a high local Courant number and an unbounded ω | adaptive Δt with `maxCo 3.0`, so the worst cell sets the step |
| Pressure solver gives NaN (coarse mesh) | a zero-determinant cell makes GAMG ill-conditioned | pressure solver switched from GAMG to PCG/DIC |
| Distributed run stalls | AMI + GAMG communication over Ethernet | run on a single node; use the cluster for throughput (one case per machine) |

## 7. HPC note

The 2-node cluster (8 + 4 cores, MPI over Ethernet) was benchmarked on a 0.78 M-cell `simpleFoam`
case: on one node the gain saturates at 4 cores (1.51× over 1 core, 1.44× on 8 cores), and 12 cores
across both machines reach 1.93×. For the AMI + GAMG fan case, distributing a single run across the two
machines gave no speed-up because network communication cancels the extra cores, so I used the
cluster to run independent cases in parallel instead.
