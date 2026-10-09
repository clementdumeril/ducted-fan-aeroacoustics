# Methodology

A technical companion to the full reports (`report_FR.pdf`, `report_EN.pdf`): the physics, the numerical
choices, and the practical problems solved along the way.

## 1. Source of the noise

The tonal noise of a fan is dominated by the unsteady aerodynamic force the blades exert on the fluid,
a compact rotating dipole source. In the stationary frame, a loaded blade rotating at frequency *n*
radiates tones at the blade-passing frequency BPF = *Z·n* (here 5 × 25.15 = 126 Hz) and its harmonics.
Capturing this requires resolving the unsteady blade loading, hence a URANS simulation with true
rotation (sliding mesh) rather than a frozen-rotor approximation.

## 2. The Ffowcs Williams–Hawkings analogy

Two obstacles prevent reading the pressure at the microphones directly:

1. The domain is about 0.3 m across; the observers (1 m) are outside it, and meshing out to 1 m finely
   enough to propagate sound would cost millions of extra cells.
2. The solver is incompressible, so the solution contains no acoustic waves (infinite sound speed).

The FW-H analogy resolves both. The unsteady flow is recorded on a surface *S* enclosing the source,
and the acoustic wave equation is integrated analytically to any observer **x** at time *t*:

- thickness (monopole) term: blade volume displacement;
- loading (dipole) term: the unsteady surface pressure, dominant for tonal fan noise;
- quadrupole term: volume turbulence, secondary at low Mach and not integrated here.

The implementation is Farassat formulation 1A, the standard retarded-time formulation for moving
surfaces (`libAcoustics`). The source surface *S* is the blade skin (`helice` patch) with
`nonUniformSurfaceMotion` (rotating surface). Reference values: c₀ = 343 m/s, ρ∞ = 1.2 kg/m³,
p_ref = 20 µPa.

Assumed limitation: *S* covers the blades only. The casing, the three support arms, the hub and the
volumetric quadrupoles are not included.

## 3. Mesh

Type: hex-dominant `snappyHexMesh`. A Cartesian background box is refined by octree around the STL
surfaces, then snapped onto the geometry. It gives a low cell count with targeted local refinement,
and it handles the sliding interface natively.

Refinement (base cell 10 mm): blades level 5 (0.31 mm), casing and arms level 3, hub level 2, plus a
refined rotor cylinder and wake region.

Sliding interface (AMI): `snappyHexMesh` creates a `faceZone` on a virtual cylinder (R = 72 mm)
separating the rotating rotor zone from the static zone. `createBaffles` splits each internal face of
this zone into a coupled cyclicAMI pair (AMI1 / AMI2); `createPatch` cannot duplicate internal faces
this way. At each step the AMI interpolates fluxes across the sliding interface (about 35 k faces).

Quality (`checkMesh` on the retained mesh):

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

y⁺: the intended low-Re prism layers could not be inserted on the thin serrated blades, so there are
no wall layers. Blade y⁺ averages about 4 (buffer layer) and reaches 445 locally, so the surface loading
is preliminary.

## 4. Numerical model

- `pimpleFoam` (incompressible, transient); air ρ = 1.2 kg/m³, ν = 1.5·10⁻⁵ m²/s.
- k-ω SST (URANS), with automatic wall treatment.
- Solid-body rotation ω = 158 rad/s, ramped from 0 over the first 5 ms to avoid an impulsive start.
- Schemes: bounded first-order convection, PIMPLE pressure-velocity coupling with `correctPhi`
  (required for a moving mesh).

### Two-phase time strategy

- Phase 1 (establish the flow, `system/controlDict`): adaptive Δt with `maxCo 3.0` and
  `maxDeltaT 1e-4`, FW-H off, run to 0.19 s.
- Phase 2 (acoustics, `system/controlDict.phase2`): restart with a fixed Δt = 3.5·10⁻⁵ s, because an
  FFT needs uniform sampling, and FW-H on. The forces in `data/forces_data/force_0.19.dat` are sampled
  every 10 steps at exactly 3.5·10⁻⁴ s, which confirms the fixed step.

## 5. Acoustic post-processing (FFT protocol)

Two 100 ms FW-H records; the first 15 ms of each are discarded and the rest concatenated into a 0.17 s
signal (N = 972 points, fs = 5714 Hz). Hann window, single block (no Welch averaging), single-sided
spectrum, Δf ≈ 6 Hz, no A-weighting. The mean is removed per segment, because the raw signals carry a
non-zero offset. A 5-point median filter removes isolated numerical spikes; the raw and filtered
spectra are compared in `figures/acoustics/raw_vs_filtered.png`. The 25 Hz and BPF tones agree within
about 3 dB; above 130 Hz the filter adds artifacts, and the OASPL depends on the filter by up to 7 dB,
so neither is used as a result.

Units: in `pimpleFoam`, `p` is a kinematic pressure (m²/s²); conversion to Pa uses p_Pa = ρ·p_kin.

## 6. Practical pitfalls solved (lessons learned)

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
machines gave no speed-up because network communication cancels the extra cores. The cluster was
therefore used to run independent cases in parallel.
