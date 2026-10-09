# Methodology

A technical companion to the full reports (`report_FR.pdf`, `report_EN.pdf`). It documents the
physics, the numerical choices, and the practical pitfalls solved along the way.

## 1. Source of the noise

The tonal noise of a fan is dominated by the **unsteady aerodynamic force the blades exert on the
fluid** — a compact, rotating **dipole** source. In the stationary frame, a loaded blade rotating at
frequency *n* radiates tones at the blade-passing frequency BPF = *Z·n* (here 5 × 25.15 = 126 Hz) and
its harmonics. Capturing this therefore requires resolving the **unsteady blade loading** — hence an
URANS simulation with *true rotation* (sliding mesh), not a frozen-rotor approximation.

## 2. The Ffowcs-Williams & Hawkings analogy

Two obstacles prevent reading the pressure at the microphones directly:

1. The domain is only ~0.3 m across; the observers (1 m) are **outside** it, and meshing out to 1 m
   finely enough to propagate sound would cost millions of extra cells.
2. The solver is **incompressible** → the solution contains **no acoustic waves** (infinite sound speed).

The **FW-H** analogy resolves both. The unsteady flow is recorded on a surface *S* enclosing the
source, and the acoustic wave equation is integrated *analytically* to any observer **x** at time *t*:

- **Thickness (monopole)** term — blade volume displacement,
- **Loading (dipole)** term — the unsteady surface pressure (dominant for tonal fan noise),
- **Quadrupole** term — volume turbulence (secondary at low Mach, not integrated here).

We use **Farassat Formulation 1A**, the standard retarded-time formulation for moving surfaces
(`libAcoustics`). The source surface *S* is the **blade skin** (`helice` patch) with
`nonUniformSurfaceMotion` (rotating surface). Reference: c₀ = 343 m/s, ρ∞ = 1.2 kg/m³, p_ref = 20 µPa.

**Assumed limitation:** *S* is limited to the blades — the casing, three support arms, hub and the
volumetric quadrupoles are not included.

## 3. Mesh

**Type — hex-dominant `snappyHexMesh`.** A Cartesian background box is refined by octree around the
STL surfaces, then snapped onto the geometry. Chosen for accuracy at a low cell count, targeted local
refinement, and native handling of the sliding interface.

**Refinement** (base cell 10 mm): blades level 5 (0.31 mm), casing/arms level 3, hub level 2, plus a
refined rotor cylinder and wake region.

**Sliding interface (AMI).** `snappyHexMesh` creates a `faceZone` on a virtual cylinder (R = 72 mm)
separating the rotating rotor zone from the static zone. `createBaffles` splits each internal face of
this zone into a **coupled cyclicAMI pair** (AMI1 master / AMI2 slave); `createPatch` cannot duplicate
internal faces this way. At each step the AMI interpolates fluxes across the sliding interface (≈ 35 k faces).

**Quality & y⁺.** On the retained mesh, max non-orthogonality ≈ 61 (good), but ~50 k "concave" cells
persist (intrinsic to hex meshing of a curved geometry — non-fatal). The intended Low-Re prismatic
layers could not be inserted on the thin, serrated blades, so **there are no wall layers**: the blade
y⁺ averages ≈ 4 (buffer layer) and reaches 445 locally → the surface loading is **preliminary**.

## 4. Numerical model

- `pimpleFoam` (incompressible, transient); air ρ = 1.2 kg/m³, ν = 1.5·10⁻⁵ m²/s.
- **k-ω SST** (URANS) — accurate near walls, robust in adverse pressure gradients; automatic wall treatment.
- Solid-body rotation ω = 158 rad/s with a **0 → 158 rad/s ramp over 5 ms** to avoid an impulsive-start shock.
- Schemes: first-order-bounded convection, PIMPLE pressure-velocity coupling with `correctPhi` (mandatory for moving mesh).

### Two-phase time strategy
- **Phase 1 (establish):** adaptive Δt (Co ≤ 5). In practice a fine cell caps Δt ≈ 3.7·10⁻⁵ s.
- **Phase 2 (acoustics):** **fixed** Δt = 3.5·10⁻⁵ s (an FFT requires *uniform* sampling), FW-H enabled.
  The fixed value is the one the adaptive controller already validated — so stability and cost are the same.

## 5. Acoustic post-processing (FFT protocol)

Combined 0.17 s record (two contiguous segments), **Hann** window, N = 972 points, fs = 5714 Hz
(Nyquist 2857 Hz), **single block** (no Welch averaging), single-sided spectrum, Δf ≈ 6 Hz, **no
A-weighting**. OASPL from the RMS signal, L = 20·log₁₀(p_rms / p_ref). The mean is removed per segment
(the raw signals carry a non-zero offset and must not be read as pure acoustic pressure); a 5-point
median filter removes numerical spikes — the raw and filtered spectra are compared for transparency.

**Units:** in `pimpleFoam`, `p` is a *kinematic* pressure (m²/s²); conversion to Pa uses p_Pa = ρ·p_kin.

## 6. Practical pitfalls solved (lessons learned)

| Symptom | Root cause | Fix |
|---|---|---|
| Solver crashes at step 1 (`sigFpe`) | ~50 k degenerate cells → 1/Vol → Inf/NaN | `export FOAM_SIGFPE=false` (tolerate; does not fix quality) |
| Parallel `snappyHexMesh` aborts | fragile load-balancing of the AMI faceZone across ranks | mesh **serially** (slower, reliable) |
| Time step collapses / diverges | one sliver cell → high local Courant → ω (turbulence) unbounded | Δt sized so the worst cell stays at Co ≈ 5 |
| Pressure solver → NaN (coarse mesh) | a determinant-0 cell makes GAMG ill-conditioned | switch p solver **GAMG → PCG/DIC** (robust) |
| Distributed run stalls | AMI + GAMG communication over Ethernet | run single-node; cluster is for *throughput* (1 case/machine), not one-case speed-up |

## 7. HPC note

A self-built 2-node cluster (12 cores, MPI over Ethernet) was benchmarked. For this AMI + GAMG case,
distributing a single simulation across nodes gave **no speed-up** (network communication cancels the
extra cores). The lesson — reported honestly — is that such a cluster is best used for **throughput**
(independent cases in parallel), not for accelerating a single communication-bound run.
