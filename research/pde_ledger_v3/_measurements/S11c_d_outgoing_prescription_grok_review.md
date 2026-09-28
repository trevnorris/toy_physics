I'll review this packet as an independent method/implementation reviewer, starting from the packet index and the two files named first.The packet index is loaded. Next I’ll read the review prompt and the two named starting files.I’ll read the implementation, inputs, and source definitions next so the method review is grounded in the supplied artifacts.The implementation is in view. I’ll read the physics spec, inputs, and the kernel reports next to check the pole sign, residue formula, and coverage joins.I’ll now read the physics specs, kernel reports, and the remaining source excerpts that define the residue and radiation rules.Checking one remaining radiation-rule section and the kernel-report language on cofactors so the residue formula is pinned to the saved objects.# Independent method/implementation review

**Verdict: CLEAR FOR THIS BOUNDED STAGE**

The supplied worker is a correctly assembled **candidate** real-axis outgoing prescription on this fixed evanescent-bulk slice. It is not an accepted outgoing Green operator. Scope clearance here does not accept scientific output.

---

## Method question

The real-normal-momentum rule is justified for **separated** source and observation on this slice. It uses the inherited decaying bulk branch, the full signed-current 2×2 blocks, and the `exp(+i k_n (z-z_p))` reconstruction. It does not claim equivalence to a retarded `ω+i0` continuation; the worker stores that claim as false (`S11c_d_outgoing_prescription.py:333–334`).

The computed object remains a candidate. Two premises still sit outside what this stage may claim. Closing them does not require a new root census, frequency campaign, or A12 witness.

---

## What this stage actually builds

The worker source-joins the accepted 5×5 inverse to the saved REFERENCE modal/pairing packets, reuses whole degenerate blocks, builds exact inverse residues from determinant/cofactor derivatives, compares them to the saved projected blocks, and writes a symmetric principal-value Fourier integral plus current-oriented delta terms. Status is `FIXED_INPUT_OUTGOING_PRESCRIPTION_CANDIDATE_BUILT` with `outgoingGreenOperatorAccepted: False` and `resultMethodAcceptancePending: True` (`S11c_d_outgoing_prescription.py:329–334`). That matches the implementation directive (`S11c_d_outgoing_prescription_implementation.md:112–117`).

---

## Checks that hold

### Pole sign for `exp(+i k (z-z_p))`

`ConstantEndPencil.strong_matrix` uses field factors `exp(+i k_n z)` (`source-definitions.md:108–109`). `EdgeReduction` Fourier mass comes from `∫ du exp(-a u^2 + i q u)` (`source-definitions.md:42–46`), which evaluates to `2π`. Reconstruction is therefore `∫ dk_n /(2π) exp(+i k_n (z-z_p)) P(k_n)^{-1}`.

The worker uses that phase and the saved `fourierMass` (`S11c_d_outgoing_prescription.py:290–292, 302`). Sokhotski–Plemelj for denominator `k-k_a-i0 s_a` is `PV + i π s_a δ`. That is the stored delta coefficient (`:273, 302–303`).

Energy-current sign is the c1 radiation rule for real normal wavenumber (`S11c_c1_SHARED_PHYSICS.md:109–111`): keep waves that carry energy away from the source. `SIGNED_CURRENT = diag(sign(eigenvalues))` of the infinite-depth current (`source-definitions.md:760–776`). Positive diagonal means flux along increasing `z`. The worker sets `s_a = +1` on that block (`S11c_d_outgoing_prescription.py:258–260`). For `z>z_p` the contour encloses `s_a=+1`; for `z<z_p` it encloses `s_a=-1`. That is outgoing.

The local identity
`iπ(s + s_a) = 2πi s` when `s=s_a`, else `0`
(`:307–315`) is the Fourier residue theorem for this phase, including clockwise lower closure. It is an identity on the real-pole indentation. It is not a large-contour vanishing theorem.

The same current convention appears in `continuum_boundary.construct_end`: `SIGNED_CURRENT` times end orientation, outgoing when the product is positive (`S11c_d_continuum_boundary.py:199–202`). For the infinite-line Green function both current signs are kept, each attached to its own half-line. That is the right split.

### Multiplicity-two exact residue

If `det = (k-k_a)^m det^{(m)}/m! + ⋯` and `adj = (k-k_a)^{m-1} adj^{(m-1)}/(m-1)! + ⋯`, the inverse has a **simple** pole with residue `m · adj^{(m-1)}(k_a)/det^{(m)}(k_a)`. Nullity two does not make a double inverse pole when `L† P' R` is invertible.

The worker checks that (`:232–234`):
- `det^{(r)}(k_a)=0` for `r<m`, `det^{(m)}≠0`
- `adj^{(r)}(k_a)=0` for `r<m-1`
- residue `= m · adj^{(m-1)} / det^{(m)}`

The independent block formula is `R (L† P_{k} R)^{-1} L†` with an explicit 1×1 or 2×2 inverse (`:251–254`). That is the same Laurent ansatz as `EndResolventAudit.pole` (`source-definitions.md:999–1002`). Projector and both Laurent residuals are stored (`:261–267`). A singular projected derivative stops the job (`:251`).

Saved cofactors are stacked as a matrix and used as the adjugate (`:176–179, 234`). The accepted kernel report calls them transposed cofactors (`S11c_d_reference_kernel_report.md:12–14`). A transpose mismatch would move `cofactorVsBlockResidue` and fail the `1e-8` join (`:282`). That check is the convention lock.

### Coordinate, normalization, and unit joins

Exact checks:
- pairing pencil equals modal `PENCIL_PLUS` (`:129`)
- full 5×5 identity `P_modal(k_n, q_out) = P_kernel` (`:162–164`)
- `q_out = i sqrt(-q_phys^2)` on the strictly evanescent quadratic (`:148–158`)
- `q_out(k_a) = RADICAL_SCALE · Q_native = PHYSICAL_Q` (`:237–241, 263`)
- `k_a`, `ω` joins (`:262–264`)
- field/row units from `s11cdCurrentPlus*` and `s11cdPairingPlusRowResidual*` versus kernel units (`:165–173`)
- normal/bulk/spectral measure units (`:157`)
- flux-frequency identity `L_flux† (∂_ω P) R_flux = I` (`:268`)
- whole-block current `s_a J = I` (`:269`)
- one-sided current reversal moves the delta coefficient and breaks that identity, with residue held fixed (`:272–283`)

`eta_bg` and `σ_W` are set to zero for the reference slice (`:135–136`). The physical JSON still carries `eta_bg=1/100`. Correctness of that override is the exact 5×5 symbol join, which fails closed on mismatch.

### Real-pole coverage transfer

Saved inventory (`S11c_d_end_normalization_reference_thickness_repair_checkpoint.json:55–61, 292–351`):
- 9 disks, 18 unique lifts, all reality labels `PROVED_REAL` or `PROVED_NONREAL`
- 2 real disks / 4 real candidates
- 2 current-normalized blocks, indices 16 and 17, `k ≈ ±0.785`, opposite current eigenvalues, both on-sheet, both nullity two
- the other real pair (disk 7, indices 14–15) is off-sheet

The worker requires finite polynomial coverage, zero certificate residuals, denominator exclusion, both lifts, and then selects **on-sheet and real**, demanding current normalization and exact `K` (`:185–212`). That selection is the physical-sheet real poles of `P(k_n, q_out(k_n))` once the 5×5 identity holds: `q_out` lies on the acoustic curve, so a real zero of `det(P_kernel)` is a modal root with `q = q_out`. Off-sheet real `k` belong to the opposite radical and stay out of the PV exclusions.

No new census is required. The worker fails closed if a selected real on-sheet block lacks current normalization or an exact low-degree lift.

### Distribution form

Symmetric exclusions around the exact real poles, `ε → 0+` in momentum units, plus `exp(i k_a x) i π s_a A_a / fourier_mass` (`:293–303`). The Abel profile regulator is stored separately (`:323`). The regular integrand is kept; this is not a discrete mode sum. Complex off-axis poles stay inside the real-line values of the inverse.

---

## Why it remains a candidate, and what closes it

**1. Coincident-position contact terms and large-`k` class**

The contour identity is written only for `side ∈ {−1,+1}` and explicitly drops `x=0` (`:305–306`). The candidate pickle lists that domain as outstanding (`:325–327`). After PV regularization of simple poles, the Fourier integral is a tempered distribution only given the actual large-`k_n` growth of the bound inverse, which still contains `q_out ∼ |k_n|`.

Evidence that would close this, using only the already bound inverse: extract the polynomial/large-`k` part of `fixed_inverse`, state the decay order, and emit any `δ(z-z_p)` or derivative contact terms in the distributional Fourier transform at coincident arguments. That is a bounded exact operation on the saved inverse.

Until that is done, the candidate is a well-defined **separated-point** kernel (`z ≠ z_p`) plus a formal PV limit. It is not a complete Green operator on the diagonal.

**2. Off-sheet real-`k` invertibility is implied, not executed**

Selected poles have `q_out = PHYSICAL_Q` and `det^{(m)} ≠ 0` with lower derivatives zero. Off-sheet real records are excluded by sheet membership only. The 5×5 identity plus complete `(k,q)` coverage makes extra physical-sheet real zeros impossible, and the checkpoint’s two real disks match that split. A direct `fixed_det(k) ≠ 0` at those off-sheet real `k` would make the transfer an executed identity on the inverse being integrated. It is a small strengthening, not a missing census.

---

## Substantive blockers

None for this bounded stage. The worker does not accept the Green operator, does not claim `ω+i0` equivalence, and does not treat algebraic joins as method clearance.

---

## Optional observations

- `momentumDeltaCorrection` is the *k*-space addend `iπ s A δ(k-k_a)` without `/fourier_mass`; the *z*-space `singularKernelCorrection` includes that factor (`:301–303`). `candidateKernel = pv + correction` is internally consistent.
- Regular density is rebuilt from inverse, phase, and saved mass. The accepted `regular-spectral-density` pickle is not hash-joined. A join would be extra evidence, not a different object.
- Python chained comparison `int(m)==int(n) in (1,2)` means `m==n` and `n∈{1,2}` (`:211`). That is the intended multiplicity check.
- Launch containment (900 s, 2 GiB, zero swap, one CPU, nice 15, 32 tasks, one thread), durable journals before guards, and `SavedCodec` restricted to sympy/numpy/builtins match the gate (`:51–65, 85–94, 341–360`; `scripts/s11c_guarded_run.py:30–47, 139–185`).
- Hash pins in `S11c_d_outgoing_prescription_inputs.json` match the modal checkpoint, pairing relocation (`complete.pickle` `9838c505…`), and kernel inverse/cofactor artifacts in `saved-evidence.json`.

---

The bounded construction may proceed under the existing guard. Full outgoing-Green, FORM, two-asymptote response, A11, and A12 remain unclaimed, as the worker already records.
