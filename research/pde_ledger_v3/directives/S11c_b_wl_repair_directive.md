# S11c-b Wolfram engine — repair directive (round 1)

## Task and authority
Patch `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (committed `0c0e9a8a`)
to fix three defects the two build legs found by FORM ablation. ⭐ **Fix only these; keep everything that
passed unchanged** (the §3a energy-basis spurion construction, the operator, the N15 new invariants, the
operator-extraction data dependence — all ablation-clean). Its only write is the flushed stdout tag stream.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (HEAD `0c0e9a8a`, with the
round-2/round-3-repaired §3c/§3d) is the sole physics authority, together with the two sibling specs it
inherits (`S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`). The defects and their ablation evidence are
in `directives/_measurements/S11c_b_wl_build_review.md` (W1–W3); read it. Add no expected value (rule 5).

⛔⛔ **This is the blind engine.** It imports nothing and re-derives every object (including the S11c-a
substrate) from the specs. The build inputs besides the specs carry no S11c-b computed physics. ⛔ Do not
import `S11c_a_exports.py` or any transcript; emit no import tag; no do-not-read list (rule 12). Keep the
existing blindness (the acceptance test copies the `.wl` alone into an empty dir and runs it).

## The three fixes (implement the corrected spec sections; the finding record is the evidence)

**W1 — Admissibility operand is identically zero (BLOCKING).** `backgroundBalanceFromModel` (L1121–1143) scales
the perturbation (ε²) energy fields by `backgroundOrder`, takes `D[·,backgroundOrder] /. backgroundOrder→0`,
then `/. zeroWaveRules` — the ε→0 limit of the bilinear §3b operator, ≡ 0 (all four branch/density cases).
Rebuild it per the corrected **§3d**: the **background-order (ε⁰) first variation of the §3a
variable-coefficient energy-and-geometry functional written in FULL fields** (`W=W_bg+δW`,
`ρ_4D=ρ_4D,bg⁰(1+θ)`, brane displacement), with the thickness/coefficient **gradient content of the full
fields** (`∇(W_bg+δW)`, …) so the profile's own gradients are present at background order — ⛔ `|∇(δW)|²`
(perturbation-only gradient) is not an acceptable background-order representative. Vary wrt the brane
configuration, THEN evaluate at `𝔅⁰` (`θ=e_W=δW=0`, `u=0`; profile+jets retained). ⭐ Retain every spatial
derivative the variation/divergence generates at the §3a background-bookkeeper order (a second spatial
derivative of `W_bg` is still first order in background amplitude — not dropped). ⛔ Not `D[perturbation
energy(scaled fields), backgroundOrder]`; no `W_bg−W_0` insertion. Emit operator/support/residual in the same
ordered pairing as `𝒮_hold⁰`.

**W2 — Coupling-kernel off-diagonal blocks do not reduce to the decoupled zero (§3c).** `extractCoupling`
(L897–925) applies the sector split (`transverseInputRules`, `localSectorRules`) only to UNDIFFERENTIATED
field occurrences — inert on the operator's gradient content — so the uniform-limit block is a nonzero
diagonal-thickness expression. Rebuild extraction per the corrected **§3c**: a **weak variational restriction**
of the operator under the §1c pairing, using independent transverse (solenoidal, `∇·u_T=0`) and longitudinal
(irrotational, `∇×u_L=0`) trial and test displacements plus θ/e_W test fields — transverse→thickness is the
transverse-trial × {θ,e_W,u_L}-test pairing, reciprocal the other. This attributes both gradient and
undifferentiated-`u` spurion (`g·u`) terms without a global projector (N5), and must reduce to §1d's decoupled
zero in the uniform limit (a solenoidal `u_T` does not enter the uniform `∇·u` coupling). Take trial/test
fields of compact in-plane support so the IBP boundary term is zero. ⛔ Do not use value-only `field→0`
zeroing; ⛔ extract from the operator, not a parallel route.

**W3 — Adjointness residual tautological (§3c).** `ADJOINTNESS_RESIDUAL = D[D[E,transverseProbe],thicknessProbe]
− D[D[E,thicknessProbe],transverseProbe]` (L945–957) ≡ 0 by Clairaut for any energy. Replace with the
pairing-based residual between the two extracted operator blocks (per §3c) — whether the thickness→transverse
block is the adjoint of the transverse→thickness block under the §1c pairing. ⛔ Not the scalar-energy mixed
second derivative. If adjoint by construction, emit the two blocks and state there is no independent second
route rather than a structural zero.

## Unchanged from the build directive
Tag grammar `WL_S11CB_<QUANTITY>`; internal symbols without `_` but emitted-name strings exactly the spec's
`S11CB_*` bases; no `VERDICT`/`PASS`/`FAIL`; booleans as `Inactive[Equal]` retaining operands. Run discipline:
one kernel at a time, `timeout 600`, 6 GB kill, `--sandbox danger-full-access`; demonstration runs under
scratch, never `mathematica/out/`. Acceptance: blindness (copy `.wl` alone, byte-identical runs, no stray
write), flush under redirected stdout. Do not restate; they are the build directive's, unchanged.

## Report (≤15 lines)
Files written, the functions changed, tasks run, runtime, all emitted tag names, any non-computable object,
the blindness/flush demonstrations. Confirm the untouched-passing objects were not regressed. No computed value.
