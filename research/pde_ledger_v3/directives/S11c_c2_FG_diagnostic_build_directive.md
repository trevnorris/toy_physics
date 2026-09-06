# Build directive — S11c-c2 F/G re-grounding DIAGNOSTIC

**This directive is orchestrator-written; the DELIVERABLE script will be astra/Codex-written** (so its two legs are
a fresh Claude agent + Grok). Working dir = repo root `/var/projects/toy_physics`; all `directives/…`, `scripts/…`,
`_measurements/…` paths are under `research/pde_ledger_v3/`.

**Role.** F (uniform-limit decoupling) and G (directionality / adjointness) of the c2 self-energy increment were
resolved in the STEP A adjudication with **orchestrator-authored** scripts (`verify_F.py`, `verify_EG.py`), now
retired as biased. A Codex-sol question-vet + a 2-leg decision review sharpened both. Build ONE diagnostic script
that computes the objects below and **PRINTS operands + residuals**; the orchestrator adjudicates. ⛔ This script
decides nothing.

**Deliverable.** `scripts/S11c_c2_FG_diagnostic_sympy.py`. Import and REUSE the verified construction of
`scripts/S11c_c2_selfenergy_fold_sympy_audit.py`: `bind_inputs`/`load_model`, `build_case`, `extract`,
`retained_shape`, `shape_coefficients`, `difference`, `wave_jet`, and — for `ORDERING_COMMUTATOR` — `regression_coordinates`.
⛔ Do NOT re-implement the construction; ⛔ do NOT call `run()`. `ROOT` resolves from the script's own location.
Env `S11CC2_FG=TRIAGE` runs one `(α,ρ)` case fast. There is **no `close()` function**: the §3a closure is the
per-face substitution map returned by `build_case` (`model['substitutions']`, per-face slot/normal-jet maps at
`audit.py:509–523,565–567`); pin to that. ⚠ `model['faces'][s]` (`audit.py:549`) is the **per-face closed carrier**
`S_{P,s}[χ_s] = (rows − rows|_{P_s=0}).subs(χ_s)` (that face's slots only) — it is the **pre-`extract` carrier** (so
`I_closed,s = extract(model['faces'][s])`), ⛔ NOT the per-face increment `I_s`. Use the PINNED per-face set below,
⛔ do not treat `model['faces']` as an increment.

## ⭐⭐⭐ THE THREE SCRIPT CLAUSES — verbatim, non-negotiable
> **1. The script PRINTS computed CAS objects (expressions, roots, symbolic booleans, machine-readable tables). It
> may NOT state conclusions or emit a prose payload** — no "decoupled"/"directional"/"adjoint"/"no route" tags.
> **2. PRINT the residual; do NOT assert it.** Compute → `emit` → *then* optionally a guard that hard-stops only on
> a **structural** impossibility (a temp Dummy escaping, a shape mismatch), ⛔ NEVER `assert residual==0`.
> **3. Interpretation belongs to the orchestrator; the builder REPORT (prose) is separate from the script.**

⭐ The ONLY hand-combination of physical symbols is in the imported operands + the trial/test ansatz. Every other
expression is REACHED BY COMPUTATION; every control re-enters at the imported operands, ⛔ never at an emitted result.

---

## Diagnostic F — uniform-limit decoupling as a WEAK-KERNEL identity (⛔ not `.doit()==0`)

**Objects.** With `P = {δp_s±, ∂_w δp_s±}` (BOTH the slot and its w-jet — §3a closes both, `SHARED_PHYSICS:170–176`),
`S_P = SLAB − SLAB|_{P=0}` (valid because the slots enter linearly), `E = extract`, `χ` the per-face closure
substitution map, per `(α,ρ)`:
```text
I_closed = E( S_P[χ] )                       (closed image of the complete pressure-slot carrier) ,
I_bare   = − E( S_P )                          (the §3c open-slot bookkeeping; free δp_s AND ∂_w δp_s slots) .
```
**F1 — independent reconstruction (⛔ not a tautology).** Obtain the increment by an INDEPENDENT route and difference:
```text
I_direct = difference( build_case(…)[closed_kernel], build_case(…)[open_kernel] )   (the §3c increment) ,
R_split  = I_direct − (I_closed + I_bare) .
```
`emit R_split`. ⚠ This is a **transcription/linearity identity**, ⛔ NOT an F finding and ⛔ NOT an acceptance target.

**Per-face objects — PINNED (⛔ not builder-invented); simultaneous substitution + the same `retained_shape`/physical-field
stage as `build_case`:**
```text
P_s        = {δp_s, ∂_w δp_s}                 (that face's slots) ,
S_{P,s}    = SLAB − SLAB|_{P_s=0} ,
I_closed,s = extract( S_{P,s}[χ_s] )           ( = extract(model['faces'][s]) ) ,
I_bare,s   = − extract( S_{P,s} ) ,
I_s        = I_closed,s + I_bare,s ,
R_add      = I_direct − Σ_s I_s                (per-face additivity residual — emit it) .
```
F4 and G2 build their per-face objects from these; the assembled objects are the sums over `s∈{+,−}`.

**F2 — the complete uniform map, ENUMERATED, applied at the imported operands BEFORE close/extract.** ⛔ Do NOT reuse
the audit's `UNIFORM_LIMIT` override (`audit.py:1081` = the retired 3-symbol proxy `{η:0,σ_W:0,W_bg:W_0}`, which
leaves `W_bg_d*` live — an M3 freeze). Build the simultaneous substitution:
- `W_bg→W_0`, `mu_R_bg→mu_R`, `e_W_bg→e_W`;
- **every** derived spatial jet `W_bg_d*`, `mu_R_bg_d*`, and the `w1_profile`/`m1_profile` jets the ansatz identifies
  with them → 0;
- both density representatives via `background_density_map`: `rho_4D,bg → rho_br/W_0`, `rho_br,bg → rho_br`, and the
  background-density gradient → 0;
- `η=0`, `σ_W=0`, and the η-tied first-jet face tilt in `n̂_s` → 0.
⛔ This map is NOT `Z→0`, NOT `Λ_I→0`, NOT `μ_θ→0` (those are the OTHER §5e reduction limits — ⛔ do not conflate).
⚠ **Pass this complete map so it reaches BOTH the slab rows AND `build_face`/`kernel_bridge`** (via `overrides=` and/or a
specialized `Inputs`) — ⛔ `input_map` **alone** is forbidden for F: it uniformizes only the slab rows, leaving
`W_bg_d*`/profile hats live inside `χ`, re-creating the M3 freeze (`audit.py:541–546`). Include the c1 Fourier hats
`s11cc1_w1_profile_hat_transfer`/`s11cc1_w1_profile_jet_hat_*` (identified with `w1`) → their uniform values.
Emit the applied map's key/value inventory so the specialization is auditable.

**F3 — the weak-kernel predicate.** For `I_closed|_uniform`: extract every retained multigrade coefficient
(`shape_coefficients`), align integral bound dummies, and read each `Integral`'s **kernel** (`.function`) **AND any
non-`Integral` remainder** as weak operators against **generic** `wave_jet`/`extract` trial/test potentials (⛔ do NOT
specialize them). ⚠ A true weak-zero can collapse to a **no-`Integral`** form (`audit.py:472–473` sends a zero
integrand to `sp.S.Zero`), so the remainder MUST be inspected — an absent `Integral` is ⛔ NOT itself weak-zero. The only
permitted integration-by-parts is the S11c-b §3c in-plane interior IBP with the boundary term fixed to 0
(`S11c_b_SHARED_PHYSICS.md:328–330`); ⛔ NO IBP or `.doit()` of the `Z`/Fourier/distributional integrals
(`audit.py:453–469`). "Weak-zero" = the kernel coefficients on generic potentials vanish — ⛔ never `.doit()==0`, an
absent `Integral`, or zero after Fourier evaluation.

**F4 — emit, per face AND assembled.** `emit` `I_closed|_uniform`, its per-grade weak-coefficient table, and the
weak-zero checker result, **each per-face `s∈{+,−}` and assembled** (so a face-parity cancellation is visible, not
mislabeled structural). `emit` the full `I_bare` (both slot families). ⛔ Assembled-only is a forbidden proxy.

**F5 — slot-inventory control.** `emit` the free-symbol inventory: `I_closed` must contain **none** of `P`; `I_bare`
must show the `δp_s`/`∂_w δp_s` slots where present. (Catches a missed `∂_w δp_s` substitution.)

**F6 — live sentinel.** Inject a source-level transverse↔thickness coupling **into the F-checked path** (the
pressure-slot carrier `S_P` or the map `χ`, so it reaches `I_closed|_uniform`), with a **polynomial test×trial
kernel carrying no `Integral`** (so its nonzero-ness does not depend on `.doit()`), run the SAME weak-zero checker,
and `emit` the checker's **result** on the injected kernel (baseline vs injected). ⛔ Do NOT assert "returns nonzero" —
emit the result, adjudicated on our side. Tag it as injected control data; ⛔ never use it as the F finding.

**F7 — mandatory FORM ablation.** Relocate a `P` slot **across sectors** to a **dimension/retained-grade-matched**
destination row, using **named existing** source + destination operands (e.g. a θ-row `δp_s` coefficient into a U-row),
AND flip a sign — a genuine FORM change of a load-bearing source object, ⛔ NOT a ±-face swap (a no-op at a true
uniform zero), ⛔ NOT a rescale. `emit` the named FROM/TO, baseline operands, mutated operands, and their literal difference. ⛔ Do NOT prescribe that it
"must change" — the orchestrator adjudicates. A COEFFICIENT rescale tests arithmetic only.

**F8 — scope note (state it in the emitted metadata).** Even a clean weak-zero is only a **secondary smoke test**:
it cannot validate the nonuniform coefficient, sign, or parity (decisions `N6` — the uniform limit is a smoke test only). ⛔ No target-value citation may be supplied.

**Separately:** `emit` the uniform limit of `ORDERING_COMMUTATOR` via the audit's path (close the imported
`coupling_kernel`, run `regression_coordinates` — the audit uses a function ansatz vs the predecessor's jet symbols,
so the commutator is a coordinate-mismatch object without it, `audit.py:1055–1062,1174–1205`). ⛔ Do NOT conflate
`ORDERING_COMMUTATOR`, `I_closed`, or `CLOSED_COUPLING_KERNEL` — they are three distinct objects.

⛔ **Forbidden F proxies:** `.doit()==0` / absent `Integral` / Fourier-evaluated zero; special/frozen trial field,
`k=k′`, odd-integrand symmetry, special boundary; dropping `∂_w δp_s`; uniformizing an emitted result; the
3-symbol `UNIFORM_LIMIT` override; assembled-only F; a `S_P` ±-face-swap ablation.

---

## Diagnostic G — directionality (Eulerian, representation-QUALIFIED) + the adjointness independence audit

**Objects.** Construct the §3c increment `difference(closed_kernel, open_kernel)` at source and `extract` all SIX
directional weak blocks per `(α,ρ)` and every retained grade — named to the `extract` keys (`audit.py:325–342`):
`T→θ = TRANSVERSE_TO_THICKNESS/THETA`, `T→e_W = …/E_W`, `T→u_L = …/DIV_U`; `θ→T = THICKNESS_TO_TRANSVERSE/THETA`,
`e_W→T = …/E_W`, `u_L→T = …/DIV_U`. Direction = trial sector → tested row.

**G1 — witness table for ALL six blocks.** For **each** of the six blocks and each retained grade, `emit` a
witness-candidate row: weak-basis label, the weak-operator **coefficient** against arbitrary independent
compact-support trial/test functions, the cleared numerator, and any denominator/domain assumption. ⛔ Do NOT
extract by zeroing only undifferentiated fields (`S11c_b_SHARED_PHYSICS.md:325–327`). ⛔ `expr!=0` / expression size
/ an unevaluated `Integral` is NOT a witness. An **empty table is data** (a computed absence), ⛔ not a conclusion —
the script must not classify a block as "live"/"blocked" itself.

**G2 — per-face AND assembled.** Build the per-face increment `I_s` (the PINNED per-face set above) and `extract` the
six blocks from **`I_s`**, ⛔ NOT directly from `model['faces'][s]` (which is only `I_closed,s`). Report per-face and
assembled, so an assembled reverse-zero from face cancellation is not mislabeled structural.

**G3 — adjointness independence audit (settle the OPEN sub-question; ⛔ do not presume).** The pairing is c2 **§3b /
S11c-b §3c** — the stored-energy/kinetic **weak variational pairing** of S11c-b §1c (`S11c_b_SHARED_PHYSICS.md:312–339`),
⛔ NOT c2 §3d.4 (that is the traction/mechanical-power face-force-sign pairing — a different object). **The formal-adjoint
convention is PINNED (⛔ NOT the builder's choice):** the pairing is **bilinear** — the `extract` test×residual form
(`audit.py:325–342`) with **no conjugation**, so `i`, `ω`, and `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)` are carried **unchanged** (⛔ no
`ω→−ω`, ⛔ no complex conjugation). The formal adjoint of a weak block = **swap the `extract` trial/test slots and the
momentum legs `(X,k_out)↔(Y,k_in)`**, with the **face label `s` PRESERVED** (⛔ not flipped); the only IBP is the in-plane
compact-support one (boundary term 0, `S11c_b §3c:328–330`), and the `Z`/Fourier integrals are ⛔ NOT evaluated. `emit`
this convention as a machine-readable table alongside the adjoint so a leg can check it.
- Construct, by **genuinely independent source routes**, (a) the direct reverse block and (b) the formal adjoint of
  the independently-extracted forward block under that pairing.
- `emit` the pairing residual (`reverse − forward†`) and a **one-sided source corruption** that MUST satisfy — and the
  builder must `emit` evidence of each, since the concrete choice needs the CAS:
  **(i) one route's source COPY only** — mutate either the direct-reverse source OR the independently-extracted forward
  (before adjoint), ⛔ NEVER a shared pre-`extract` operand (that moves BOTH routes = a fake independence test);
  **(ii) a NON-INERT coefficient** — the builder selects, on face `s=+`, a named existing `δp_s` source coefficient that
  **demonstrably appears in the tested reverse block after `extract`** (⛔ e.g. a curl-free coefficient relocated inside
  the curl-projected `u_L→T` source is a no-op), and `emit`s which coefficient / FROM-row / TO-row it chose;
  **(iii) a genuine FORM change** — relocate that coefficient to a **dimension/retained-grade-matched** destination row
  AND flip its sign, ⛔ not a rescale. Recompute baseline + corrupted for BOTH routes and `emit` both literal differences
  (showing **only one route moved**) + a **machine-readable route-provenance/availability object**. ⛔ NEVER `C − B†` from two
  aliases of one precomputed block (an `A−A`); ⛔ "not implemented in the audit" is NOT proof no route exists; ⛔ no
  prescribed outcome target. The one-sided corruption **is** G's FORM ablation (⛔ do NOT add a trial/test
  sector-collapse — that degrades the probe, a special-field proxy, not the physical source).

**G4 — representation scope.** Construct G in the **Eulerian** representation only; `emit` results labeled by
object/grade/face/representation. ⛔ Do NOT claim representation-independent directionality — the material-coordinate
route is the not-yet-built N6 correction (`_measurements/S11c_c2_N6_spec_adjudication_record.md`); leave the
cross-representation comparison to it.

⛔ **Forbidden G proxies/overclaims:** directionality from `CLOSED_COUPLING_KERNEL`; reversed label orientation;
frozen/special trial fields, special momentum, uniform/zero-DtN limits, zeroing undifferentiated fields; an
assembled zero from face cancellation; a representation-dependent zero shown as invariant; a canonicalization zero
under unstated denominator/boundary assumptions; ⛔⛔ inferring reciprocity, passivity, **dissipativity**, stability,
or an energy sign (the construction is balance-laws + virtual work, not an ordinary action for the irreversible
response, `SHARED_PHYSICS:124–131`).

---

## Emit + acceptance discipline
- ⛔ **No acceptance criterion referencing an expected value is supplied.** Whether `I_closed|_uniform` weak-vanishes,
  which G blocks are live, and whether an independent adjoint route exists are **findings adjudicated on our side** —
  ⛔ never a builder exit condition, ⛔ never "iterate to zero." Every ablation/corruption `emit`s baseline + mutated
  operands + literal difference; ⛔ it does not assert "must change."
- ⛔ Use ONLY the named audit machinery; ⛔ do not reference any deleted `reduction/` helper. The FORM ablations
  (F7, G3) + the weak-coefficient/witness tables are the defect-catchers.
- Serialize heavy CAS (the increment machinery is memory-heavy); one job at a time; watch `free -h`.
