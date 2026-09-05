# Independent physics review — S11c-c1 SHARED_PHYSICS (RETROSPECTIVE spec gate, review-until-clear)

You are one of two independent legs reviewing an **orchestrator-written physics spec**. This is a **retrospective**
gate: the S11c-c1 spec was folded once from its decision review and **never re-reviewed to clear**, while everything
downstream of it (both engines, the comparator, the reconcile, the step record) was built and reviewed. The engine
legs checked **fidelity to this spec** and **cross-engine agreement** — so **a physics error in the spec itself would
be invisible to them: both engines would faithfully execute it and agree on the same wrong thing** (CLAUDE.md rule 7).
Your job is to find any such error, or confirm there is none. ⛔ Do not rubber-stamp; ⛔ do not soften to agree.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md`

## What S11c-c1 is
The curved-interface bulk closure. It solves the perturbed curved two-face outgoing bulk acoustic problem for: (1) the
nonlocal **DtN/impedance** operator `Z_s(ω;k,k′)` — a **two-momentum** operator carrying both branch legs
`q_out(k),q_out(k′)`; (2) the **permeable curved face response** `(δp_s,J_s,t_s)(V_s,μ_θ)` as an operator inverse
`[I+(Λ_A/ρ_m²)Z]⁻¹`; (3) the Fredholm noninvertibility loci; (4) a three-object dissipation audit. It generalizes
S11b's flat-face B0b/B0c to the tilted faces S11c-a shape-differentiated.

## What you are handed (derive the physics YOURSELF first, then read the spec, then cross-check the outcome)
Read in this order (a method instruction, not a blindness control):
1. Parent decision list `directives/S11c_decisions.md` (the ratified N-series scope, esp. N3/N4/N11).
2. The inherited physics: S11b spec `directives/S11b_SHARED_PHYSICS.md` + step `steps/S11b_interface_coupling_law.md`
   (flat B0b/B0c: impedance, three regimes, permeable response); S11c-a spec + step (the tilted-face shape
   derivatives T-a..T-i); S11c-b spec + step (the slab operator, the μ_θ operand).
3. **Derive independently, from those sources (NOT from the c1 spec):** the curved-face DtN as a two-momentum
   operator, the permeable face response by the operator inverse, and the dissipation objects. Save your derivation
   script(s) + literal stdout (see method).
4. **Then** read the c1 spec under review §§1–6.
5. **Then** cross-check against the achieved outcome: the c1 step record `steps/S11c_c1_curved_bulk_closure.md` and the
   reconcile `directives/_measurements/S11c_c1_comparator_reconcile.md` (what cross-engine AGREED vs was left
   UNDECIDED/DEFERRED). Use these to test whether the spec's decisions were sound *given what the build found* — e.g.
   was leaving seal-5 (background density) UNDECIDED the right call, or did the spec force a rule-17 freeze?

## Required method (DOCUMENT branch + computation)
Form your own independent view from the inherited sources first, then read the spec. Quote both sides for every
finding. ⛔ **A prose derivation is worth nothing.** Where a claim is computational, **write your own SymPy script and
save both the script and its literal stdout to named absolute paths** (`/tmp` or the `_legs/` dir); without them your
derivation claims are discarded. Do NOT spawn Mathematica kernels. Copy anything you run to `/tmp`; never modify the
working tree.

## What to hunt (the spec's DECISIONS — is the physics right?)
Check the spec's **SUPPLIED setup (§§1–2)** for a physics error both engines would inherit, and its **objects-to-
compute + controls (§§3–5)** for defects. In particular settle, with computation where possible:

1. **The two-momentum DtN operator (§3a).** Is `Z_1 ∝ Ŵ_bg(k−k′)·[q(k)q(k′)+k·k′−ω²/c_s0²]/(q(k)q(k′))` the correct
   first-shape-order curved-face DtN? Is the composition `N₀∘M_{h}∘N₀ + Div(h∇) + κ²h` right? Does the radiation
   branch / outward-normal orientation (`n̂_s = (−½∇W_bg, s)`, ⛔ not `sign(∇₄F)`) give the correct sign, and is the
   rigid-shift (`k=k′`, `Ŵ(0)`) cancellation genuine?
2. **The permeable face response (§3b).** Is the operator inverse `[I+(Λ_A/ρ_m²)Z]⁻¹` the correct elimination of the
   coupled `{δp_s = Z·v_bulk,s ; v_bulk,s = V_s + J_s/ρ_m ; J_s = Λ_A𝒜_s+Λ_V V_s ; t_s = −(δp_s+Λ_X𝒜_s)n̂_s}` system?
   Is the Λ-channel placement (`Λ_X` in traction only; `Λ_A,Λ_V` in flux) physically correct?
3. **The dissipation audit (§3b).** Are the three objects genuinely distinct and correctly defined — the bulk-radiation
   Hermitian part `H_a[Z]`; the two-port permeable Hermitian form; the INDEPENDENT traction-vs-far-field-Poynting
   energy balance (face operand the true-area traction pairing, ⛔ not `½Re(δp V*)`)? Is the O(η²) evanescent-nullspace
   caveat correct, or does it mask a real first-order error?
4. **The rest-frame limit + grazing (§2b).** Is dropping the convective operator (`N11a`) legitimate at first shape
   order? Is the non-uniform grazing validity domain correct?
5. **Rule-17 / freezes.** Does the spec anywhere freeze a varying field (background density, `W_bg`, a jet, a branch
   leg) that should stay live? Was leaving seal-5 UNDECIDED correct, or did the spec's `ρ` handling force a freeze?
6. **Rule 5 / rule 3 / rule 2.** Any leaked expected value; any manufactured recipe (specifying a derivation path
   instead of naming the object); any tautological residual or A−A control (the build found several emit-only A−A
   controls — were any rooted in a SPEC defect, not just an engine defect?).
7. **Inheritance correctness (rule 16).** Does the spec correctly inherit S11b's B0b/B0c and S11c-a's T-a..T-i, or
   does it silently re-point a standard name or import a superseded relation?

## Physics filter
Report a finding only if it catches a way the spec's **decisions** make the physics wrong, under-specified, or
un-reviewable. Rank most-load-bearing first; mark each MUST-FIX / SHOULD-FIX / NIT with the exact spec line and the
source line it contradicts. For any MUST-FIX, state the **downstream impact**: does it change what was COMPUTED (c1
must re-open) or only what may be CLAIMED (step-record caveat)?

## Output
Settle items 1–7 with your verdict + the computation/quote. Then the ranked findings. End with an explicit verdict:
**CLEAR (the c1 spec's physics decisions are sound; no MUST-FIX)** or **BLOCK** with the must-fix list and each one's
downstream impact. If it is genuinely clean, say so plainly — a real confirmation from an independent first-principles
derivation is the goal.
