# Independent physics review — the S11c-c1 RECORD corrections (accuracy + completeness)

## Artifact (what to review)
Two committed c1 RECORDS, just edited with correction blocks (⛔ the ENGINES/EXPORTS are NOT edited and are NOT
under review — c1's computed objects STAND, adjudicated + independently Codex-verified as NO REOPEN). Review the
CORRECTIONS for accuracy, scope, and completeness:
- `research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md` — the c1 step record. New/changed: the *Method notes*
  rule-17 bullet (density O(εη)), a new *Exact grazing — claim SCOPED* bullet, a new *Independence is SCOPED* bullet,
  a scoping parenthetical on the "agreement is independent construction" line in *The two blind engines*, and a new
  *## Carry-forward corrections* subsection (7 lower-severity items).
- `research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md` — §3.5 (seal 5): the "harmless
  O(η²)" claim replaced with an O(εη)-recoverable-representational framing.

The exact edit set is the working-tree diff (`git -C /var/projects/toy_physics diff` against HEAD `f8509b7a`).

## What to check (this is a DOCUMENT / physics-bearing-prose review)
These are RECORD corrections earned by a retrospective review of the c1 SHARED_PHYSICS spec. Your job: confirm each
correction ACCURATELY states the physics and correctly scopes what may be CLAIMED, does not over- or under-claim,
and that the set is COMPLETE (no claim-level MUST-FIX from the two source reports was dropped or misstated).

Settle these load-bearing points with your OWN computation where a computation is cheap (SymPy), and quote the
source both ways:

1. **Density channel order.** Is `d(μ_s)/dη|₀ = −μ_θ·w₁/ρ_br` with `μ_s=μ_θ/ρ`, `ρ=ρ_br(1+η·w₁)`, `μ_θ=O(ε)` — i.e.
   a channel that is **FIRST-order in the shape η** and overall **O(εη)** — and is the earlier "harmless O(η²)"
   reading therefore WRONG? (Both source legs claim so; verify independently.)
2. **Recoverability / no-reopen.** The corrections assert PY carries `rho_br_bg_rho4_constant` **opaquely** (a bare
   `1/ρ` symbol, **0 derivatives**, never combined with `w1_profile`), so re-binding it to the live
   `background_density_map` relation recovers the EXACT O(εη) channel downstream in c2 ⇒ c1 need NOT reopen. Is that
   correct? Inspect the REAL export: fold it via
   `ledger_fold.load_model("scripts/S11c_b_exports.py","scripts/S11c_c1_exports.py")` (module
   `research/pde_ledger_v3/reduction/ledger_fold.py`), read `s11c_c1_face_response_coeffs` for the RHO4_CONSTANT
   MU_THETA_COEFFICIENT and confirm the density is an opaque reciprocal with no `w1_profile` and no derivative, and
   that `background_density_map` carries `rho_br_bg_rho4_constant = W_bg·ρ_br/W_0 = ρ_br(1+η·w₁)`; then substitute
   and expand to O(η) and confirm the recovered coefficient is `−μ_θ·w₁/ρ_br`. ⭐ Is anything about PY's exported
   RHO4_CONSTANT response **irrecoverable** by a downstream re-bind? If YES → the no-reopen framing is wrong, say so
   plainly. If NO → the record's "STANDS, no reopen, recoverable in c2" framing is right. (This override was already
   independently Codex-sol-verified with residual 0; ⛔ do not merely trust that — but do not spin on it either: one
   clean recovery computation settles it.)
3. **Grazing scoping.** Is the face-response inverse `[I+(Λ_A/ρ_m²)Z]⁻¹` genuinely **nonanalytic as `q_out→0`**
   (~`1/η`, not a Taylor series in η), so that marking exact grazing `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER` — `Z₁` a
   valid non-grazing asymptotic coefficient only, `‖N₀⁻¹N₁‖≪1` on both legs — is the correct and sufficient claim
   scope (⛔ not an over- or under-statement)?
4. **Independence scoping.** The c1 spec supplied the composition recipe + some expected values (rigid-shift
   cancellation, flat `Z₀=ρ_m ω/q_out`, zero-jet) ⇒ for those objects the cross-engine agreement is partly fidelity
   to supplied structure, ⛔ not independent discovery — while the two-momentum DtN kernel IS independently confirmed
   (both blind engines + both retro legs re-derived it). Is the record's scoping of the "agreement is independent
   construction" claim accurate — neither erasing the genuine kernel confirmation nor letting the supplied-value
   agreement pose as discovery?
5. **Carry-forward completeness + accuracy.** The 7 lower-severity items (energy-residual orientation `P_face+P_∞=0`
   / minus outgoing Poynting; `h_s` graph-vs-outward + `N`/`Z` terminology; density-as-multiplication-operator
   naming for c2; `K_a=(Z−Z†ₐ)/(2i)` is **Hermitian** not anti-Hermitian; evanescent caveat = η²/η·σ_W/σ_W²;
   drain-projection `O(σ_W²)` applies to drain-TILT not convection; flat-`Z₀`/rigid-shift rule-5 leakage). Is each
   dispositioned ACCURATELY against the two source reports, and is this the COMPLETE set — does either source report
   carry a claim-level (SHOULD/NIT/MUST) item the corrections dropped or misattributed?

## What you are handed (read these SOURCES first, form your own view, THEN read the corrected records)
- The two source leg reports: `research/pde_ledger_v3/directives/_measurements/S11c_c1_spec_retro_review_grok.txt`
  and `…_codex_sol.txt` (one CLEAR, one BLOCK-reopen; the corrections adopt the record-only resolution).
- The adjudication that resolved them: `…/_measurements/S11c_c1_spec_retro_review_adjudication.md`.
- The real engine + fold: `scripts/S11c_c1_exports.py`, `scripts/S11c_b_exports.py`, `reduction/ledger_fold.py`.
- The c1 SHARED_PHYSICS spec (context for the source line-refs): `directives/S11c_c1_SHARED_PHYSICS.md`.
⛔ There is no do-not-read list; anything you must not use, you have not been given. In particular ⛔ this review is
about RECORD ACCURACY, not re-deriving the whole closure — but where a cheap SymPy check settles points 1–3, run it.

## Required method (DOCUMENT branch)
Read the source of truth FIRST (the two source reports + the adjudication +, for points 1–2, the real export), form
your own view of what the corrections SHOULD say, and only THEN read the corrected records. Quote both sides for
every finding. ⛔ A prose "I re-derived it and agree" is worth nothing — for points 1–3 save your SymPy script AND
its literal stdout to named absolute paths under `/tmp` and cite them; without them a derivation claim is discarded.
⛔ Do not modify the working tree (copy anything you run to `/tmp`).

## Physics filter
Report a finding ONLY if it catches a way the corrected RECORD (a) misstates the physics, (b) over- or under-scopes
what may be CLAIMED, or (c) drops/misattributes a claim-level item one of the two source reports actually raised.
⛔ Do not report "the record would be wrong for a different problem," and ⛔ do not re-litigate settled decisions
unless you have a concrete computation showing the override (no-reopen) is wrong.

## Output
Ranked findings (MUST-FIX / SHOULD-FIX / NIT), each with the quoted record text, the quoted source, and (for 1–3)
your script+stdout paths. End with an explicit verdict: **CLEAR** (the corrections are accurate, correctly scoped,
and complete) or the exact list to fix. If you conclude the no-reopen override is actually wrong, say **REOPEN c1**
and show the irrecoverable channel.
