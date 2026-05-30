---
stage: 139
unit_id: 139
batch: IV.5
created_at: 2026-05-27T00:00:00Z
rewritten_at: 2026-05-29T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-29T22:53:03Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 139 (REWRITTEN to encode codex_review findings)

Authoritative source: `redteam/codex_reviews/stage_139.md` (verdict FINDINGS, 4 findings R1–R4).
This directive SUPERSEDES the prior 2026-05-27 version, whose F1 "outlet consistency"
and F2 "R_q^comp = 1/4" checks the review flagged as still-tautological. The prior
banner item (old F4) is folded into F-BANNER below. The F2 math direction was settled by a
Claude+Codex read-only consult (all CONCUR; record `redteam/codex_reviews/_consult_batch5.md`,
Q4/Q5).

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append
`## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new numeric literals. Every retained literal traces to its owning upstream
stage (`r_F1` → Stage 121; `Pi_*`, `S_q(Pi_*)` → Stage 134; `g_-^F1` is the canonical
cross-stage constant owned at 127/142/144/164/169). Do NOT generate any target by evaluating
this script's own definition of the quantity it checks.

Do NOT introduce unrelated features, refactors, or stylistic changes. Edit only the regions named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>`
for Mathematica) and iterate until they exit 0 with all in-file checks passing. Use
`timeout 600`. Getting the scripts to run cleanly is your job; the orchestrator independently
re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — tautological_check (R1: outlet-consistency residual is X−X)

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:54-55`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:77-78`

**Issue:**
The two outlet-consistency assertions check `Pi_* - (M_s + M_q * S_q)` while `M_s` is DEFINED
as `Pi_*/(1 - R_q*S_q)` and `M_q = -R_q*M_s` in the SAME script (sympy lines 12-13, 17-18;
mathematica lines 35-36, 51-52). Substituting those definitions:
`M_s + M_q*S_q = M_s(1 - R_q*S_q) = [Pi_*/(1 - R_q*S_q)](1 - R_q*S_q) = Pi_*`, so the residual
collapses to `Pi_* - Pi_*` exactly, for BOTH the natural and compensated branches. The check
cannot fail and does not exercise the imported `Pi_*` / `S_q(Pi_*)` literals. (Claude+Codex
consult Q4: CONCUR — these are identities by the gain definitions.)

**Required change:**
Replace the two structurally-zero outlet residuals with a GENUINELY INDEPENDENT check that
re-derives `S_q(Pi_*)` from the closed-form mouth-response kernel that Stage 134 owns, and
confirms it matches the imported literal `Sq_star`. Stage 134 fixes
`S_q(Pi) = S(Pi, pi/2)` with
`S(Pi, kappa) = Pi*(kappa*tanh(kappa) + Pi*(exp(-Pi)/cosh(kappa) - 1)) / ((1 - exp(-Pi))*(kappa^2 - Pi^2))`
(see `notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md` and
`scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:22-30`). Evaluating
that kernel at `Pi = Pi_*` with `kappa = pi/2` is a path that does NOT reuse the gain
definitions, so the comparison to `Sq_star` can fail if either the imported literal or the
kernel transcription is wrong.

In the SymPy script, REPLACE lines 51-55 (the comment block and the two
`assert abs(Pi_star - (Ms_* + Mq_* * Sq_star)) < tol_algebraic` lines) with:

```python
# R1 (was tautological outlet residual): independently reconstruct S_q(Pi_*) from
# the Stage 134 closed-form mouth-response kernel S(Pi, kappa) at kappa = pi/2, and
# confirm it equals the imported Stage 134 literal Sq_star. This exercises the
# imported value via a route that does NOT reuse M_s = Pi_*/(1 - R_q S_q).
kappa_q = sp.pi/2
S_kernel = (Pi_star * (kappa_q*sp.tanh(kappa_q)
            + Pi_star*(sp.exp(-Pi_star)/sp.cosh(kappa_q) - 1))
            / ((1 - sp.exp(-Pi_star)) * (kappa_q**2 - Pi_star**2)))
Sq_recon = sp.N(S_kernel, 30)
assert abs(Sq_recon - Sq_star) < tol_literal, (Sq_recon, Sq_star)

# Documented identity (NOT independent): with M_s = Pi_*/(1 - R_q S_q) and M_q = -R_q M_s
# the outlet form Pi_* = M_s + M_q S_q is true by construction; kept only as a structural
# sanity print, NOT as a verification of the imported literals (see R1 in codex_review).
print('outlet form residual (nat, structural) =',
      sp.N(Pi_star - (Ms_nat + Mq_nat * Sq_star), 5))
print('outlet form residual (comp, structural) =',
      sp.N(Pi_star - (Ms_comp + Mq_comp * Sq_star), 5))
```

In the Mathematica script, REPLACE lines 77-78 (the two
`expectApprox["outlet consistency ...", piStar, m*Star + m*Star*sQStar, tolAlg];` lines) with:

```mathematica
(* R1 (was tautological outlet residual): independently reconstruct S_q at piStar
   from the Stage 134 closed-form kernel S(p, kappa) at kappa = Pi/2, and confirm it
   matches the imported Stage 134 literal sQStar via a route that does NOT reuse
   mS = piStar/(1 - rQ sQStar). ASCII-safe names only in comments/strings. *)
kappaQ = Pi/2;
sQRecon = N[piStar*(kappaQ*Tanh[kappaQ]
            + piStar*(Exp[-piStar]/Cosh[kappaQ] - 1))
            / ((1 - Exp[-piStar])*(kappaQ^2 - piStar^2)), 30];
expectApprox["S_q recon from Stage 134 kernel", sQRecon, sQStar, tolLit];
(* Structural sanity only (true by construction, NOT a literal check): *)
Print["outlet form residual nat structural = ", fmt[N[piStar - (mSNat + mQNat*sQStar), 5]]];
Print["outlet form residual comp structural = ", fmt[N[piStar - (mSComp + mQComp*sQStar), 5]]];
```

**Anti-tautology guard:**
The new `Sq_recon` / `sQRecon` is computed from the closed-form kernel only (tanh, sech of
`pi/2`, and `Pi_*`); it never references `Ms_nat`, `Ms_comp`, `Mq_*`, or `R_q`. Therefore the
assertion `|Sq_recon - Sq_star| < tol_literal` can fail if `Sq_star` is mistyped or the kernel
is wrong — it is NOT X−X. The old `Pi_* - (M_s + M_q S_q)` residuals are demoted to non-assert
`print` lines explicitly labelled structural, so they no longer masquerade as verification.
Do NOT re-introduce them as asserts. CONFIRM the kernel expression character-for-character
against `scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:22-30`
(SymPy) and `mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:32`
(Mathematica) — note Stage 134 writes the denominator as `(kappa^2 - p^2)`, equal to
`(pi^2/4 - Pi^2)` at `kappa = pi/2`; use the `(kappa^2 - p^2)` form to mirror the owner exactly.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 139` and `redteam exec-mathematica 139`
(one at a time per single-seat policy) and confirms the sympy transcript contains a passing
`S_q recon` assertion (no AssertionError) and the mathematica transcript contains
`PASS: S_q recon from Stage 134 kernel`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`
- summary: Replaced the tautological outlet-consistency assertions with independent Stage 134 kernel reconstruction checks and demoted the outlet residuals to structural prints.
- deviation: none

---

## F2 — tautological_check (R2: R_q^comp = 1/4 is DEFINITIONAL and branch-blind)

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:15-16,58`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:38-49,79`

**Issue:**
The compensated `R_q^comp = 1/4` claim is tautological in BOTH engines, AND it is branch-blind:
- SymPy (line 15) sets `g_minus = rF - sqrt(1 + rF**2)/2`, so
  `R_q^comp = (g_minus - rF)^2/(1 + rF^2) = (sqrt(1+rF^2)/2)^2/(1+rF^2) = 1/4` ALGEBRAICALLY for
  any `rF`; the line-58 assert cannot fail.
- Mathematica (lines 38-49) "derives" `g_minus` by `Solve[(gc - rF)^2 == (1 + rF^2)/4 && gc < rF]`.
  But `R_q^comp = (gc - rF)^2/(1 + rF^2)`, so the solving condition `(gc - rF)^2 == (1 + rF^2)/4`
  IS `R_q^comp == 1/4` rearranged; the line-79 check is X−X with the answer pre-baked in.
- Critically, `R_q = 1/4` holds for BOTH roots `g_± = rF ± sqrt(1+rF^2)/2` (the `±` is squared
  away), so `R_q^comp = 1/4` does NOT even discriminate the lower (compensated) branch from the
  upper one. A sign typo `g_+ = rF + sqrt(1+rF^2)/2` would still pass `R_q = 1/4`.

**RESOLUTION (Claude+Codex consult, all CONCUR — `redteam/codex_reviews/_consult_batch5.md` Q5):**
`R_q^comp = 1/4` is the DEFINITION of the compensated branch (notes stage139 lines 91-100: the
compensated branch is the one on which `g_c = r - (1/2) sqrt(1 + r^2)`, whence `R_q = 1/4`). It
CANNOT "emerge" non-tautologically — there is no independent route by which `1/4` is a derived
surprise; it is definitionally true for every `r`, and for both branches. Do NOT manufacture a
fake "emergence" (deriving `g_c` from the `1/2`-coefficient family still gives `1/4` for any `r`
because `(1/2)^2 = 1/4` analytically). Instead:
1. RELABEL `R_q^comp = 1/4` in both engines as a DEFINITIONAL-CONSISTENCY check (kept, explicitly
   NOT the anchor).
2. The LOAD-BEARING falsifiable checks are the ones that CAN fail for a wrong `rF` or branch:
   - `r_F1` itself vs the Stage 121 closed form (already asserted, sympy py:38 / wl:71) — the
     root input the whole branch hangs on;
   - the NATURAL-branch `R_q^nat = 0.145454452260421 ≠ 1/4` (already asserted, py:41 / wl:72) —
     a genuinely DIFFERENT computation `(1 - rF)^2/(1 + rF^2)`, NOT forced to any value; this is
     what proves the branch SELECTION is real and non-vacuous;
   - the `g_-^F1` VALUE vs the canonical cross-stage literal `0.758035078944662826919680890414`
     (owned at 127/142/144/164/169) — this DISCRIMINATES the lower branch (`g_- ≈ 0.758`) from
     the upper (`g_+ ≈ 2.79`), which `R_q = 1/4` cannot; a sign/branch typo gives `≈2.79` and FAILS;
   - the compensated mouth gains `M_s^comp,*`, `M_q^comp,*` (already asserted, py:48-49 / wl:75-76)
     which depend on `R_q^comp` and the imported `Pi_*`, `S_q(Pi_*)`.
3. For the Mathematica transliteration: COMPUTE `g_minus` DIRECTLY as the closed form
   `rF - Sqrt[1 + rF^2]/2` (a sanctioned mirror of the SymPy closed form — accepted under the
   mirror policy because it is an algebraic expression in the independently-anchored `rF`, NOT a
   transcendental root, AND it is cross-checked against the independent numeric branch literal),
   NOT by `Solve`-ing the `(gc - rF)^2 == (1 + rF^2)/4` condition (which re-bakes `1/4`).

**Required change (SymPy):**

REPLACE line 15:
```python
g_minus = sp.N(rF - sp.sqrt(1 + rF**2)/2, 30)
```
with (computation unchanged; provenance/definitional comment added):
```python
# g_- is the LOWER compensated branch g_c = rF - (1/2) sqrt(1 + rF^2) (notes stage139 lines
# 91-100). R_q^comp = 1/4 is DEFINITIONAL on this branch (true for any rF) AND branch-blind
# (g_+ gives 1/4 too). The falsifiable content is the g_-^F1 VALUE (discriminates the lower
# branch, below) and the natural-branch R_q^nat != 1/4 (line 41), NOT R_q^comp = 1/4.
g_minus = sp.N(rF - sp.sqrt(1 + rF**2)/2, 30)
```

After the `R_q^nat` literal assert (line 41), ADD the branch-discrimination anchor:
```python
# (R2 anchor) g_-^F1 VALUE vs the canonical cross-stage literal (owned at 127/142/144/164/169).
# This DISCRIMINATES the lower compensated branch (g_- ~ 0.758) from the upper (g_+ ~ 2.79),
# which R_q = 1/4 cannot. Falsifiable: a sign/branch/rF typo gives ~2.79 and FAILS.
assert abs(g_minus - sp.Float('0.758035078944662826919680890414', 30)) < tol_literal, (g_minus,)
```

REPLACE the line-58 block:
```python
# Compensated R_q closed form
assert abs(Rq_comp - sp.Rational(1, 4)) < tol_algebraic
```
with (relabeled as definitional-consistency, kept):
```python
# (R2 definitional-consistency, NOT the anchor) R_q^comp = 1/4 holds by construction of the
# compensated branch (g_- = rF - sqrt(1+rF^2)/2 ⇒ (g_- - rF)^2/(1+rF^2) = 1/4 for any rF, both
# branches). The falsifiable checks are the g_-^F1 value (line ~42) and R_q^nat != 1/4 (line 41).
assert abs(Rq_comp - sp.Rational(1, 4)) < tol_algebraic, (Rq_comp,)
```

**Required change (Mathematica):**

REPLACE the Solve block at lines 38-49 (the `(* Derive g_c ... *)` comment, the
`gMinusSolutions = gc /. Solve[(gc - rF)^2 == (1 + rF^2)/4 && gc < rF, gc, Reals]` line, the
`If[Length[gMinusSolutions] =!= 1, ...]` guard, `gMinus = N[First[gMinusSolutions], 30]`, and
the X−X `expectApprox["g_minus closed form", gMinus, rF - Sqrt[1 + rF^2]/2, 10^-25]`) with the
direct closed form ONLY:
```mathematica
(* g_- is the LOWER compensated branch g_c = rF - (1/2) Sqrt[1 + rF^2] (notes stage139 lines
   91-100). R_q^comp = 1/4 is DEFINITIONAL on this branch (true for any rF) AND branch-blind.
   Compute g_- DIRECTLY as the closed form (a sanctioned mirror of the SymPy route; rF is
   independently anchored at line 71 and the branch is discriminated by the g_-^F1 value check
   below) — NOT by solving (gc - rF)^2 == (1 + rF^2)/4, which would re-bake 1/4. *)
gMinus = N[rF - Sqrt[1 + rF^2]/2, 30];
```

Then ADD the branch-discrimination anchor INTO THE ASSERTION GROUP (after the
`expectApprox["R_q^nat", rQNat, 0.145454452260421, tolLit];` line at line 72, where `tolLit` is
already defined at line 68):
```mathematica
(* (R2 anchor) g_-^F1 VALUE vs the canonical cross-stage literal (owned at 127/142/144/164/169);
   DISCRIMINATES the lower branch (g_- ~ 0.758) from the upper (g_+ ~ 2.79), which R_q = 1/4
   cannot. Falsifiable: a sign/branch/rF typo gives ~2.79 and FAILS. *)
expectApprox["g_-^F1 value", gMinus, 0.758035078944662826919680890414, tolLit];
```

Line 79 — relabel the definitional-consistency check (keep the assert, change ONLY the label):
```mathematica
expectApprox["R_q^comp - 1/4 (definitional-consistency)", rQComp, 1/4, tolAlg];
```

**Anti-tautology guard:**
`R_q^comp = 1/4` is acknowledged as definitional/branch-blind and relabeled in BOTH engines.
The genuine falsifiable content is: (a) `r_F1` vs the Stage 121 closed form (independent); (b)
`R_q^nat ≠ 1/4`, a different computation that proves the selection is non-vacuous; (c) the
`g_-^F1` value, which discriminates the lower branch from the upper (a `g_+` sign typo gives
`≈2.79` and FAILS). The Mathematica `gMinus` MUST be the direct closed form
`rF - Sqrt[1 + rF^2]/2` — do NOT re-introduce `Solve[(gc - rF)^2 == (1 + rF^2)/4 ...]` (that is
the transliteration R2 flagged) and do NOT write `1/4` or `(1 + rF^2)/4` into the `g_minus`
derivation. The `g_-^F1` literal is the cross-stage canonical value (owned at 127/142/144/164/169),
NOT this script's own frozen output — it is the branch-discriminating anchor, not an independent
re-derivation of `g_-` (g_- has no independent definition; it IS `rF - sqrt(1+rF^2)/2`).

**Verification command:**
Verifier confirms the sympy transcript shows the new `g_-^F1 value` assert and the relabeled
`R_q^comp` assert both passing, and the mathematica transcript shows
`PASS: g_-^F1 value` and `PASS: R_q^comp - 1/4 (definitional-consistency)`, AND inspects that the
Mathematica `gMinus` is computed as the direct closed form `rF - Sqrt[1 + rF^2]/2` (NO `Solve`,
no `(1 + rF^2)/4` condition anywhere).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`
- summary: Added the branch-discriminating `g_-^F1` literal anchor, relabeled `R_q^comp = 1/4` as definitional consistency, and replaced the Mathematica branch solve with the direct lower-branch closed form.
- deviation: none

---

## RESOLVED: R3 (paper_misalignment) — NO SCRIPT FIX

The codex_review R3 says the paper card requires checking the "self-matched susceptibility
closure before using the one-scalar branch law," but the scripts never evaluate a
susceptibility-closure residual.

RESOLUTION: This is forward-reference boilerplate on the Stage 139 card. The self-matched
susceptibility closure is the deliverable of **Stage 140** (`paper/stages/stage_140.tex`:
"Self-Matched Mouth Susceptibility Closure", boxed result `Sigma_0 = M_s = (20/9) That_m^2`),
which carries `\StatusReduced / \StatusExactClosure`. Stage 139 only EVALUATES the gain pair on
the Family-1 branch; it does not establish or use that closure. The shared `Checks` checklist on
the 139 and 140 cards is identical boilerplate. Logged as paper-card cleanup item **P4-42
(139 -> 140)**: the closure Checks line belongs on the 140 card, not 139.

DO NOT add any susceptibility-closure residual check to the Stage 139 scripts. Only add at most a
one-line provenance COMMENT noting where the closure is verified. Concretely:

In the SymPy script, after the `Sq_star = sp.N('0.658075937605429', 30)` line, ADD the comment:
```python
# (Self-matched susceptibility closure Sigma_0 = M_s is established at Stage 140, not here;
#  Stage 139 only evaluates the gain pair on the Family-1 branch. See P4-42.)
```
In the Mathematica script, after the `sQStar = SetPrecision[...]` line, ADD the comment:
```mathematica
(* Self-matched susceptibility closure (Sigma_0 = mS) is established at Stage 140, not here;
   Stage 139 only evaluates the gain pair on the Family-1 branch. See P4-42. *)
```
ASCII-safe names only; do NOT use `_*)` or `Pi_*)` or a bare `*)` substring inside the comment
body — the Mathematica comment would close early (this exact pitfall hit Stage 139 in a prior
batch). The text above is already safe.

## Applied: R3

- files_changed:
  - `scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`
- summary: Added the requested provenance comments clarifying that self-matched susceptibility closure is established at Stage 140, not in the Stage 139 gain evaluation.
- deviation: none

---

## RESOLVED: R4 (other / provenance stage numbers) — SCRIPTS ALREADY CORRECT

The codex_review R4 says the prior directive demanded provenance comments citing "Stage 223"
(for `r_F1`) and "Stage 236" (for `Pi_*`, `S_q(Pi_*)`), while the edited scripts cite "Stage 121"
and "Stage 134" with no note explaining the deviation.

RESOLUTION: The scripts' **Stage 121 / Stage 134 are the CORRECT post-renumber canonical
numbers.** A global renumber applied a -102 offset to body citations in the 220-251 range:
`223 - 102 = 121` and `236 - 102 = 134`. The prior directive's 223/236 were the STALE
pre-renumber numbers. Verified against the canonical notes:
- `r_F1`: Stage 121 owns it.
  `notes/stages/moving_throat_pde_stage121_geometric_r_selection.md:64-72` boxes
  `r_F1 = sqrt(12/pi^2 (37/20)^2 - 1) = sqrt(4107 - 100 pi^2)/(10 pi) ~= 1.77799353547498`.
  (The notes' own body still says "before Stage 223" at line 97 — that is the pre-renumber
  inline citation, confirming the offset.)
- `Pi_*` and `S_q(Pi_*)`: Stage 134 owns them.
  `notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md:71-88` give
  `Pi_* ~= 1.50882951349316` (from carried "Stage 233" pre-renumber) and
  `S_q(Pi_*) ~= 0.658075937605428`, with the closed-form kernel at lines 49-61. Stage 134's own
  codex_review corroborates its `S_q` checks.

The scripts ALREADY carry the right provenance comments (sympy lines 5,7; mathematica lines
28,30 say Stage 121 / Stage 134). NO CHANGE to those comments. This block documents that the
121/134 numbers are authoritative and that the prior directive's 223/236 demand is void. Codex
must NOT "correct" the scripts back to 223/236.

## Applied: R4

- files_changed: []
- summary: Verified the scripts retain the canonical Stage 121 and Stage 134 provenance comments and did not reintroduce stale Stage 223 or Stage 236 labels.
- deviation: none

---

## F-BANNER — stale_output (banner-typo guard; folds prior F4)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:26`

**Issue:**
A prior batch found the Mathematica banner mislabelled as a different stage number. The current
file at line 26 reads `banner["STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS"];` which is ALREADY the
canonical Stage 139 label.

**Required change:**
VERIFY line 26 (and any banner/title string in either script) reads exactly
`STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS`. If it already does (expected), make NO change and
record `deviation: none (banner already canonical)`. If a stale `STAGE 1xx` (e.g. 122) banner
is present in either engine, relabel it to `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS`. Do NOT
put any `*)` / `_*)` / `Pi_*)` substring inside a Mathematica comment while editing nearby lines.

**Verification command:**
Verifier confirms the mathematica transcript banner line reads
`STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS` and the SymPy script carries no conflicting stage label.

## Applied: F-BANNER

- files_changed: []
- summary: Verified the Mathematica banner already reads `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS` and found no conflicting SymPy stage label.
- deviation: none (banner already canonical)

---

## Combined verification (all findings)

After Codex applies and iterates to clean exit, the verifier runs `redteam exec-sympy 139` then
`redteam exec-mathematica 139` (never >2 concurrent `math -script`; one at a time here) and
confirms:
- both scripts exit 0;
- the sympy transcript ends with `all assertions passed` and shows the new `S_q recon` assert
  passing (no AssertionError), the new `g_-^F1 value` assert passing, and the relabeled
  `R_q^comp` definitional-consistency assert passing;
- the mathematica transcript banner reads `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS`;
- the mathematica transcript contains, at minimum, the PASS lines:
  `PASS: S_q recon from Stage 134 kernel`, `PASS: g_-^F1 value`,
  `PASS: R_q^comp - 1/4 (definitional-consistency)`, plus the boxed-literal PASS lines (`r_F1`,
  `R_q^nat`, `M_s^nat,*`, `M_q^nat,*`, `M_s^comp,*`, `M_q^comp,*`);
- the Mathematica `gMinus` is computed as the direct closed form `rF - Sqrt[1 + rF^2]/2` —
  there is NO `Solve` for `gMinus` and NO `(1 + rF^2)/4` condition anywhere;
- no `Stage 223` / `Stage 236` provenance string was reintroduced (R4 — must stay 121 / 134).
