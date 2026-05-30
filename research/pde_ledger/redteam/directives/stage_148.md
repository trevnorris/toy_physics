---
unit_id: 148
batch: IV.5
created_at: 2026-05-29T00:00:00Z
supersedes: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-30T04:53:38Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 148 (rewrite encoding Codex review of 2026-05-29 + consult batch 7)

## What this rewrite does

The 2026-05-27 directive (now superseded) is STALE: its prose, quoted code,
claim-manifest, and self-test all carry the wrong radical constant `168π²` and an
under-strength tolerance, and its F3 prescribes a Mathematica `dSigmaOfDeltas`
restructure that — as actually applied to the LIVE script — silently DROPS the
S-follows-Π chain term, so the two engines now compute DIFFERENT first-order
T-corrections while both still pass. A Codex read-only review
(`redteam/codex_reviews/stage_148.md`, findings_count: 2) plus a Claude+Codex consult
(`redteam/codex_reviews/_consult_batch7.md`, Q1/Q2 — both DECIDED) fix the math. This
directive re-anchors every file:line to the LIVE scripts and encodes exactly two
findings.

The stage's mathematical content lives in the paper card and its appendix. The operative
external anchors used below are:

- First-order traction-shift model
  (`paper/appendices/stage_appendix_part04.tex:838-841`,
  `eq:app-part04-deltaT-firstorder`):
  `δT̂_m = ε[ A_T(ḡ_ς - 𝔤_*) + B_T(S̄_ς - 𝒮_*) ]`.
- Coefficient literals
  (`paper/appendices/stage_appendix_part04.tex:846,848`,
  `eq:app-part04-AT-BT-values`):
  `A_T ≈ -4.27263956256927`, `B_T ≈ 0.134875005736706`.

These literals are an EXTERNAL anchor (verbatim from the published appendix), NOT a
SymPy mirror — they are the independence witness that pins the Mathematica derivation.

Apply each finding below in order. After applying, append an `## Applied: F<n>` block
under that finding with: `files_changed`, `summary` (one sentence), and `deviation`
(or "none"). If a finding's required change is ambiguous or unsafe to apply
mechanically, append `## Blocked: F<n>` with a question instead — skip that finding,
continue with the rest.

Do NOT introduce features/refactors/stylistic changes beyond the named edits. Do NOT
touch paper.tex, notes/, or any prose document. After editing, RUN the affected scripts
(`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until
they exit 0 with all in-file checks passing, all within `timeout 600` (a timeout is a
failure — reformulate the math, never raise the cap). Getting the scripts to run cleanly
is your job; the orchestrator independently re-runs afterward.

---

## F1 — paper_misalignment / insufficient_assertion (the ξ_* bridge: radical constant + tolerance)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:86-92`
  (`xi_star_closed`, `xi_star`, `residual`, the `assert`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:84-86`
  (`xiStarClosed`, the `Print`, the `expectZero`)

**Issue (review R1, refined by consult Q2):**
The check `(1 - λ_{Π,0}) - ξ_* == 0` compares the convex-interpolation bias-neutral
point `1 - λ_{Π,0}` against the positive-family compensation closed form `ξ_*`. Per
consult Q2, this is KEPT as a **same-source closed-form consistency check** — both sides
descend from `gMinus`, and that is acceptable for this stage; it is NOT required to be an
independent bridge. Two real problems remain:

1. **Radical constant.** The superseded directive's prose/code/manifest/self-test carry
   `168π²` under the radical. That is WRONG. The LIVE scripts correctly use `100π²`, and
   `100` is MATHEMATICALLY FORCED, not a free choice:
   - `rF1 = sqrt(12·(37/20)²/π² - 1)` (py:22, wl:34);
   - `12·(37/20)² = 12·1369/400 = 4107/100`;
   - therefore `rF1² = (4107/100)/π² - 1 = (4107 - 100π²)/(100π²)`,
   so the `4107 - 100π²` and the `100π²` denominator come directly out of `rF1`. The
   `168π²` form in the old directive was a transcription error and is now PURGED here.
   **GUARD:** `100` is forced by `rF1`; do NOT substitute any other constant under the
   radical (in particular NOT `168`), and do NOT loosen the assert to mask a residual.
   (Unrelated: do NOT touch wl:36 `MaxIterations -> 100`, which is a FindRoot iteration
   cap, not this constant.)

2. **SymPy tolerance.** The SymPy `assert` (py:92) currently reads
   `abs(residual) < sp.Float("1e-15")` on a residual of ~`7.82e-17` — too loose to be a
   meaningful check (consult Q2). Mathematica already gets EXACT `0` via `FullSimplify`
   at WP80 (wl:86) — KEEP it as-is.

**Required change (SymPy py:86-92):**

Prefer the EXACT symbolic route (consult Q2: cleaner, does not need the nsolve float
`Pi_star`). `1 - λ_{Π,0}` is fixed by `gLam == gMinus` (py:78 in Mathematica;
SymPy gets `lam_Pi_zero` numerically at py:80). Build the EXACT bias-neutral point from
the symbolic `dPi_lam`/`g_lam` already in scope and compare symbolically:

```python
# Stage 126 positive-family compensation closed form (see notes section 3).
# GUARD: 100 under the radical is FORCED by rF1 (12*(37/20)^2 = 4107/100, so
# rF1^2 = (4107 - 100*pi^2)/(100*pi^2)); do NOT substitute any other constant.
xi_star_closed = (-37*sp.sqrt(3) - 5*sp.pi**2 + 2*sp.sqrt(4107 - 100*sp.pi**2)) / (5*(8 - sp.pi**2))
xi_star = sp.N(xi_star_closed, 30)
print("xi_* (Stage 126 closed form) =", xi_star)
# Preferred (EXACT, no Pi_star / no nsolve): 1 - lambda_(Pi,0) has the closed form
# (pi/4 - gminus)/(pi/4 - 2/pi), since the convex family is g_lam = (1-lam)*(2/pi)+lam*(pi/4)
# and lambda_(Pi,0) solves g_lam = g_* = gminus. Build it from the EXACT gminus (NOT the
# sp.N-numericized g_star/lam_Pi_zero, which would defeat the symbolic reduction). Note
# 1 + rF1^2 = 4107/(100 pi^2) so sqrt(1+rF1^2) = 37*sqrt(3)/(10 pi), and sqrt(4107)=37*sqrt(3),
# so SymPy collapses the nested radicals and reduces the difference to exact 0.
rF1_exact = sp.sqrt(12*sp.Rational(37, 20)**2/sp.pi**2 - 1)
gminus_exact = rF1_exact - sp.sqrt(1 + rF1_exact**2)/2
one_minus_lam_exact = (sp.pi/4 - gminus_exact) / (sp.pi/4 - 2/sp.pi)
exact_resid = sp.simplify(one_minus_lam_exact - xi_star_closed)
print("exact (1-lambda_(Pi,0)) - xi_* =", exact_resid)
residual = sp.N((1 - lam_Pi_zero) - xi_star, 30)   # numeric, for transcript display
print("(1-lambda_(Pi,0)) - xi_* =", residual)
assert exact_resid == 0, f"Stage 148 D4 consistency (exact) failed: residual = {exact_resid}"
```

If SymPy will NOT reduce `exact_resid` to a clean `0` (e.g. it leaves a nested-radical
expression), fall back to the high-precision NUMERIC route instead: raise the nsolve
working precision so the residual honestly drops below `1e-25`, then tighten the assert.
Concretely, at py:24 change `Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5), 40)` to
`Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5, prec=50), 50)` (so all downstream `sp.N(...,
40/30)` carry honest precision), and at py:90-92 use:

```python
residual = sp.N((1 - lam_Pi_zero) - xi_star, 50)
print("(1-lambda_(Pi,0)) - xi_* =", residual)
assert abs(residual) < sp.Float("1e-25", 50), f"Stage 148 D4 consistency failed: residual = {residual}"
```

The EXACT-symbolic route is PREFERRED; only use the prec-50 numeric fallback if SymPy
does not reduce the symbolic difference cleanly within `timeout 600`. Either route closes
the review's tolerance finding. Record which route you used in the `## Applied: F1` block.

**Required change (Mathematica wl:84-86):** purge any stale `168` — confirm the LIVE
line already reads `Sqrt[4107 - 100*Pi^2]` (it does) and leave the exact `expectZero`
(wl:86) UNCHANGED. No numeric edit; this is a no-op confirmation that the LIVE constant
is `100` and the check stays exact-zero via `FullSimplify`.

**Claim manifest:**

M1: `(1 - λ_{Π,0}) = ξ_*` where `λ_{Π,0}` is the convex-interpolation bias-neutral point
(`dPi_lam = 0`), and `ξ_* = (-37·√3 - 5·π² + 2·√(4107 - 100·π²))/(5·(8 - π²))` is the
positive-family compensation closed form. Both sides descend from `gMinus`
(same-source consistency check, consult Q2). `100` under the radical is forced by `rF1`.

**Verification command:**
`redteam exec-sympy 148` and `redteam exec-mathematica 148` — both must exit 0; the SymPy
transcript must show the new assertion (exact `0` residual, or `< 1e-25` numeric), the
Mathematica transcript must show `(1-lambda_(Pi,0)) - xi_* = 0` / PASS.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py`
- summary: Replaced the loose SymPy numeric xi residual assertion with the preferred exact symbolic consistency check while confirming Mathematica keeps the forced `100*pi^2` radical and exact-zero check unchanged.
- deviation: Used the exact symbolic route; Mathematica required no code change for F1 because the live file already used `Sqrt[4107 - 100*Pi^2]`.

---

## F2 — insufficient_verification (LIVE cross-engine divergence: the two engines compute DIFFERENT dT)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:43-48`
  (`dSigmaOfDeltas`, `dTOfDeltas`), `:53` (`dTU`), `:63` (`dTD`), `:73` (`dTLam`)
- cross-reference (do NOT edit) SymPy `AT`/`BT` at
  `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:32-41`

**Issue (review R2, root-caused by consult Q1 — DECIDED):**
The two engines disagree on the first-order T-correction. The LIVE transcripts show
SymPy `dTU = 0.5087563022150839...`, `dTD = -0.1169438021518107...`, while Mathematica
prints `dTU = 0.4976692633908835...`, `dTD = -0.1144451420715406...` (and the `dTLam`
slope differs likewise). **ROOT CAUSE (consult Q1):** Mathematica's
`dSigmaOfDeltas`/`dTOfDeltas` (wl:43-47) routes `dG` only through
`dPi = -dG/gPrimeStar` and then `dPi/(1 - sStar/4) + pStar·dS/(4·(1 - sStar/4)²)` — it
NEVER uses `sPrimeStar`, so it DROPS the S-follows-Π chain term that SymPy's `AT`
(py:32-37) carries through `Pi_star·Sp_star/(4·gp_star·(1-S_star/4)²)`. SymPy's `AT` is
the TOTAL derivative of `T_m` along the `S = Sformula(Π)` curve = the paper's `A_T`.
**SymPy is correct; the Mathematica `dSigmaOfDeltas`/`dTOfDeltas` path is the bug.**

Two engines passing while disagreeing is possible ONLY because nothing currently asserts
cross-engine agreement. Both problems are fixed together below.

**Required change (Mathematica wl:43-48):** make Mathematica compute its `aT`/`bT`
INDEPENDENTLY via its own symbolic differentiation (`D[]`) of `T_m` expressed as a
function of `Π`, with `S` substituted as the Π-function (mirroring how Stage 147 derived
`A_T`) — NOT a hand-port of SymPy's closed-form `AT` algebra. Replace the
`dSigmaOfDeltas`/`dTOfDeltas` block (wl:43-48) with:

```mathematica
(* aT, bT derived by Mathematica's OWN symbolic differentiation (NOT a port of SymPy AT). *)
(* T_m as a function of p with S = sFormula(p): the total p-derivative of T_m along the   *)
(* S-follows-Pi curve gives aT; bT is the explicit S-sensitivity at fixed p. Per consult  *)
(* Q1 this restores the sPrimeStar chain term the dSigmaOfDeltas route dropped.            *)
Tm[pp_] := Sqrt[(9/20)*(pp/(1 - sFormula/4))] /. p -> pp;
(* total dT_m/dPi along S = sFormula(Pi), projected onto dg via dPi = -dg/g'_*:           *)
aT = N[-(D[Tm[p], p] /. p -> pStar)/gPrimeStar, 30];
(* explicit S-sensitivity bT = dT_m/dS at fixed Pi = (9/(40 T_*)) * Pi_*/(4 (1 - S_*/4)^2) *)
bT = N[(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2), 30];
(* first-order correction in the paper's two-moment form: *)
dTOfDeltas[dG_, dS_] := N[aT*dG + bT*dS, 30];

(* --- Independence anchor 1: aT, bT vs the published appendix literals (EXTERNAL). --- *)
(* part04:846/848 -- these are paper values, NOT a SymPy mirror. *)
expectZero["aT vs paper literal A_T",
  If[Abs[aT - (-4.27263956256927)] < 10^-11, 0, aT - (-4.27263956256927)]];
expectZero["bT vs paper literal B_T",
  If[Abs[bT - 0.134875005736706] < 10^-11, 0, bT - 0.134875005736706]];
```

Note the sign/structure of `aT`: `D[Tm[p], p]` is the TOTAL p-derivative of
`Sqrt[(9/20)(p/(1 - sFormula/4))]`, so it already carries the `sFormula'` chain term;
dividing by `gPrimeStar` (with the `-` for `dPi = -dg/g'_*`) projects it onto `dg`. This
reproduces SymPy's `AT` WITHOUT typing SymPy's hand-expanded factors. If the Stage-147
`Tm[p_] := Sqrt[(9/20)*(p/(1 - sFormula[p]/4))]` form (with `sFormula` as a function)
matches the local script's conventions better, use that — the requirement is that `D[]`
regenerates the `(1 - S/4)²` and `sFormula'` factors itself.

**Required change (downstream uses):** the `dTU`/`dTD`/`dTLam` call sites already use
`dTOfDeltas[...]` (wl:53, wl:63, wl:73) and do NOT need to change — the redefined
`dTOfDeltas` now carries the correct `aT*dG + bT*dS`. Keep them as-is.

**Required change (SymPy py:32-41) — symmetric paper-literal anchor (NEW):** SymPy already
computes `AT`/`BT` via the chain-rule formula but does NOT assert them against the published
values. Add, immediately after `BT` (py:41), the EXTERNAL paper-literal asserts:

```python
# Paper-literal anchor (EXTERNAL, appendix part04:846/848): AT/BT reproduce the published
# A_T/B_T. A wrong AT formula (e.g. a dropped chain term) FAILS here.
assert abs(sp.N(AT, 30) - sp.Float("-4.27263956256927")) < sp.Float("1e-11"), f"AT vs paper A_T: {AT}"
assert abs(sp.N(BT, 30) - sp.Float("0.134875005736706")) < sp.Float("1e-11"), f"BT vs paper B_T: {BT}"
```

**Cross-engine agreement — how it is guaranteed (do NOT bake a cross-engine literal):** with
BOTH engines now anchored to the SAME external paper literals `A_T`/`B_T` (Mathematica via
"Independence anchor 1" above, SymPy via the asserts just added), the two derivations are
pinned to the same published truth and therefore agree — the `dTU`/`dTD`/`dTLam` values
(`= aT·dg + bT·dS`, with `g_u, g_d, S_u, S_d, g_*, S_*` computed identically in both engines)
follow. **Do NOT** hardcode SymPy's numeric `dTU`/`dTD` into the `.wl` as a literal: that
would (a) couple the two scripts via a fragile transcribed number and (b) re-introduce a
hardcoded-literal smell. The original divergence was ~`1e-2` (`0.4976` vs `0.5087`); the
`1e-11` paper-literal anchor on `AT` in BOTH engines catches a dropped chain term with
enormous margin. Final cross-engine numeric agreement of the printed `dTU`/`dTD` is
independently confirmed by the orchestrator's exec re-run + transcript diff (the reliability
gate), NOT by an in-`.wl` literal.

**Anti-transliteration / independence guard (consult Q1):**
- Mathematica's `aT`/`bT` MUST be derived by Mathematica's OWN `D[]` of `Tm[p]` (with `S`
  as the Π-function), NOT by copying SymPy's closed-form `AT`/`BT` expression. The
  paper-literal asserts (anchor 1, EXTERNAL) and the SymPy-agreement asserts (anchor 2)
  together are the independence witnesses: anchor 1 proves the value is the published one
  without copying SymPy's algebra; anchor 2 proves the two independent derivations
  converge.
- Consult Q1b CONFIRMED there is NO double-counting between `A_T` (which carries the
  `sFormula'`/`Sp_star` chain term) and `B_T·(S̄ - S_*)`: `A_T` is the total-derivative
  g-projection coefficient and `B_T` multiplies the full projected S-deviation in the
  SAME validated coordinate/projection model (already verified at Stage 147 / batch 6).
  This is the paper's validated projection model — do NOT treat it as a model error and
  do NOT "fix" it by removing the chain term.

**Verification command:**
`redteam exec-sympy 148` and `redteam exec-mathematica 148` — both must exit 0. SymPy
transcript must show `AT vs paper A_T` / `BT vs paper B_T` asserts passing (no AssertionError);
Mathematica transcript must show `aT vs paper literal A_T` PASS and `bT vs paper literal B_T`
PASS. Crucially, the printed Mathematica `dTU`/`dTD` must now MATCH the SymPy values
(`dTU ≈ 0.5087563...`, `dTD ≈ -0.1169438...`), NOT the old divergent `0.4976692...` /
`-0.1144451...` — the orchestrator confirms this cross-engine match by diffing the two
refreshed transcripts (reliability gate).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl`
- summary: Added SymPy paper-literal `AT`/`BT` anchors and changed Mathematica to derive `aT`/`bT` by symbolic differentiation of `T_m` along the `S(Pi)` curve with matching paper-literal checks.
- deviation: none

---

## Self-test notes (auditor)

- F1 identity to reduce:
  `(1 - λ_{Π,0}) == (-37√3 - 5π² + 2√(4107 - 100π²))/(5(8 - π²))`,
  with `λ_{Π,0}` the root of `dPi_lam = 0`, and `gMinus = rF1 - √(1+rF1²)/2`,
  `rF1 = √(12·(37/20)²/π² - 1)`. Because `12·(37/20)² = 4107/100`, the radical constant
  is `100` (NOT `168`). The SymPy side must reduce to exact `0` (preferred) or a residual
  `< 1e-25` at prec-50; Mathematica already reduces to exact `0` at WP80.
- F2: after the rework, BOTH engines anchor `AT`/`BT` (resp. `aT`/`bT`) to the paper
  literals `A_T ≈ -4.27263956256927`, `B_T ≈ 0.134875005736706` (part04:846/848), so the
  two independent derivations are pinned to the same published truth and the printed
  `dTU`/`dTD` agree across engines. Old divergent values `dTU = 0.4976692...`,
  `dTD = -0.1144451...` are the SIGNATURE of the dropped chain term — they must be GONE.
- Do NOT "fix" a failing paper-literal assert by editing the literal to match a wrong
  computed value, nor by re-introducing the chain-term-dropping `dSigmaOfDeltas`. Do NOT
  bake SymPy's numeric `dTU`/`dTD` into the `.wl`. A failure is a genuine signal — if it
  persists, append a `## Blocked: F2` question.
- Codex must not edit notes or paper. Out-of-scope label inconsistencies
  (Stage 250/127/126 references in those documents) are NOT this directive's concern.

---

## Already applied / no-op (bookkeeping from the superseded directive)

The superseded directive carried three sub-items that are now MOOT and require NO action:

- Its **F2 stub** ("subsumed by F1", the hardcoded `sp.Float("0.183918405511538")`): the
  LIVE SymPy already uses the symbolic `xi_star_closed` (py:87) — the hardcoded float is
  GONE. No-op.
- Its **F3(a) banner fix** (`STAGE 131` → `STAGE 148`): the LIVE Mathematica already reads
  `banner["STAGE 148 — ..."]` (wl:26). No-op.
- Its **F3(b) `dSigmaOfDeltas` restructure**: this is precisely the LIVE bug — it is
  CORRECTED here under F2, not preserved.

Net findings carried forward: **2** (F1, F2).

---

## RESOLVED (consult batch 7)

Both findings encode DECIDED outcomes from
`redteam/codex_reviews/_consult_batch7.md` (5/6 unconditional CONCUR; no disputes; no
change to any published paper-card conceptual claim):

- **Q1 (= F2, the live bug)** — CONCUR. SymPy is paper-consistent; Mathematica's
  `dSigmaOfDeltas`/`dTOfDeltas` (wl:43-47) drops the `sPrimeStar` S-follows-Π chain term.
  Fix: Mathematica derives `aT`/`bT` independently via `D[]` autodiff of `Tm[p]` (à la
  Stage 147); anchor `A_T`/`B_T` to the paper literals (part04:846/848) in BOTH engines
  (symmetric external anchor — guarantees mutual agreement without baking a cross-engine
  literal; orchestrator exec-diff confirms `dTU`/`dTD` now match). Q1b CONCUR: no
  double-counting between A_T (carries Sp) and B_T·(S̄-S_*) — validated projection model.
- **Q2 (= F1, the ξ_* bridge + tolerance)** — CONCUR (with scope). Keep the
  `(1-λ_{Π,0}) - ξ_*` check as a same-source closed-form consistency check (both sides
  from `gMinus`). `100π²` is CORRECT (`12·(37/20)² = 4107/100`); purge the stale `168π²`.
  Do NOT touch wl:36 `MaxIterations -> 100`. Tolerance: prefer the EXACT symbolic
  comparison (no nsolve float needed); acceptable fallback is `nsolve(..., prec=50)` so
  the residual honestly drops below `1e-25`. Either closes the review's tolerance finding.
