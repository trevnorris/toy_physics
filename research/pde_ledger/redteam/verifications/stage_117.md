---
unit_id: 117
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 117

Status-consolidation card; four findings, all resolved via Claude+Codex consult `019e74f7`
(all CONCUR). Verified the post-fix scripts under their own assertions, plus an adversarial
structural read of each finding. Both engines exit 0; outputs are fresh.

PASS-line count per engine: **SymPy 12 `expect_zero` assertions** (the SymPy log shows 13
lines ending `= 0`, but one — `standalone mixed-pole sigma match = 0` at log L22 — is an
informational `print` of an intermediate value, not an assertion; the real assertion
`standalone mixed pole disappears` is the next line). **Mathematica 12 `PASS:` lines.**
The two counts reconcile exactly (12 ↔ 12).

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved (accepted policy-mirror; no script change for F1)

**What changed:** Nothing for F1. Per the consult and the Mathematica mirror policy, the
`.wl` is accepted as a sanctioned cross-engine port of a pure rational series-coefficient
classification, not as an independent-route second engine.

**Assessment:** The `.wl` diff (patch L1–54) contains ONLY the F2 substitution
(`lWForward`/`kappa0FromTube`, deltaCore rule change) and the F3 comment-mirror edit. No
independent-route algebra was added — `rC`, `rhoC`, `sigmaC`, `kappaC`, `gammaC`,
`gqSolutions`, `deltaCoreExpected` are unchanged from the pre-fix mirror. This is exactly the
agreed disposition. Resolved as recorded.

### F2 — tautological_check → upstream-import de-tautologization (THE CRUX)

**Classification:** resolved

**What changed:**
- SymPy (`scripts/...stage117..._sympy_audit.py:177-188`): the target-inverting
  `Lw_required = sp.solve(sp.Eq(4*Lw**2/(pi**2*a**2),(1+r_c)/3), Lw)[0]` is GONE. Replaced by
  the forward closed form `L_W_forward = sp.pi*a*sp.sqrt((1+r_c)/3)/2` then
  `kappa0_from_tube = sp.simplify(4*L_W_forward**2/(sp.pi**2*a**2))`. The `delta_core` subs
  dict now uses `kappa0: kappa0_from_tube` (was hardcoded `(1+r_c)/3`); `gamma0: (1+r_c)/9`
  retained.
- Mathematica (`mathematica/...stage117..._mathematica_audit.wl:128-135`): mirror — `lWRequired`
  Solve-inversion GONE; `lWForward = Pi a Sqrt[(1+rC)/3]/2`,
  `kappa0FromTube = FullSimplify[4 lWForward^2/(Pi^2 a^2)]`; `deltaCore` rule now
  `kappa0 -> kappa0FromTube`.

**Assessment (adversarial):**
1. CONFIRMED the inverted-target value is no longer fed into `delta_core`. `grep` for
   `Lw_required`/`lWRequired` returns NO MATCH in either script — the X−X pitfall variable was
   fully removed, not merely bypassed.
2. CONFIRMED `kappa0_from_tube`/`kappa0FromTube` is built from the explicitly-written forward
   `L_W` closed form in BOTH engines (sqrt prefactor `pi*a/2`, argument `(1+r_c)/3`,
   exponent 2 on `L_W`).
3. Adversarial coefficient test: the residual `delta_core - delta_core_expected` (SymPy L190 /
   `.wl` L137-139) is the load-bearing assertion. If the law were written `sqrt((1+r_c)/2)`,
   `kappa0_from_tube` would equal `(1+r_c)/2` → `kappa_c = 1/2 ≠ 1/3`, and the O(z²) term of
   `delta_core` would no longer match `delta_core_expected`'s `z²/3` pole structure, so the
   series residual would be nonzero and the assertion would FAIL. Same for a wrong prefactor
   (`pi*a` instead of `pi*a/2` gives `kappa0_from_tube = (1+r_c)` → `kappa_c = 1`). So the
   check genuinely EXERCISES (and could falsify) the in-script tube-length coefficient.

HONEST strength assessment: this ties `κ_c = 1/3` to the L_W closed form *written in-script*
— a consistency check that fails on a wrong/typo'd coefficient. It is NOT a from-scratch
re-derivation of the D/N half-wave BVP eigenvalue `k_W = pi/(2 L_W)` (that lives upstream at
stage 116; stage 117 does not re-solve the BVP). For a status-consolidation card whose declared
job is to carry forward upstream results, that is an acceptable de-tautologization: it converts
a pure `1/3 → 1/3` identity into a falsifiable check of the forward law's coefficient. Verdict:
acceptable, not still-tautological. (A strictly stronger check would re-solve the D/N
eigenproblem in-stage; that is stage 116's role and out of scope for a consolidation card.)

### F3 — tautological_check → comment-only provenance fix

**Classification:** resolved

**What changed:** Section-5 comment block (SymPy L141-155 / `.wl` L98-108) rewritten to label
κ₀ FORWARD-DERIVED at stage 116 and γ₀ a pure-scale ANSATZ citing the stage-116 note. The
γ₀ status print (SymPy L180 / `.wl` L131) now reads "...pure-scale ANSATZ ... not derived ...
consistency-of-assumption check".

**Assessment:** (a) CONFIRMED the false "Stage 119"/"derived upstream" attribution for γ₀ is
gone — `grep "119"` returns NO MATCH in both scripts AND both regenerated outputs. (b) γ₀ now
explicitly labeled an ansatz citing the stage-116 note "Bare outgoing normalization". (c) κ₀'s
"FORWARD-DERIVED ... at Stage 116" attribution is retained (genuinely derived upstream).
(d) CONFIRMED no assertion logic touched: the diff under F3 is comments + two `print`/`Print`
strings only; `kappa_c`/`gamma_c` assignment lines (SymPy L159-160 / `.wl` L112-113) and all
`expect_zero`/`expectZero` calls are byte-unchanged. The `gamma_c = 1/9` consistency check is
unchanged. Comment-only as required.

### F4 — tautological_check → reconcile already-applied

**Classification:** resolved

**What changed:** No diff under F4 (flags were wired by a prior pass; F4 only reconciled +
re-ran). I read the live scripts to confirm the wiring stands.

**Assessment:** CONFIRMED the survivability flags are computed from sections 1-5, not hardcoded
`True`:
- `even_ok_*`/`odd_ok_*` reference genuinely-upstream-computed variables — `beta_solutions`
  (L58), `robin_even_solutions` (L74), `kappa_match` (L90), `sigma_match` (L91), `chi_mix`
  (L95), `hybrid_solutions` (L108), `chi_cancel` (L120), `chi_comp` (L121), all confirmed to
  exist via grep.
- `nontrivial_scale/robin/hyb_cancel = False` are correct physics (trivial classes), not
  hardcoded survivability passes; `nontrivial_standalone = (simplify(sigma_match) != 0)` is
  computed.
- `nontrivial_compensated` (SymPy L214-216) is tied to the section-5 series residual
  `delta_core - delta_core_expected == 0`; `.wl` L163-165 mirrors this.
- The capstone (SymPy L228-230 / `.wl` L175-178) prints EXACTLY ONE nontrivial survivor —
  `compensated Robin-mixed core realization` — in both outputs (sympy log L58, mathematica log
  L70).

F2 consistency: since `kappa0_from_tube` simplifies to `(1+r_c)/3`, the residual stays 0 and
`nontrivial_compensated`/`nontrivialCompensated` stays `True` — confirmed by both fresh logs
showing the residual `= 0` / `PASS` and the single-survivor capstone passing. Consistent.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L47 `carrying forward (Stage 116): L_W = pi a sqrt((1+r_c)/3)/2 -> kappa_0_bare = 4 L_W^2/(pi^2 a^2) -> kappa_c = 1/3` (forward-law print, F2)
- L48 `gamma_0_bare = (1+r_c)/9 is a pure-scale ANSATZ ... not derived ...` (F3; no "Stage 119")
- L49 `concrete core collapses to the compensated hybrid class = 0` (load-bearing residual, F2 path)
- L58 capstone single survivor `('compensated Robin-mixed core realization', True, True, True, ...)`

**Mathematica:** exit=0. Notable lines:
- L58 `carrying forward (Stage 116): L_W = Pi a Sqrt[(1+r_c)/3]/2 -> ... -> kappa_c = 1/3`
- L59 `gamma_0_bare = (1+r_c)/9 is a pure-scale ANSATZ ...`
- L61 `PASS: concrete core collapses to the compensated hybrid class`
- L70 capstone single survivor; L77 `Stage 117 Mathematica audit passed.`

**Output freshness:** CONFIRMED. Both `.txt` outputs regenerated 2026-05-29 13:25:44, newer
than both scripts (`.py` 13:23:26, `.wl` 13:20:34). Logs dated 13:24:37 / 13:25:25.

## Diff-scope check

Full patch is 124 lines; `git diff --stat` = 2 files, 48 insertions / 28 deletions. I scanned
every hunk:
- `.wl`: two hunks — (1) F3 comment block rewrite at L95-110; (2) F2 substitution
  (`lWForward`/`kappa0FromTube`, two `Print` strings, `deltaCore` rule) at L120-138.
- `.py`: three hunks — (1) F3 comment block at L138-155; (2) F2 substitution at L161-188;
  (3) removal of a single trailing blank line at EOF (patch L120-124). The trailing-blank-line
  deletion is cosmetic and harmless.
No stray or unauthorized edits; the diff is exactly the F2 substitution + F3 comment edits
(plus the trailing-newline trim). Engines use genuinely different mechanics (`sp.solve` /
`sp.series` vs `Solve`/`Series`, distinct variable-naming conventions) — accepted as a
policy-mirror per F1, not blind transliteration of a bug.

## Material-change assessment

`material_change`: **false**. F2 reroutes an already-true substitution through its forward
source; the substituted VALUE is unchanged (`kappa0_from_tube` simplifies to `(1+r_c)/3`,
identical to the old hardcode). Every printed residual, the `delta_core` collapse, the
capstone survivor set, and `kappa_c = 1/3` / `gamma_c = 1/9` are bit-for-bit the same as
pre-fix. No carried result moved; no downstream unit's inputs change. (F3/F4 are
comment/already-wired, no value change.)

## Side observations (non-blocking)

- The banner/heading numbering quirks the auditor noted (script banner historically "STAGE 100",
  paper section "Stage 134" vs `\label{stage:117}`) are cosmetic and outside scripts-math scope;
  the current banner correctly reads "STAGE 117". Not a verification blocker.
- Two paper "Checks" items (Schur-complement positivity, parent-overlap ratios) remain
  untested in-script per the original audit; these were not findings in this directive and are
  out of scope for this verification.

## Verdict justification

All four findings resolve. F2 (the crux) genuinely removes the X−X inversion (`Lw_required`
fully deleted) and routes `delta_core` through the in-script forward tube-length closed form in
both engines; a wrong coefficient would break the load-bearing series residual, so the check is
falsifiable rather than tautological — acceptable for a consolidation card, with the honest
caveat that it is a consistency check against the written law, not an in-stage BVP
re-derivation. F3 is a clean comment-only provenance fix (γ₀ relabeled ansatz, "Stage 119"
gone everywhere, κ₀ attribution and all assertion logic untouched). F4's flags are computed
from real section-1–5 residuals and the capstone yields exactly one survivor, consistent with
F2's unchanged value. Diff scope is exactly F2+F3 (plus a harmless trailing-newline trim); both
engines exit 0 with fresh outputs. `material_change: false`.
