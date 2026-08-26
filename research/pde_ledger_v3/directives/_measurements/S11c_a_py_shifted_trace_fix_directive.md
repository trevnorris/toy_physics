# S11c-a PY shifted-trace fix — measurement record (rule 2)

The commands and their literal results behind the fix. Interpretation is in the step record; this file is
the evidence.

## The finding (adjudicated)
T7 cross-engine reconciliation: the two blind engines agreed on the 24 pure-geometry shape-derivative cases
and diverged on every traced-bulk one. Root cause (adjudicated): the SymPy (PY) engine mis-applied §3c —
(D1) it injected free background-normal-jet PREMISE symbols the spec does not supply
(`d_w_v_bulk_0`, `d_w_delta_p_0`, `d_w_j_0`, `d_w_rho_4D_0`), producing a spurious term-2; and (D2) it froze
the traced perturbation at the flat face `sW₀/2`, dropping the mandated ε·η term-1 correction. WL was correct.

## Adjudication legs (independent derivation from spec alone; both residual-0 CAS)
- Fresh Claude agent + Grok, each composed `f=f⁰+ε δf` at `h_s=sW_bg/2+ε δh`, expanded to first order in
  ε and η, residual 0. Both concluded: for the bulk velocity/pressure/current/density the background is
  zero/`w`-independent ⇒ term-2 drops; the term-1 η-correction is mandated (§3c:382-383, §2a:195-198).
  Verdict: keep-only-η-correction ⇒ PY wrong on both, WL correct. My own spec/source check agreed
  (spec supplies no background bulk-velocity profile; PY line 87 imports `v_dr` but never uses it, and
  fabricates `d_w_v_bulk_0` at line 157).

## The fix (Codex build; directive `S11c_a_py_shifted_trace_fix_directive.md`, rev 2 folded from 2 legs)
Removed the four fabricated background-jet premises; every background normal jet now comes from
`sp.diff(supplied_background, w)` (line ~597) and is genuinely zero; the traced perturbation is evaluated at
`h_s⁰=sW_bg/2` (line ~844/601). §3c clarified in `S11c_a_SHARED_PHYSICS.md` (structural premises only).

## Fix-review legs (Codex-written ⇒ fresh Claude agent + Grok; both PASS with mandatory FORM ablation)
Each independently derived (residual 0), then ablated a /tmp copy:
- Ablation A (perturbation eval face `sW_bg/2`→flat `sW₀/2`): velocity/pressure/density lose the η-correction
  → load-bearing. Ablation B (a supplied background made nonzero/`w`-dependent): term-2 `δh·∂_w f⁰` reappears
  and matches exactly, and vanishes under the true zero/`w`-independent background → the machinery DERIVES
  term-2 from the real background, not deletes it. Conormal generic probe correctly retained. No assert-before-emit.
- Grok artifacts: `research/pde_ledger_v3/_measurements/s11c_a_sympy_shifted_trace_review_grok/`.

## Cross-engine reconciliation of the FIXED engine (the confirmation)
Run fixed engine → `~/.s11_build/S11c_a_py_fixed_run.out`. TWO reconciliation commands (distinct outputs):

(1) `python3 ~/.s11_build/S11c_a_reconcile_fixed.py` (measure_reconcile over the fixed PY `.out` + the
mechanical `d_w_X`↔`X_dw` perturbation-jet rename). Columns = TAG join zero nonz unmap err. Literal:
```
FACE_VELOCITY                     8    8    0     0   0
FACE_NORMAL                       8    8    0     0   0
FACE_MEASURE_SHAPE_DERIV          8    8    0     0   0
RELATIVE_FLUX                     8    8    0     0   0      (WAS 0 zero / 8 nonz before the fix)
TRACTION                         16    0   16     0   0      (nonz is naming/CAS-form only — resolved in (2))
EVOLUTION_MASS_BALANCE            8    0    8     0   0      (comparator single-face gap — see OPEN below)
```
(2) `python3 ~/.s11_build/S11c_a_reconcile_traction_check.py` (adds the mechanical `mu_theta_L/M→mu_theta`
rename + `sp.cancel(sp.together(·))` for the λ_X `1/(1−iωτ_X)` complex-denominator factoring — both
CAS-form/naming, NOT physics). Literal:
```
RELATIVE_FLUX            join=8 zero=8 nonzero=0
TRACTION                 join=16 zero=16 nonzero=0
```
⇒ Both traced-bulk primaries that carried the finding reconcile EXACTLY once the naming/CAS-form folds are
applied. The fixed PY agrees with WL. (The 24 pure-geometry cases in (1) stay at residual 0.)
⇒ Both traced-bulk primaries that carried the finding now reconcile EXACTLY. The fixed PY agrees with WL.

## OPEN (harness, NOT a finding, NOT a fix defect): EVOLUTION_MASS_BALANCE reconciliation
EVOLUTION still shows a residual under the instrument. Classified: the raw WL EVOLUTION
`SHAPE_DERIVATIVE.EXPRESSION` is complete (2528 chars, both faces: 8×`{-1/2*W0}` + 8×`{1/2*W0}`,
brackets balanced 19/19,98/98,24/24), but EVOLUTION's case key has NO single FACE axis, so
`canon_wl(expr, face)` applies ONE face to a two-face-summed expression and collapses the minus-face field
evaluations. Fix belongs in the comparator (detect the face from the `{±1/2*W0}` evaluation-point argument),
not in the engine. Recorded for the comparator-hardening step.
