---
unit_id: 111
batch: IV.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 111

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
A new independent re-derivation block was inserted into the Mathematica script at
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl:54-62`:

```
(* Independent re-derivation of chi_Q^mix: bypass the L0/L5 extraction *)
(* and compute directly from the geometric-series form of the pole. *)
poleSeries = Series[sigma/(1 - kappa*z^2 - I*gamma*z^5), {z, 0, 5}];
imagPart5  = Coefficient[Normal[poleSeries], z, 5]/I;
chiMixAlt  = FullSimplify[27*(1/9 - imagPart5)/(3 + sigma),
                          Assumptions -> $Assumptions];
Print["chi_Q^mix (independent route) = ", fmt[chiMixAlt]];
```

A fifth assertion was added at L75:
`expectZero["chi_Q^mix routes agree", chiMix - chiMixAlt];`

The existing four `expectZero` calls (L71-74), the existing `chiMix` definition (L52),
and the existing `lambdaMix`/`l0..l5`/`kappaMatch`/`sigmaMatch` pipeline were preserved
exactly as before. The directive's required block matches the inserted text verbatim.

**Assessment:**
The new route is algebraically independent of the SymPy choreography in three concrete
ways: (i) `Series` is applied to the pole `sigma/(1 - kappa z^2 - I gamma z^5)` alone,
not to `lambdaMix = lambdaOut - pole`; (ii) the `1/9` contribution from `Lambda_out`'s
imaginary `z^5` coefficient is hard-coded rather than extracted via `Coefficient[lambdaMix, z, 5]/I`;
(iii) the `-L0 = 3 + sigma` denominator is hard-coded rather than extracted via
`Coefficient[lambdaMix, z, 0]`. So a hypothetical bug in the `lambdaMix` construction
or in the `l0`/`l5` extraction in the original pipeline would not propagate into
`chiMixAlt`. The new `expectZero` is non-tautological — it is the residual of two
algebraically distinct derivations of the same closed form, and it would catch a
choreography error in either pipeline. Mathematica transcript shows both routes
yield `(3 - 27*gamma*sigma)/(3 + sigma)` and the residual is `0`.

Cluster C banner sweep cross-check:
- `.py` L36: `print('stage111: PASS')` — confirmed (was `stage94: PASS`).
- `.wl` L26: `banner["STAGE 111 — MIXED SIDE-CHANNEL POLE"];` — confirmed (was `STAGE 094`).

No collateral edits beyond the directive scope (the existing four assertions, the
SymPy script's algebraic body, and the four pre-existing chi_mix/kappa_match/sigma_match
identities are untouched).

Side note: the directive specified a `## Applied: F1` block be appended to
`redteam/directives/stage_111.md`. That block is absent in the directive file as
of verification time. The code edits themselves match the directive exactly and
both engines pass, so this is a logging gap rather than an edit-correctness issue;
flagging it as a side observation rather than a verification blocker.

## Exec log assessment

**SymPy:** exit=0 (script reaches `print('stage111: PASS')` at L36; all four
`assert` statements at L32-35 evaluated cleanly). Notable lines from
`scripts/output/moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.txt`:

```
kappa from z^2 matching = -1/9
sigma from z^4 matching = 0
chi_Q^mix = 3*(-9*gamma*sigma + 1)/(sigma + 3)
chi_Q^mix linearized = -9*gamma*sigma - sigma/3 + 1
stage111: PASS
```

**Mathematica:** exit=0. All five `expectZero` calls printed `PASS:`. Notable lines
from `mathematica/output/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.txt`:

```
chi_Q^mix (independent route) = (3 - 27*gamma*sigma)/(3 + sigma)
chi_Q^mix = (3 - 27*gamma*sigma)/(3 + sigma)
PASS: kappa_match + 1/9
PASS: sigma_match
PASS: chi_Q^mix - 3 (1 - 9 sigma gamma)/(3 + sigma)
PASS: chi_Q^mix linearized - (1 - sigma (1/3 + 9 gamma))
PASS: chi_Q^mix routes agree
```

The new "chi_Q^mix routes agree" assertion and the matching independent-route print
line are both present and pass.

**Output freshness:** both `.txt` outputs are newer than their source scripts.
- `.py`: mtime 15:14, output mtime 15:18 (post-edit).
- `.wl`: mtime 15:14, output mtime 15:24 (post-edit).
Confirmed re-generated post-fix.

## Material-change assessment

`material_change`: false.

No derived numerical or symbolic result changed. The original four identities
(`kappa_match = -1/9`, `sigma_match = 0`, `chi_Q^mix = 3(1 - 9 sigma gamma)/(3 + sigma)`,
linearization) are unchanged in both engines. The new block is a strict addition:
a redundant cross-check of `chi_Q^mix`, plus a banner string update. Downstream
units cannot depend on any newly-published quantity because no new quantity was
published; `chiMixAlt` only appears inside the new `expectZero` residual.

## Side observations (non-blocking)

1. The Codex directive boilerplate asks for a `## Applied: F1` block to be appended
   to `redteam/directives/stage_111.md`. That block is absent. The code edits
   themselves are correct and complete, so this is a metadata/logging gap, not a
   correctness issue. The orchestrator may want to add the missing block manually
   or relax the requirement.
2. The new `imagPart5` quantity is named ambiguously — it represents the imaginary
   part of the z^5 coefficient of the pole alone, divided by `I`, i.e. a real
   symbolic scalar `gamma*sigma`. Cosmetic; does not affect correctness.

## Verdict justification

The single low-severity finding F1 (mathematica_transliteration) is closed. Codex
inserted the directive's specified independent re-derivation block verbatim,
added the fifth `expectZero` comparing the two routes, preserved all four
pre-existing assertions, and applied the Cluster C banner sweep on both `.py` L36
and `.wl` L26. Both engines exit 0; the SymPy transcript shows the four-assertion
PASS line and the Mathematica transcript shows all five PASS lines including the
new "chi_Q^mix routes agree" with residual 0. The new check is non-tautological
because it computes `chi_Q^mix` via a path that avoids both the `lambdaMix`
construction and the `l0`/`l5` Coefficient extraction. Verdict: `verified`.
