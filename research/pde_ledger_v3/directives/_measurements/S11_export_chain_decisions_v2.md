# Measurements behind `S11_export_chain_decisions_v2.md`

Commands and their literal output. Rule 2: a claim about an artifact carries the command that produced
it. Regenerate with the commands as written; do not transcribe.

## `B_comp > 0` is a declared premise — kills F9 round 1's supersession branch

```
$ grep -n "B_comp" directives/S11_SHARED_PHYSICS.md | grep -iE "> 0|posit|premise|force"
69:`μ_R > 0`, `B_comp > 0`, `μ_br > 0` and `β` are real constants. ⚠ **`β` carries no sign premise**; §7 says
119:ρ_br > 0 ,  μ_R > 0 ,  B_comp > 0 ,  Σ_m k_m² > 0 ,  every k_m real ,  every a_j real ,
1027:⚠ **`XCOEF_BSIGN` flips a sign in the action, ⛔ not a premise.** `B_comp > 0` stays in force; the minus
1117:| `ρ_br > 0`, `μ_R > 0`, `B_comp > 0`, `Σ k_m² > 0` | |
```

## Most of the imported namespace is not an expression — F9's totality requirement

```
$ python3 <census of payload types on which a subtraction residual raises>
rows total                 : 617
subtraction residual OK    : 129
subtraction residual RAISES: 488
     368  Tuple
      51  BooleanTrue
      24  Str
      17  And
      16  Equality
      12  AppliedPredicate
```

## The class vocabulary is not in the shared spec — the citation defect both legs found

```
$ grep -c "KNOB" directives/S11_SHARED_PHYSICS.md
0
$ grep -rn "Classes: .KNOB\|KNOB . STRUCTURAL" --include=*.md .
./directives/_measurements/S11_export_chain_decisions_v2.md:36:$ grep -rn "Classes: .KNOB\|KNOB . STRUCTURAL" --include=*.md .
./directives/S10_export_chain_decisions.md:35:   `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`. ⛔ An annotation never states what the
./directives/S9_export_structural_review_prompt.md:36:   class vocabulary `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`. ⚠ If some object is
./directives/S9_structural_review_round2_prompt.md:47:   `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`, and say where a tag is being used as a
./REBUILD_HANDOFF.md:40:   `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`.
./S9_REWRITE_PLAN.md:41:⭐ Classes: `KNOB` · `STRUCTURAL` · `COORDINATE` · `CONTROL` · `PREMISE` · `DERIVED`.
```

## The controls the OWED checklist protects — all in committed S10 code

```
$ grep -n 'class_residual = \|duplicate S10 export key\|def exact_reconstruction_match\|!= "S9"' \
    scripts/S10_brane_mode_spectrum_sympy_audit.py
1906:def exact_reconstruction_match(live_value: object, reconstructed_value: object) -> bool:
2073:            raise RuntimeError(f"duplicate S10 export key: {name}")
2085:        class_residual = sp.Integer(upstream["class"] != record.class_tag)
2087:            not record.overwrites_upstream or upstream["step"] != "S9"
```

---

## F10 added (2026-09-03) — the model-level-register export rule

F10 resolves the "model-level vs step-level split" that this file's "What F9 does NOT decide" left open
"for a third step's evidence." That third-step evidence is S11c-c1's, committed at d1fe1bf0
(_measurements/S11c_c1_sympy_build_directive.md): the pre-trim full-primary export was 152574695 bytes
(145.5 MB), of which two dissipation-AUDIT rows were 84.7 MB (56%), exceeding GitHub's 100 MB plain-git cap
(and an *_exports.py must stay plain-git, never annexed).

F10: the plain-git LEDGER carries the chain-CONSUMED model-definition set (the emitting spec's declared
consume-set) plus its free-symbol closure; step-level diagnostics (dissipation/energy/loci/structural-views)
are emit-only, in the annexed .out. It refines D1 by a SPEC-declared boundary (not builder judgement), so
D1's anti-under-export default stands. Safe because the cross-engine comparator, the review legs, and the WL
engine read the .out tag streams, not the LEDGER (F9 is PY-only; WL writes no ledger). Gated by the S11c-c1
build-directive review (Grok CLEAR + Codex must-fix, saved under directives/_legs/); ratified by the user
2026-09-03.

---

## F10 SUPERSEDED (2026-09-03) — by the bind-closure design c04e071f

A two-engine deliberation + two-leg design review found F10 insufficient: it used category membership (not the
bind test), carried the whole inherited LEDGER forward (leaving the 61.6 MB dump), wrongly claimed D3 catches
the closure gap, and kept an export-every-primary fallback. F10's correct half (the EMIT-vs-EXPORT channel
split) survives in the design's §D1/§D3. F10 is marked SUPERSEDED in place (retained as a record); the
"model-level vs step-level" open note now points at the design's bind test. Authority:
export_ledger_bind_closure_design.md (committed c04e071f). Gate: the deliberation + design-review leg outputs
under directives/_legs/export_ledger_*.
