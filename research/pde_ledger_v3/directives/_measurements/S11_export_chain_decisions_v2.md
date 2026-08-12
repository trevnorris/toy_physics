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
