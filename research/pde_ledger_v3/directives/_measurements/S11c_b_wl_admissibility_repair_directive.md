# Measurements — S11c-b STEP 3 WL repair directive (v2, one item)

v2 has one item: the §3d full-field thickness-gradient content of the background-order admissibility operand.
v1's second item (thickness kinetic) was withdrawn by the decision-list legs — WL's `μ_W W_0²` is spec-correct
(§1a `e_W ≡ δW/W_0`); see `_measurements/S11c_b_wl_admissibility_repair_directive_review.md` and §9 of
`_measurements/S11c_b_step2_adjudication.md`. Every premise below carries its command (rule 2). Reads at HEAD
`fb0ff3ce`.

## 1. The disagreement that motivates the item (committed instrument run)
`~/.s11_build/S11c_b_row_residual_fixrun.out` (exit 0, 26.76 MB). `sed -n '363,372p' | cut -c1-400`
(LAB_HELD/RHO4, BODY_FORCE): the θ admissibility component differs cross-engine while the E_W component agrees
(`ROW_RESIDUAL = Integer(0)` for E_W). The v2 directive does not state either engine's θ value — it names the
§3d object and the structural non-uniformity in WL's construction.

## 2. WL's background-order admissibility construction (the structural fact the item names)
```
$ awk 'NR>=528 && NR<=577' mathematica/S11c_b_brane_operator_mathematica_audit.wl   # constructFullFieldBackgroundEnergy
$ sed -n '1329,1351p' mathematica/S11c_b_brane_operator_mathematica_audit.wl        # backgroundBalanceFromModel
```
- L540 `localEw = (WZero/anchoredWidth) ewVariation` — thickness **perturbation**, order-scaled.
- L541 `fullWidth = anchoredWidth + WZero ewVariation` — the **full** thickness field.
- L543 `gradLocalEw = gradient[localEw]`.
- L554 invariant [[7]] `anchoredWidth^(-2) Dot[gradient[fullWidth], gradient[fullWidth]]` — the pure-thickness
  gradient invariant, built on the **full** field `fullWidth` (this sources the agreed E_W body force).
- L555 invariant [[8]] `Dot[gradTheta, gradLocalEw]` — the **mixed** thickness-gradient invariant, built on
  the **perturbation** `gradLocalEw`, not the full-field thickness gradient.
- L1334–1335 `firstVariation = D[fullFieldEnergy, backgroundOrder] /. backgroundOrder -> 0`; L1340
  `eulerScalar[firstVariation, thetaField]`.

⇒ §3d's full-thickness-field gradient `∇(W_bg+δW)` is applied to [[7]] but not to the mixed [[8]]. (Rule-13
orchestrator derivation, matching both decision-list legs: the mixed invariant on the perturbation-only
`gradLocalEw` is `order²` in the perturbations, so its background-order first variation vanishes; on the
full-field thickness gradient its `order⁰` background piece makes it linear-in-θ and it survives — the
mechanism the item repairs.)

## 3. Spec grounding
- §2b (L209–235), §3d (L325–356), §3a (L242–270), §1c (L93–149). §3d quote fixing the gradient content:
  "Its thickness and variable-coefficient gradient content must be the gradient content of the full fields,
  including `∇(W_bg+δW)`…". §3a names `∇θ·∇e_W` among the independent invariants and retains a second spatial
  derivative of a background profile at first background-amplitude order. §3d/§3a prohibit a `ρ_4D` density
  multiplier and the `W₀→W_bg` substitution into the uniform `U`.

## 4. Leak gate (v2) — grep + semantic read
```
$ D=directives/S11c_b_wl_admissibility_repair_directive.md
$ grep -n 'kappa_theta\|K_thetaW\|w1_d1d1\|covector\|Integer(0)\|∇²' $D      # → no matches
$ grep -ni 'do not freeze\|currently zero\|should be\|must be nonzero\|is zero' $D  # → only L26 (the withholding clause: "currently zero or nonzero. Those are outputs.")
```
- No withheld value (θ body-force coefficient/form/sign) appears.
- The direction is explicitly withheld: L26 states the directive "does not tell you whether the object is
  currently zero or nonzero. Those are outputs."
- Semantic read of every `W_bg`/`δW`/`full-field` line: each is either the §3d spec quote naming the object
  (`∇(W_bg+δW)`), a §3d/§3a prohibition, the §1a kinetic scope note (`μ_W W_0²`, a term the build leaves
  byte-identical), or the structural code-fact (L541/L554 vs L555/L543). None states the θ body-force value.
- The v1 grep-only gate was incomplete (it missed negation-bounds); v2 adds the semantic read above per the
  decision-list legs' finding.
