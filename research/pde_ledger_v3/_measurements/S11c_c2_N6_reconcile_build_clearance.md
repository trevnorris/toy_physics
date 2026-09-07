# S11c-c2 N6 RECONCILE instrument build — CLEARED (both build legs) (2026-09-06)

Instrument `scripts/S11c_c2_N6_reconcile_sympy.py` (astra-built per the CLEAR-TO-BUILD directive `08d72d46`; baseline
committed `36d59b95`). Two build legs (astra/Codex-written → **fresh Claude agent + Grok**), identical prompt
`directives/_legs/S11c_c2_N6_reconcile_build_review_prompt.md`. **Both legs: BUILD CLEAR.** Reports: the fresh-Claude
agent transcript + `scratchpad/grok_N6rc_buildreview.log`. Convergent, independently-derived + ablated evidence.

## Both legs, convergent (independent derivation + FORM ablation + one-sided bites + PIT + AST)
- **Affine identity + three-way split verified on the REAL objects** (each leg's OWN sampler + OWN sig-{6,9,12} `B`,
  not the artifact's `closed_response`): `I(C,S)=−C·p+B(C,S)` (sig-0 is source-independent), and
  `R_N6 = I(ΔC,ms) + B(C_M,ΔS) + B(ΔC,ΔS)` exact (0 nonzero split residual). The artifact's `closed_response` == the
  sig-{6,9,12} restriction of `build_increment` (0 differ). The forbidden `build_increment(C_M,ΔS)` re-adds 12 sig-0
  columns — the artifact does NOT do that.
- **FORM ablation (mandatory): SPLIT_CHECK is load-bearing** — forcing SOURCE through `build_increment` (spurious
  sig-0) moves SPLIT_CHECK 0→12; swapping `ms`→`es` in SOURCE moves it 0→18. On the unmodified instrument SPLIT_CHECK
  is samplewise-zero across the full PIT (288/288, resp. 480/480) on genuine `add`/`mul` nodes (not `number(0)`; the
  structurally-empty columns are zero on both sides).
- **One-sided bites clean:** source bite (`a_ρ: t→0`, RHOBR) moves SOURCE_BRIDGE/SOURCE_CHANNEL/MATERIAL_OPERAND/R_N6,
  leaves CARRIER_*/C_M/es/EULERIAN byte-identical (⇒ `C_M` is μ-independent, routes independent); carrier bite (mix the
  material `normal_exact`, `ms` held) moves CARRIER_BRIDGE(0→4)/CARRIER_CHANNEL/C_M/R_N6, leaves
  C_E/SOURCE_BRIDGE/es/EULERIAN byte-identical. RHO4 emits computed absence (`a_ρ=0`, `∂μ/∂t=0` via `sp.diff`), ⛔ not
  an A−A.
- **PIT sound:** shared samples (one `Sample` per (prime,cell,draw), all roots); 3 primes ≡1 mod4 isprime-checked;
  draws 8→10 by degree bound; joint singular rejection; on-shell `q` from `k` (all dispersion numerators 0); honest FN
  bound `family·max(per_prime)` worst-per-prime (δ≈3.97e-20), `conditional_good_prime=True`; residual-zero never an
  exit. **0 `assert`s.** Every stage built from `C_E,C_M,es,ms`, ⛔ never from `R_N6`; tag names name objects not
  values; frozen relations (incl. jet bridge `theta_d1↔grad_theta_1`, fixed-anchoring `h_α=0` at MATERIAL_ADVECTED)
  predeclared + emitted before any residual; `es`/`ms` are the diagnostic's `source_terms` circuits, ⛔ not a recoded Λ;
  carrier bridge is coefficient-level.

## Two caveats — both NON-BLOCKING (do not change what the instrument computes or may be claimed)
1. **Identity-map carrier bite is false-negative-prone (BOTH legs).** Corrupting `material_inverse_transpose → eye`
   (the target the directive/prompt named) perturbs only the first-λ-jet, which the retained pressure carrier does not
   keep (route-2 §4), so it does NOT bite `C_M`. The carrier bridge is nonetheless able-to-fail via a FORM change of
   the material `normal_exact` (both legs verified). ⇒ **the carrier bridge is a LIVE control** — a certified-zero
   carrier bridge is genuine reconciliation, ⛔ not a dead knife. (Directive note for reuse: the carrier-bite target is
   the material normal at retained grade, ⛔ not `material_inverse_transpose`.)
2. **Cosmetic (agent only):** the reconcile `emit` override stamps `target=(1,-1,0)` dimension METADATA on `SOURCE_*`;
   metadata only, ⛔ never touches the payload/residual. Non-blocking.

⇒ Instrument ACCEPTED (both legs CLEAR; no repair). NEXT = the orchestrator adjudicates the `R_N6` disposition from the
`.out` (does `R_N6` vanish modulo the defining relations — i.e. does the residual localize to the sanctioned N4
constitutive SOURCE channel with the geometric CARRIER bridge reconciled), ⛔ NOT pre-judged.
