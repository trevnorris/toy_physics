# P2b §3a coefficient bridge — decision-leg review record (2 legs, computation-backed, convergent → v2)

Directive `directives/S11c_b_p2b_gamma_bridge_directive.md`. Legs Codex + Grok (orchestrator-written, rule 7). Logs
`~/.s11_build/S11c_b_p2b/decision_{codex,grok}.log`. ⚠ Grok's 1st + 2nd attempts died on an xAI capacity 500; the
leg-of-record is the recovered 3rd run. **v1 REJECTED by both (6/5 defects, file:line + literal probe).** Both
CONFIRMED B1's scale direction (`γ_WL = W_0·γ_PY`); both independently converged on the same holes. v2 folds once.

## Convergent defects folded into v2
1. **B1 — the string bridge cannot carry a scale** (both): apply path emits only a `Symbol`
   (`_reconcile_basic:279/310`, `_extra_basic:396`); `"W_0*gamma_…"` becomes a NEW ATOM (Grok probe), a `Mul` throws
   `TypeError` in `sp.Symbol(...)` (Codex probe). ⇒ B1 = an EXPRESSION-valued substitution `Mul(scale, PY_symbol)`
   via the Bridge-A / `s11ca.field_symbol` mechanism; the scale EMITTED, the fold made UNREPRESENTABLE.
2. **B3 — invariant-only certificate passes a folded scale / string rename** (both): `I_PY − W_0·I_WL = 0` holds
   even with NO coefficient scale (Grok `TERM_RES_noscale False`), and a fold with printed factor `1` is a false
   zero. ⇒ B3 must (a) LOCK the factor ∈ `{W_0,μ_R}` by source, computed as that symbol (not `I_PY/I_WL`, not `1`);
   (b) print the ENERGY-TERM residual `γ_PY·I_PY − subst(γ_WL)·I_WL` under the SAME substitution the residual applies
   (not a sidecar); (c) walk `ENERGY_BASIS_NEW_INVARIANTS` (currently `STRUCTURE_INCOMPLETE`), bijection BOTH
   directions (15/15/src, 30/30/branch, 60/60). Codex's decisive control: `RIGHT_PAIR_RESIDUAL_ZERO True /
   SWAPPED_PAIR_RESIDUAL_ZERO False` on the invariant relation.
3. **B2 — WL does NOT emit in quotient order; positional currently coincides with the pairing** (both): so a
   positional zip passes today (`POSITIONAL_EQ_SCALED True`) but is not the pairing. ⇒ B2 = the pairing METHOD (match
   unfolded `I` up to the family factor) + a PERMUTATION-INVARIANCE control + duplicate-record rejection; ⛔ the
   directive must not claim a "quotient order" object exists.
4. **B4 — the collapse is at the FULL parse + multiple sites + must not regress S11c-a** (both): the first collapse
   is `canonical_value`→`_extra_basic` (Grok: `widthBackground ∉ s11ca FIELD` so canon keeps args, `_extra_basic`
   drops them; `full_parse W_base-shift is_zero=True`); a canon-only fixture passes without fixing `_extra_basic`.
   ⇒ B4 fixture on the FULL parse for `widthBackground` AND `modulusBackground` AND `eWBackground`; fix the S11c-b
   tables, ⛔ NOT s11ca canon/`FIELD` (strips `thetaWave` legitimately → S11c-a regression); classify every
   `BARE_APPLIED` head (inert / live-background / trial-test-support-bare); ⛔ keep the two-site jet decode working.
5. **B5 — retire PROTECTED with a planted-residual test** (both): `07/10` + dead `gammaWidth*` divert a nonzero
   residual to `PROTECTED_UNREDUCED` (`:953`) and exclude ablation-touch (`:1478`) ⇒ B3 exercises `07/10` with a
   planted residual that must NOT route to `PROTECTED_UNREDUCED`. Scoped removal is safe (S11c-b-local, reverse-degree
   `{1:30}`).

## Status
v2 folded once (rule 7). ⇒ Codex build (after — or with — P1-WL's fresh WL `.out`; sequenced BEFORE P2a's final
validation) → 2 build legs (fresh Claude + Grok). The one genuine implementation change: the bridge gains an
expression-valued scale substitution + a factor-locked, energy-term, bijection certificate on the residual apply path.
