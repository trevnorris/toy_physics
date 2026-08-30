# Build review — S11c-b WL admissibility repair (Codex-written; two legs, both SOUND)

Codex-written engine change → legs = **fresh Claude agent + Grok** (rule 7). Both FORM-ablated the `.wl` on
`/tmp` copies (working tree untouched), ran real Mathematica kernels under `timeout 600`, one at a time
(serialized for the 2-seat licence), and derived the normalization independently in SymPy. Both verdict:
**SOUND**, no defects.

The change (3 lines, `constructFullFieldBackgroundEnergy`): mixed invariant [[8]]
`Dot[gradTheta, gradLocalEw]` → `Dot[gradTheta, gradFullEw]`, `gradFullEw = anchoredWidth^(-1)
gradient[fullWidth]` (was `gradLocalEw = gradient[localEw]`); `git diff` = L531/L543/L555 only.

## 1. FORM ablation (mandatory) — PASS, the lift is load-bearing
Both legs reverted `gradFullEw` → `gradient[localEw]` (perturbation-only) on a `/tmp` copy and re-ran:
- Baseline θ body force (both anchorings): `gradientThetaEwCoefficient·σ_W·(−1+η·w1)·(∇²w1)/(L_W·W_0)` — nonzero.
- FORM-ablated θ body force: **`0`** (both anchorings).
The θ background force exists ONLY because of the full-field lift; the mixed `∇θ·∇gradFullEw` invariant is the
sole background-order θ contributor (every other θ-bearing invariant is a perturbation bilinear, order²→0 at
𝔅⁰). A coefficient rescale would not do this — the structural revert does. θ computed at `eulerScalar[
firstVariation, thetaField]` (L1339–1340), `firstVariation = D[fullFieldEnergy, backgroundOrder]/.0`
(L1334–1335); mixed-invariant site L555, `gradFullEw` L543.
- Grok: `/tmp/S11cb_form_baseline_harness.wl`→`.stdout`, `/tmp/S11cb_form_ablate_harness.wl`→`.stdout`.
- Claude agent: `/tmp/s11cb_review/baseline.wl`→`.out`, `/tmp/s11cb_review/ablate_pertonly.wl`→`.out`.
Both cross-checked that the engine's OWN emitted `ADMISSIBILITY_OPERATOR_OPERAND` THETA equals their probe
byte-for-byte (probe faithful).

## 2. Normalization `/W_bg` — PASS, mandated and matched (orchestrator-verified, rule 13)
Both legs derived independently (SymPy, no engine/run read) that `∇W/W_bg` is the §3a/§3d-correct
representative, not `/W_0` (they differ at retained order `σ_W·η`; the code's `(−1+η·w1)` factor is exactly
the `σ_W·η` piece `/W_0` lacks — Grok confirmed with a second live FORM ablation swapping `/W_bg`→`/W_0`,
which drops the `etaBg` factor).
- Decisive non-circular determinant (Claude agent, verified by orchestrator): the committed pure-thickness
  invariant [[7]] = `½·kappaW·anchoredWidth^4·anchoredWidth^(-2)·|∇fullWidth|²` = `½κ_W·W_bg²·|∇W|²`, whose
  representative is `|anchoredWidth^(-1)∇fullWidth|² = |∇W/W_bg|²`. So [[7]] already uses `∇W/W_bg`; the
  repaired [[8]] shares that primitive (MATCH under `/W_bg`; FAIL under `/W_0`).
- Uniform limit (η→0, ∇W_bg→0): both reduce to the S11b `∇θ·∇e_W` invariant (residual 0).
- ⚠ Honest caveat (both legs): §3d prose + the uniform limit alone do not UNIQUELY pin `/W_bg` vs `/W_0`; the
  determination rests on internal consistency with the committed [[7]] and the local-strain convention
  `localEw = e_W,bg = δW/W_bg` (established/reviewed in prior rounds, unchanged by this 3-line diff). Given
  [[7]] as committed, `/W_bg` is forced.
- Grok: `/tmp/S11cb_independent_mixed_lift_derive.py`→`.stdout`, `/tmp/S11cb_form_ablate_W0norm_harness.wl`.
  Claude agent: `/tmp/s11cb_review/normalization_derivation.py`→`.out`.

## 3. N15 untouched — PASS
`newInvariantExpressions` byte-identical to HEAD (both: md5 match). No spurion `∇e_W` lifted ⇒ no second
background jet injected (no double-count).

## 4. Scope byte-identity — PASS
`git diff` = 3 lines in `constructFullFieldBackgroundEnergy`; file line count unchanged (1896). Byte-identical
to HEAD (md5): both `kineticEw = muW WZero^2 …` (L838/L923), `constructEnergyData`,
`backgroundBalanceFromModel`, §3b slab operator+origins, §3c coupling kernel. The single sink of
`constructFullFieldBackgroundEnergy` is `backgroundBalanceFromModel`; its effect reaches only the
admissibility-operator family (operand + its REP_INVARIANCE/CONTROL_INDEPENDENCE controls + residual). §5c
`UNIFORM_LIMIT_RESIDUAL`, slab operator, coupling kernel, `MU_THETA` do not consume it — unaffected.

## 5. Blindness — PASS
No `Import`/`Get`/`Needs`/`ReadString`/`<<`/`.py`/`sympy` in the engine; imports nothing, re-derives from the
specs. The diff adds no sibling read.

## Deliverable
Regenerated transcript `mathematica/out/S11c_b_brane_operator_mathematica_audit.out` (97,055,766 B, 29 tags,
none missing, no kernel error, 201 s run, SHA-256 `976481da…31ea3`).

⇒ Build SOUND. Commit the reviewed engine + transcript. The θ body force now matches the sibling's structure
(`(η·w1−1)·∇²w1·σ_W/(L_W·W_0)`) — but agreement is CONFIRMED by the step-4 comparator re-run + re-adjudication
(rule 13), not asserted here.
