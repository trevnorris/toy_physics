# Claude+Codex consult — red-team remediation batch 6 (stages 143, 144, 146, 147)

Read-only codex-chat consult (ephemeral) on the four rewritten directives, run BEFORE
the fix loop. Prompt: `redteam/tmp_prompts/consult_batch6.md` (raw transcript was bloated
with file dumps and deleted per the consult-bloat caveat). Verdict scale:
CONCUR / DISPUTE / CONCEPTUAL-ESCALATE.

**SUMMARY (Codex): DISPUTE on Q3's claimed `gminus` sensitivity and Q5's claimed ability
to catch missing centering terms; all others CONCUR. No CONCEPTUAL-ESCALATE.**

## Q1 — 143 exp-remainder global positivity — CONCUR
The Taylor block (`directives/stage_143.md:72-82`) proves global positivity by repeated
integration from R(0)=R'(0)=R''(0)=0 and `R'''=exp(Pi)>0`; `Pi` is `positive=True` at
`...stage143...py:31`. Fails for the counterexample `Pi^3/6-Pi^4` (third-deriv check vs
`exp(Pi)` fails), unlike the old coefficient-only check at py:52-54. Mathematica `Reduce`
feasible — the script already runs a harder `Reduce[num>0,piM,Reals]` at wl:60-65.

## Q2 — 144 independent cleared-denominator Pi_* — CONCUR
`gThresholdResidual` (`directives/stage_144.md:101-103`) is exactly numerator minus
`gMinus·`(positive denominator from `...stage144...wl:35`); sign MUST be `(Exp[p]-1)` as
stated at directive:94-99. Genuinely different Mathematica route from the mirrored
`FindRoot[gPi==gMinus,{piM,1.5}]` at wl:46; anchoring to stage 131's owning value
(`scripts/output/...stage131...txt:2`) is the right non-self-anchor. No conceptual concern.
→ This settles the `needs_user_resolution` MIRROR-policy call: adopt the independent route.

## Q3 — 146 affine-law tolerance/precision fix — DISPUTE (justification overclaim; fix is sound)
The collapse `g_eps residual = (1-eps)(gPi(Pi_*)-gMinus)` is CORRECT (Sigma_eps + physical
moments built at `...stage146...py:101-105`; residuals at py:109-110 reduce to root +
quadrature residuals). Raising `nsolve` precision (`directives/stage_146.md:121-131`) is
mathematically sound and **closes the tolerance finding** (`codex_reviews/stage_146.md:26-34`).
**DISPUTE:** the check is NOT an independent guard against a wrong `gMinus`, because that
same `gMinus` is used to compute `Pi_star` at `...stage146...py:70-73` — in that coupled
failure mode the residual stays small. (It DOES still catch a wrong `c`/`Kq`/`Sigma_*`,
which do not feed the root.)
**Resolution (orchestrator):** keep the code edits as drafted (they close R1/R2). Correct
the directive's anti-tautology text to drop the overclaimed `gMinus` from the
"a bug in … perturbs the residual" list and record the coupled-mode caveat so the verifier
knows the affine check tests intercept-vs-direct-integral and the kernels, NOT `gMinus`
independently. Non-conceptual.

## Q4 — 147 source-moment definitions from the appendix — CONCUR
`paper/appendices/stage_appendix_part04.tex:798-824` is the operative source for `Sigma_*`,
normalized deformations, `gbar`, `Sbar`, `c`, `Kq`; `stage_147.tex:19` points back to it.
Retuning identity is exactly `deltaPi=-eps(gbar_sigma-g_*)/g'_*` (appendix:825-831).
Script's `Sigma_star=Pi_*/(1-S_*/4)`, `T_star=sqrt(9 Sigma_star/20)` (py:30-31) match the
branch law / `R_q(Pi_*)=1/4` at appendix:763-775.

## Q5 — 147 rigidity-kernel projection identity — DISPUTE (add a centering assertion)
The projection identity is analytically valid (appendix:833-850) and the implementation
(`directives/stage_147.md:151-174`) is a real can-fail check for coefficient sign / c↔Kq
swaps. **DISPUTE:** it is BLIND to an additive constant — omitting the centering terms in
`W_*` gives the same projection against `(sigma-Sigma_*)` because both profiles are
normalized (appendix:802-810). The added `D[...,x]==0` check (directive:220-225) only proves
the offset is CONSTANT in x, not that it equals `-A_T g_* - B_T S_*`. **A source-centering
assertion `∫ Sigma_* W_* dx == 0` is still needed.**
**Resolution (orchestrator):** ADD `∫_0^1 Sigma_*(x) W_*(x) dx == 0` to R2 (SymPy) and R3
(Mathematica). Verified sensitive: `∫Sigma_* W_* = A_T(gbar_{Σ*}-g_*) + B_T(S̄_{Σ*}-S_*) = 0`
since `gbar_{Σ*}=g_*`, `S̄_{Σ*}=S_*`, `∫Sigma_*=1`; dropping the `-A_T g_*-B_T S_*` constants
makes it `A_T g_* + B_T S_* ≠ 0`. This closes the centering gap the projection identity
cannot see. Non-conceptual (it strengthens the centering test the stage already claims).

## Q6 — 147 A_T autodiff + ratio-anchor honesty — CONCUR (with a label fix)
Autodiff of `T_m(Pi)=sqrt((9/20)Pi/(1-Sformula/4))` (`directives/stage_147.md:87-104`)
independently regenerates the chain-rule factors hand-expanded at `...stage147...py:33-38`,
so sign/power errors are caught. Keeping the `31.6785` ratio as a redundant computed
cross-check is fine — **but** the live script currently LABELS it as paper-quoted at
`...stage147...py:62-64`, while the paper only quotes `A_T`,`B_T` (appendix:843-848).
**Resolution (orchestrator):** add to the 147 directive an instruction to CORRECT the
misleading "paper-quoted" label on the `31.6785` ratio (script comment/print-label only;
describe it as the script's own computed ratio cross-check). Keep the check.

---

## Net actions before the fix loop
- 143, 144: apply as drafted (CONCUR). 144 `needs_user_resolution` → clear (MIRROR-policy
  settled: independent route).
- 146: keep code edits; rework the anti-tautology TEXT (drop `gMinus` overclaim + coupled-mode
  caveat for the verifier).
- 147: ADD the `∫ Sigma_* W_* = 0` source-centering assertion to R2/R3; ADD an instruction to
  correct the `31.6785` paper-quoted mislabel.
All Codex-agreed, none conceptual → resolved Claude+Codex per the standing rule.
