---
batch_id: III.3
label: Part III.3 — Microclosure, gain thresholds, equilibrium, walls
started: 2026-05-22
completed: 2026-05-22
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: true
---

# Batch III.3 — Final summary

All 12 stages reached `verified` status. Both engines (SymPy + Mathematica) exit 0 on every unit; all 27 findings closed substantively across the 10 dirty stages; the remaining 2 (061, 066) passed the auditor's transliteration / tautology / hardcoded probes on first read. One stage (068) flagged `material_change: true` — the resonance-corrected `W_res`/`P_res` relations are now derived via `Solve` from the matched-threshold premises rather than postulated, but the symbolic content of the derived expressions is identical to the previously postulated forms; no downstream cascade.

## Per-stage outcomes

| Stage | Audit findings | Codex iterations | Verifier verdict | Notes |
|---|---|---|---|---|
| 061 | 0 (clean) | 0 | verified | Both engines pass 15 corresponding assertions with residual 0; soft-support / `eta -> oo` / stiff-support limit checks substantively exercise the claimed regimes |
| 062 | 3 (taut×2, math_translit) | 1 | verified | Hand-built `sigma_star` and EOS-monomial assertions replaced with a real parent-action Gaussian-elimination `sp.solve` for `sigma_star(rho)` plus an `n=5 EOS identity` wrong-exponent probe that prints a nonzero residual `K*rho^3*(5 - 6*rho)`; Mathematica side now derives `sigma_star` via `Series`/`SeriesCoefficient` rather than mirroring the SymPy choreography |
| 063 | 3 (taut, math_translit) | 1 | verified | All 8 SymPy + 8 Mathematica tautologies replaced with `sp.solve` / `Reduce` derivations and a Cauchy-saturation check that would fail under plausible bugs (e.g. `N_ss`-vs-`N_pp` swap in `G_max`); Mathematica side switched to `Reduce[... && gphiSq>0, ..., Reals]` to break the SymPy mirror |
| 064 | 4 (taut×3, math_translit) | 1 | verified | `chi_phi(y)` / `H(y)` profile actually instantiated; `I1`/`I2` integral overlap reductions now derived rather than carried as `const_subs` literals; F2 included one orchestrator-acceptable deviation (extra `Npp -> I1*Hw` substitution required for F1's independently-verified integral identity to canonicalize); Mma routed independently |
| 065 | 3 (taut, math_translit) | 1 | verified | Concrete-Gaussian shell integrals `J2=0` (by parity), `I1` polynomial coefficient expansion, and `gphi`'s `1/ell` scaling from `V_conf` differentiation now anchor docstring claims (1)-(3) and (6); `K_X`-cancellation identities preserved as remaining substantive checks |
| 067 | 4 (taut×2, hardcoded_translit, insuff_verif) | 1 | verified | Transverse norms `2 w_f` and `w_g sqrt(pi/2)` now integrated explicitly in SymPy via `.rewrite(sp.cosh)`-workaround `Integrate[sech^2]` and `Integrate[Gaussian^2]`; two "duality implies stationarity" blocks made non-tautological; numerical duality/monotonicity checks preserved |
| 068 | 3 (taut×2, math_translit) | 1 | verified — **material_change: true** | `Wfail_res := Pe_req/(C2*Delta_inf)` and `Wfail_match` cross-relation now derived via `Solve` from explicit resonance-corrected premises rather than postulated; `M3/M4` cross-relations `Wfail_res * C2 - Wfail_match` and `Wsuff_res * C2 - Wsuff_match` are the genuine load-bearing assertions; derived expressions match the prior postulated forms symbolically |
| 069 | 3 (taut, math_translit, prov) | 1 | verified (CHECKPOINT) | `Cres2`/`Wfail_res`/`delta_fail` definitional identities replaced by a parameterized `W_match` generator + monotonicity check (SymPy) and a `Cres2Prim` primitive + `Pres = 1/Cres2` derivation + `PresGap` via `Solve` (Mathematica); upstream `Pres`/`Wfail_match` carry-forward annotated with provenance comment (documentation fallback the directive explicitly permitted) |
| 070 | 1 (math_translit) | 1 | verified | Mathematica `.wl` inlined `1/Hw -> rhoW/(m*cSw^2)` and removed mirrored intermediates `J_1`, `gphi`, `I_1` while retaining the 3 `expectZero` calls; the substantive cancellations remain (assembled forms built from primitives; closed forms typed directly) |
| 071 | 1 (taut) | 1 | verified | `eta - L/ell` tautology replaced with a pin `K_m = pi a^2 hbar^2/(3 m rho_w)` and an `eta`-reconstruction from the closed-form `K_m` in both engines |
| 072 | 2 (insuff_verif, math_translit) | 1 | verified | F1 added 4 ratio-limit checks per engine comparing the full `Delta_0`/`Delta_inf` closed forms to the shell- and compression-dominated leading-order forms; SymPy collapses both shell ratios to `1` directly, Mathematica's `DeltaInf` shell ratio surfaces as the algebraically-equivalent surd `2/Sqrt[5] + 1/(5 + 2 Sqrt[5])` — divergent presentations confirm the cross-engine independence. F2 closed by orchestrator as "won't-fix-here, mitigated by F1" (the closed forms are derived in upstream stages; this unit only packages them) |

## Findings breakdown

- `tautological_check`: 14
- `mathematica_transliteration`: 9
- `insufficient_verification`: 3
- `hardcoded_result`: 1

Transliteration share back up to 9/12 vs III.2's 6/12. Stage cluster 062-068 (gain-thresholds / resonance / wall-shell) was particularly mirror-heavy; the auditor's signature III.3 pattern was "both engines retype the same hand-built leading-order form" for the asymptotic checks. Codex's substantive rewrites uniformly replaced the mirror with independent symbolic-limit / Series / Solve routes.

## Toolchain hardening landed mid-batch

None. The `ConditionalExpression` and `Limit` idioms documented after III.2 in `codex.md` carried III.3 without surfacing. The only fix-loop halt (072) was a directive-preordained Blocked, not a Mathematica idiom failure — the directive itself instructed codex to mark F2 Blocked with a specific reviewer question.

## Process observations

1. **Second batch with `material_change: true`** (stage 068, after stage 060 in III.2). Both flags trace to the same pattern: codex replaced a postulated relation with a `Solve`-derived one whose symbolic output matches the prior form. Downstream stages III.4+ are still `pending`; no `upstream_stale` cascade triggers. Worth tracking whether this becomes a recurring III.X pattern.

2. **Zero codex iter-2 fixes (matches II.1, III.1, III.2).** All 10 dirty stages cleared codex's iteration loop on first pass. The single halt (072 Blocked) was a directive-preordained block, not a codex failure; the orchestrator closed F2 as "won't-fix-here, mitigated by F1" after confirming both engines' new ratio-limit checks PASS exit 0 via independent symbolic-limit machinery.

3. **`codex.md` pitfalls section paid off.** The `ConditionalExpression`/`Limit` idioms patched into `codex.md` after III.2 (commit `2718222`) carried III.3 without surfacing. Codex's directive-mandated `Solve` calls in stages 062, 063, 064, 067, 068 either avoided the trap or codex itself stripped the wrapper as part of its iteration.

4. **One directive-preordained Blocked** (072 F2). The auditor's directive itself instructed codex to mark F2 Blocked with a specific reviewer question, since the `delta0`/`deltaInf` closed forms are derived in upstream stages and an independent re-derivation route can't be prescribed without that upstream context. Codex correctly followed the directive's instruction. Orchestrator resolved by accepting option (a) — F1's per-engine native-limit ratio checks provide the cross-engine independence guarantee F2 sought. Pattern worth noting: when the auditor knows in advance that a finding will be unblockable in-unit, the directive can preordain the Blocked block and let the orchestrator decide.

5. **One controlled codex deviation** (stage 064 F2). Codex added an extra `Npp -> I1*Hw` substitution beyond the directive's explicit edit; verifier confirmed this is necessary and grounded in F1's independently-verified integral identity, not a re-introduced tautology. Acceptable deviation.

6. **Exec_logs absent across all 10 verifications.** Verifiers noted that the orchestrator's `redteam/exec_logs/stage_NNN_{sympy,mathematica}.log` files were not captured by `fix_loop.sh` (only the diff was captured). All verifiers fell back to the freshly-regenerated canonical output `.txt` transcripts (mtimes post-edit) and confirmed exit 0 + PASS lines via that route. Non-blocking — the saved outputs are the substantive evidence — but worth tracking whether the next batch wants `fix_loop.sh` to also copy the sanity-exec output into `redteam/exec_logs/` for the verifier's audit trail.

## Prompt hardening landed this batch

None. The III.2 `codex.md` additions held.

## Tracker updates landed in this commit

- `notes/MATHEMATICA_MIRROR_POLICY.md`: 10 new Independent-Mirror Set entries (062-065, 067-072); snapshot bumped to "2026-05-22 (batch III.3 close)".
- `notes/CHECKPOINT_TRUST_AUDIT.md`: snapshot date updated with note on stage 069 (the one checkpoint in range 061-072); tier unchanged because the findings only affected verification non-tautologicity, not substantive claims.
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: snapshot date bumped; no new constants surfaced.
- `notes/PAPER_CLEANUP_TRACKER.md`: new P4-36 (III.3 batch) row and change-log entry dated 2026-05-22.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: III.3 is out of the linear projected-EM core range (004-021); only a date bump + completed-checks bullet noting the `material_change: true` flag on 068.
- `notes/STAGE_VERIFICATION_COVERAGE.md`: snapshot bumped; cumulative count updated to 72 / 253; new III.3 row in the per-batch coverage table.
