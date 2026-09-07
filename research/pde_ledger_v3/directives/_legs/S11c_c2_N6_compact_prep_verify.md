# Compact-prep verify — are the c2 N6-RESOLVED state docs accurate + in the clear?

You are Codex (`gpt-5.6-sol`, xhigh), verifying the session's state records before a compaction. ⛔ Document-only;
⛔ do NOT modify the tree; run no CAS. Read the sources + `git log`, confirm each claim, flag any over-claim/error
with `file:line`. Working dir `/var/projects/toy_physics`; repo paths under `research/pde_ledger_v3/` unless absolute.
Relevant commits (`git log --oneline -16`): `50bc5d92` STATUS N6-resolved, `d21c8ff5` N6 RESOLVED, `326e6123` covariance
BUILD, `123d9a18` covariance directive CLEAR, `4333ecb1` framing settled, `dd9f0fe3` reconcile over-clear CORRECTED,
`af343237` reconcile build clear + (over-cleared) adjudication, `36d59b95` reconcile BUILD, `08d72d46` reconcile
directive CLEAR.

## The claimed arc (verify each against the commits + records)
1. **Reconcile:** E1 question-vet (Codex-sol) corrected "R_N6 = J" → "R_N6 vanishes MODULO the defining relations".
   Reconcile directive CLEAR (3 rounds) `08d72d46`; astra `scripts/S11c_c2_N6_reconcile_sympy.py` 2 build legs CLEAR
   `36d59b95`. Established: geometric CARRIER reconciles (`C_E=C_M`) + residual localizes to the SOURCE channel.
   Records `_measurements/S11c_c2_N6_reconcile_{directive_review_adjudication,build_clearance,adjudication}.md`.
2. **The over-clear (verify this is recorded honestly, not buried):** my "N6 SATISFIED" adjudication `af343237` was an
   OVER-CLEAR; 2 adj-review legs SPLIT (Grok SOUND / Codex-sol OVER-CLEAR); I verified Codex RIGHT and CORRECTED
   `dd9f0fe3` — given `ΔC=0`, `R_N6=B(C_M,ΔS)` is an algebraic tautology (localization ≠ quotient reduction).
   `_measurements/S11c_c2_N6_reconcile_adjudication.md` should carry the correction prominently.
3. **The sufficient test:** framing vet → BOTH engines converged on `R_cov = ms − source_terms(μ_E.subs(Φ), V_E)`
   (source-naturality, non-circular, Φ prolonged through μ_E's rank-2 jets). Directive CLEAR (2 rounds) `123d9a18`;
   astra `scripts/S11c_c2_N6_covariance_sympy.py` 2 build legs CLEAR `326e6123`. Records
   `_measurements/S11c_c2_N6_{sufficient_test_vet_adjudication,covariance_directive_review_adjudication,covariance_build_clearance}.md`.
4. **The result + interpretation:** `R_cov=0` all 4 cases (δ≈2.6e-22), knives bite (84,4) ⇒ material builder
   implements Φ ⇒ per-engine N6 PASSES as OPERATOR COVARIANCE (Reading B, Codex-astra verdict, user-adopted `d21c8ff5`).
   `_measurements/S11c_c2_N6_RESOLVED.md`.

## Verify the STATE DOCS are accurate + not over-claimed
- **STATUS.md** new top clause (the "c2 N6 per-engine ✅ RESOLVED" clause) — accurate vs 1-4; the prior clause marked
  `[SUPERSEDED]`; no stale "NEXT = adjudicate R_N6 at reconcile".
- **Memory** `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md` (the R_N6
  RECONCILE bullet, ≈ frontmatter + body ~237-270) + `…/MEMORY.md` (the S11c-c pointer line) — accurate + consistent.
- The three CARRY-FORWARD premise caveats are present in `S11c_c2_N6_RESOLVED.md` + STATUS + memory: (1) is Φ itself
  physically correct; (2) does V transform (`V_E≡V_M` = builder-agreement, prediction uses `V_E` not `Φ(V_E)`); (3)
  extracted-block leakage. And the terminology-fix (`I_{M→E}`) + preserve-both-findings for the step record.

## Flag specifically
- Any **over-claim**: N6 stated as fully DONE/cross-engine-verified (it is PER-ENGINE SymPy covariance-passes; the
  **blind Wolfram N6 + cross-engine comparator + reconcile + step record are OWED**); `R_cov=0` stated as an absolute
  certificate (it is "no nonzero found" at conditional δ≈2.6e-22); the covariance verdict stated as strict `R_N6=0`;
  the 3 premise caveats dropped or mis-stated; the over-clear buried rather than recorded.
- Any **factual error** vs commits/records (SHAs, what cleared, what is owed).
- Any **owed item omitted**: blind Wolfram N6; comparator/reconcile (matched representations); step record; the
  ephemeral ~499 MB `/tmp` .out. (F/G re-grounding is PAUSED INDEFINITELY per user — should NOT be listed as owed.)
- **MEMORY.md compaction integrity:** the S11c-c pointer is accurate; ⛔ flag if any obviously-live pointer looks
  dropped or a hook now misstates its topic (I trimmed hooks 22→20 KB, kept all 151 links; 0 broken files verified).

## Output
End with **CLEAR TO COMPACT** (docs accurate + open items correctly recorded) or the exact fix list. Brief.
