---
unit_id: 011
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T17:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 011 (post-trim cycle)

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim)

**Classification:** resolved

**Resolution context:**
The original directive contained a `## Resolve before fix_loop` block routing
to the user, who chose **direction (c) "split across stages"**: keep only the
core K-eliminated compatibility surface in stage 011, and rely on the fact that
the destination stages (022, 023, 024) already publish the relevant signatures
so no paper move is required — just trim the stage 011 scripts.

**What changed:**

`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py`
was rewritten end-to-end (now 86 lines vs. ~211 prior). The retained content is
exactly the three core deliverables enumerated by `stage_011.tex`:

- `K_one_pole_p` and `K_norm_p` Solve-derived from the one-pole / fixed-target
  normalization equations (L40-47);
- `compat_direct_p = (N0+eps n0)/Ptarget - 3 (S+eps z2)^2/(T+eps z4)` and the
  three-line K-elimination identity (L54, L59);
- first variation `d_compat_direct - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2)`
  (L64) and `sp.diff(d_compat_direct, z0) == 0` (L65).

Removed (relative to the audited version):
- transported-target K-surface + compatibility variant (was sympy L132-155);
- constant-prefactor branch factorizations (P2=0, P4=0) (was L168-176);
- real-Y20 weak-axisymmetric lane lambdas (1, 1/2, -1) (was L179-185);
- grouped trace + b=3a (was L188-196);
- static Xi1 prefactor slope na/N0+za/D0 (was L199-204);
- per-lane u2 slope (was L207-211);
- dP0/dP2/dP4 closed-form readback prints (was the prelude to all of the above).

`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl`
was correspondingly trimmed to 64 lines. Retained: M1 fixed-target K-eliminated
surface vs closed form (L38-41), M2 one-pole K shift (L43-47), M3 normalization
K shift (L49-53), M4 fixed-target compatibility first variation (L55-59), M5
z0 independence via `D[compatFixedShift, z0]` (L60). Removed: transported-target
K-surface (was L97-116), Y20 lane lambdas + same-sign selection (was L118-133),
grouped trace + b=3a (was L141-142), Xi1 slope (was L146), u2 lane slope
(was L150-153), and the entire Gaunt/ThreeJSymbol scaffolding.

**Assessment:**

1. **Trim is complete and surgical.** All five clusters called out in the
   directive are gone from both engines. No "extra" code remains in either
   script that goes beyond what `stage_011.tex` enumerates.

2. **Three core deliverables remain in both engines.** Cross-check:
   - K-eliminated compatibility surface: sympy `assert_zero` at L59,
     Mathematica `check["M1 …", compatFixed - compatFixedClosed]` at L38-41.
   - First variation: sympy L64, Mathematica L55-59.
   - z0 cancellation: sympy `sp.diff(d_compat_direct, z0)` at L65,
     Mathematica `D[compatFixedShift, z0]` at L60.
   The sympy script also keeps the intermediate K-shift chain (`dK_one_pole`,
   `dK_norm`, `d_compat = dK_norm - dK_one_pole`) which constitutes the
   non-trivial algebraic content backing the z0-cancellation claim — i.e., the
   z0-drop is observed in the actual K-elimination chain at L60-63 rather than
   merely being asserted on a closed form built without z0. This is the right
   thing to keep: it is precisely the substrate of A4/A7/A8 that the original
   audit's self-test notes called out as load-bearing.

3. **Non-tautology.** None of the surviving assertions compare `lhs - lhs`.
   The compatibility-surface identity at sympy L59 derives `K_norm_p` and
   `K_one_pole_p` independently from `sp.solve(...)` and then asserts the
   difference equals the closed form — both sides come from independent paths.
   The z0-cancellation at L65 is checked by differentiating an expression
   whose construction at L54 nominally still contains z0-bearing summands
   (because `compat_direct_p` is built from `(N0+eps n0)/Ptarget - 3 (S+eps z2)^2/(T+eps z4)`
   which does not carry z0 explicitly, so this single line is structural; but
   the load-bearing algebra is at L62-63 where `d_compat` is derived from the
   K-elimination chain that *does* carry z0 and then shown equal to
   `d_compat_direct` at L63). Same chain on the Mathematica side via `kPole`/
   `kNorm`/`compatFixed`.

4. **Engine parity preserved.** Both engines now verify the same five
   identities. No surviving claim is checked by only one engine. The Mathematica
   side is no longer a strict transliteration: SymPy goes through explicit
   `dK_one_pole = (lin(K_one_pole_p) - K_one_pole)/eps` extraction whereas
   Mathematica uses `firstShift[expr] := Coefficient[Normal[Series[...]], eps, 1]`.
   Same target, distinct mechanism.

5. **Paper alignment.** Now aligned. The trimmed scripts verify exactly what
   `stage_011.tex` enumerates as `Output`: K-eliminated compatibility surface,
   its first variation, z0-cancellation. No script-only extras remain, so the
   paper card's enumeration is now a faithful index of what the audit covers.
   `paper_alignment: aligned`.

## Exec log assessment

**SymPy:** exit=0 per the orchestrator's report in this verification request
(the directive's `Codex apply` summary confirms "Both scripts exit 0"). No
post-trim `stage_011_sympy.log` was captured at
`/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/`; the
directory contains only `stage_011_mathematica.log` from the May 21 cycle.
**Cannot independently quote post-trim sympy log lines.** Verifier rule forbids
re-execution.

**Mathematica:** exit=0 per the orchestrator's report. The captured
`stage_011_mathematica.log` at
`/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_011_mathematica.log`
is from 2026-05-21T11:43 — i.e., it predates the May 25 trim and still lists
M6-M11 (transported-target, Y20 lanes, grouped trace, Xi1 slope, u2 lane
slope). **The captured log does not reflect the post-trim script** and would
have to be regenerated by the orchestrator to be quotable here.

**Output freshness:** the saved transcript files
`scripts/output/moving_throat_pde_stage011_*.txt` and
`mathematica/output/moving_throat_pde_stage011_*.txt` are dated 2026-05-21
(pre-trim); the scripts are 2026-05-25 (post-trim). Saved outputs are STALE
relative to the trimmed scripts. The 2026-05-21 sympy output still contains
the full "1) Base isotropic … 8) Reading the result back …" sections, and the
2026-05-21 Mathematica output still contains M6-M11 residual lines. This must
be refreshed before the next manifest sync; flagging as an orchestrator-side
freshness gap, not a verification blocker (the scripts themselves are correct).

**Diff:** `redteam/exec_logs/stage_011_diff.patch` is 0 bytes (stale from
May 21). The actual May 25 trim is not captured in any diff under
`redteam/exec_logs/`.

## Material-change assessment

`material_change`: **false**.

The trim removes verified content from stage 011 but no previously-published
value changes. Per the resolution rationale, every removed cluster lives
downstream in a stage that already publishes the relevant signature:

- **Spot-check confirmed: stage 024 publishes the Y20 lane signature.**
  `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py`
  at L310-312 (and again at L343-345) sets
  ```
  lam20 = sp.Integer(1)
  lam21 = sp.Rational(1, 2)
  lam22 = -sp.Integer(1)
  ```
  and at L329-331 prints "lambda_20 = 1 / lambda_21 = 1/2 / lambda_22 = -1".
  Identical signature, now anchored at the stage that consumes it. Per the
  resolution rationale ("destination stages 022, 023, 024 ALREADY publish the
  relevant signatures"), the move is purely an attribution adjustment.

- Constant-prefactor branch (P2=0, P4=0 factorization), the static Xi1 slope,
  the u2 lane slope, and the transported-target K variant are stated by the
  user/Codex to be published at stages 022/023/024. Per verifier rules I have
  not opened those stage cards (scripts-only scope), but stage 024's Y20 spot
  check is sufficient evidence that the resolution rationale is sound for at
  least one of the moved clusters.

No downstream unit can have depended on stage 011 producing these results
*here* (rather than at 022/023/024) because the audit pipeline does not import
across stages — each stage's script is self-contained. The verified surface
of stage 011 shrinks; nothing else changes.

The orchestrator should still mark units > 011 as `upstream_stale: true` as a
matter of protocol, but no narrow re-audit is warranted on account of this
trim.

## Side observations (non-blocking)

- **Saved transcript outputs are stale.** The May 21 `.txt` outputs still
  enumerate the removed sections (constant-prefactor branch, Y20 lanes,
  Xi1, u2). These will mislead any downstream reader who reads the saved
  output instead of the script. Orchestrator should regenerate.
- **Exec logs not refreshed by the orchestrator.** The captured
  `stage_011_mathematica.log` is the pre-trim run (still showing M6-M11);
  no post-trim sympy log was captured at all. I am relying on the
  verification request's assertion that "Both scripts exit 0".
- **Empty diff.** `stage_011_diff.patch` is 0 bytes despite a substantial
  May 25 rewrite; the orchestrator's diff-capture step appears not to have
  run for this trim cycle.
- **Note for tracker sync.** Per the post-batch tracker convention, the
  STAGE_VERIFICATION_COVERAGE entry for 011 should be updated to reflect
  the shrunk surface (5 fewer claim clusters; M1-M5 only on the Mathematica
  side; sympy A4/A7/A8/A9/A13 chain only). The other 5 trackers
  (MATHEMATICA_MIRROR_POLICY, CHECKPOINT_TRUST_AUDIT,
  CHECKPOINT_CONSTANT_PROVENANCE, PAPER_CLEANUP_TRACKER,
  EM_PROJECTED_INTEGRATION_TRACKER) likely also need a note that the Y20 /
  Xi1 / u2 / constant-prefactor / transported-target content is now sourced
  at 022/023/024.

## Verdict justification

The single F1 paper_misalignment finding was routed to user resolution
(direction c, "split across stages"). Codex applied the corresponding trim to
both engines, removing the five script-only clusters that lacked paper anchors
at stage 011 while preserving the three core deliverables (K-eliminated
compatibility surface, its first variation, z0-cancellation) that the
`stage_011.tex` card enumerates as Output. Both surviving engines verify the
same content via non-transliterated paths (`(expr_lin - K_*)/eps` on the
SymPy side vs `Coefficient[Series[...], eps, 1]` on the Mathematica side).
Non-tautology is preserved by the K-elimination chain at sympy L60-63 and the
parallel `kPole`/`kNorm`/`compatFixed` chain at Mathematica L25-41.

The destination of the removed clusters is confirmed by spot-check at stage
024 (lane lambdas 1, 1/2, -1 present at sympy L310-312, L343-345). The trim
does not change any previously-published value: it relocates attribution.
`material_change: false`.

The verification is necessarily contingent on the orchestrator's assertion
that both trimmed scripts exit 0, because (i) the captured Mathematica exec
log predates the trim, (ii) no post-trim sympy log was captured, and (iii) the
saved `.txt` outputs and `.patch` diff are stale. Verifier rules forbid
re-execution; these are flagged as side observations for the orchestrator to
refresh, not as verification blockers.

stage 011: verified
