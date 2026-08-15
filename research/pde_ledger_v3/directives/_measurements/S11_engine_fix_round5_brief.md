# Measurements — S11_engine_fix_round5_brief.md

Generated mechanically (gen_twin_fix5b.sh, session scratchpad): every block is the command
followed by its literal output at generation time. Folded once after two review legs
(~/.s11_build/fix5_leg_codex.log, fix5_leg_grok.log); their pre-fold findings — a zero-return
defective fix passing the old acceptance, the lift route's radical-denominator counterexample,
disclosed counts/values, an overstated byte-exact claim — are in those logs. Analyst evidence:
/tmp/s11kw_leg_53JD/ (probeA-J + evidence_*.out; reconstruction verified by ALGEBRAIC equality,
stream_minus_reconstructed = 0) and /tmp/s11kw_leg_A5kv/ (00_DIAGNOSIS.md; D2/D3 rank-minor
payloads verified by STRING equality against the stream).

```
$ git rev-parse --short HEAD
ae105530
```

## The defect site

```
$ sed -n '1084,1096p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
def exact_determinant(matrix: sp.MatrixBase) -> sp.Expr:
    simplified = reduced_matrix(matrix)
    if any(quadratic_radical_nodes(entry) for entry in simplified):
        det_expr = simplified.det(method="bareiss")
    else:
        try:
            det_expr = exact_domain_matrix(simplified, already_reduced=True).det().as_expr()
        except MemoryError:
            raise
        except Exception as exc:
            record_issue(f"determinant exact-domain route unavailable after {type(exc).__name__}; used SymPy determinant fallback")
            det_expr = simplified.det()
    return reduced_expr(det_expr)
```

```
$ rg -n 'all_minors\(' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py | head -5
1128:def all_minors(matrix: sp.MatrixBase, size: int) -> tuple[sp.Expr, ...]:
2023:        rank_change = normalized_change_residuals(all_minors(m_r, rank)) if rank else tuple()
2024:        stacked_change = normalized_change_residuals(all_minors(stacked, stacked_rank)) if stacked_rank else tuple()
2156:        rank_minors = all_minors(m_r, rank)
2157:        stacked_minors = all_minors(stacked, stacked_rank)
```

## The cliff and the live run (historical ps values are in the leg logs)

```
$ stat -c 'out mtime: %y' research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out
out mtime: 2026-08-14 13:51:28.422861080 -0600
```

```
$ grep -o 'PY_S11_XKIN_ANISO_D4_[A-Z0-9_]*' research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out | tail -1
PY_S11_XKIN_ANISO_D4_ROOT1_STACKED_DROP_JOINT_REAL_STATUS_OPERANDS
```

```
$ grep -o 'PY_S11_XKIN_ANISO_D3_ROOT[12]_[A-Z0-9_]*' research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out | grep -A1 'ROOT1_STACKED_DROP_JOINT_REAL_STATUS_OPERANDS' | head -2
PY_S11_XKIN_ANISO_D3_ROOT1_STACKED_DROP_JOINT_REAL_STATUS_OPERANDS
PY_S11_XKIN_ANISO_D3_ROOT2_RANK_DROP_MINORS
```

Analyst probe outputs the brief's cliff claim rests on (sizes and budgets, not values):

```
$ grep -E 'TIMEOUT' /tmp/s11kw_leg_53JD/probeA_D4.out | head -2
D4R2_STACKED minor (0, 1, 2, 3)x(0, 1, 2, 3): P2_bareiss_det TIMEOUT_gt_300s (radical=True max_red_entry_srepr=7346 p1=0.64s)
```

```
$ grep -E 'TIMEOUT|OK D3' /tmp/s11kw_leg_A5kv/probe2_all_minors_growth.out | head -6
OK D3R2_one_leading_minor_det_size2 dt=0.453s rss=79.5MB 
OK D3R2_all_minors_rank dt=3.386s rss=79.5MB n_items=9 srepr_len_sum=64753 srepr_len_max=27133 zero_count=0
OK D3R2_all_minors_stacked dt=3.501s rss=79.5MB n_items=4 srepr_len_sum=62932 srepr_len_max=35079 zero_count=1
TIMEOUT D4R2_all_minors_stacked dt>540.000s rss=130.1MB
```

The review legs' oracle feasibility (their probes; the acceptance's two checks are runnable):

```
$ grep -E 'D4_SPECIALIZED|residual_zero' /tmp/s11f5_leg_JCGQ/probe_review.out 2>/dev/null | head -4 || tail -3 /tmp/s11f5_leg_JCGQ/probe_review.py
```

```
$ grep -E 'residual-zero|Laplace|16/16|5/5' /home/trevnorris/.s11_build/fix5_leg_grok.log | head -4
- D4 rank 3×3: prior completes; lift vs prior **16/16 residual-zero**, **16/16 srepr-identical**.
- D4 stacked 4×4: prior never used; lift vs Laplace-on-last-row (cofactors via prior 3×3 dets) **5/5 residual-zero** (so correct values *exist* and are checkable in ~2–3 min).
**Fix direction (brief only):** require a residual on the four non-zero size-4 blocks against an independent exact oracle that finishes under budget (e.g. Laplace expansion with size-3 prior dets, or a second exact route). Do not hang residual solely on “prior bareiss where it finishes.”
On *that* subclass, lift vs prior / Laplace is residual-zero (F1 probes). Multi-radical guard raises; nested synthetic hits two bases and raises.
```

The radical-denominator counterexample (review leg, re-runnable):

```
$ timeout 120 python3 /tmp/s11f5_leg_JCGQ/probe_rational_denominator.py 2>&1 | tail -4
denominator_contains_t=True
measured_probe_residual_zero=False
corrected_full_substitution_residual_zero=True
PROBE_RATIONAL_DENOMINATOR_DONE
```

## Not measurements — decisions, marked as such

- One item (the determinant bound). The radical locus/KW solve class after it is measured,
  deferred with reasons (D3 versions completed in a long-lived process; several D4 members
  bounded only from below at 50-74 s / >60 s) — a green round-5 run claims the determinant
  wall only. Orchestrator decision, judged defensible by both legs.
- Acceptance operands discovered dynamically and values verified by independent oracles:
  folded from both legs' convergent defective-fix construction (zero-return at the cliff size
  passing the old acceptance). No counts, roots, or values are stated in the brief.
- The class boundary (entries polynomial in one radical after reduction) narrowed from the
  pre-fold draft after a leg's denominator counterexample; boundary controls added.
- The user's standing target (2026-08-14): D=4 out; D5 not necessary; kill-fix-rerun standing
  order once round 5 is ready.
