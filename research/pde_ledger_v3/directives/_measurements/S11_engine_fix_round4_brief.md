# Measurements — S11_engine_fix_round4_brief.md

Generated mechanically by `gen_twin_fix4.sh` (session scratchpad): every block below is the
command followed by its literal output, captured at generation time. Nothing is transcribed.
Analyst probe suites: /tmp/s11d4_leg_0mS5/ (agent; summary with supersession note at
~/.s11_build/d4_analysis_agent.md) and /tmp/s11d4_leg_X8uC/ (Grok; report at
~/.s11_build/d4_analysis_grok.log). Review legs on this brief: ~/.s11_build/fix4_leg_codex.log,
~/.s11_build/fix4_leg_grok.log.

```
$ git rev-parse --short HEAD
60dc266f
```

## Item 1 — the zero test and its consumers

```
$ sed -n '944,948p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
def algebraic_zero_test(entry: sp.Expr) -> bool | None:
    if entry == 0 or entry.is_zero is True:
        return True
    if not getattr(entry, "free_symbols", set()):
        return False if entry != 0 else True
```

```
$ sed -n '1090,1096p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
        return [sp.ImmutableMatrix(cols, 1, lambda i, _j, col=col: sp.Integer(1) if i == col else sp.Integer(0)) for col in range(cols)]

    undecided_candidates: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for row_set in combinations(range(mat.rows), rank):
        for col_set in combinations(range(cols), rank):
            minor = exact_determinant(mat.extract(row_set, col_set))
            status = algebraic_zero_test(minor)
```

The empty-basis hazard's return site:

```
$ sed -n '1115,1117p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    detail = f"; last undecided candidate failed after {type(last_error).__name__}" if last_error is not None else ""
    record_issue(f"generic nullspace basis unavailable: no usable source minor found{detail}", root=root)
    return []
```

Analyst measurements (paths above): all candidate minors undecided at D3 (9/9,
out_d3_null_whole.txt) and both probed D4 minors (out_d4_minor3.txt); the rational root's scan
decides at its first minor. Multi-radical stream content bounding the class (preserved stream):

```
$ rg -c 'ROOT_COINCIDENCE.*K_SOLUTION' /home/trevnorris/.s11_build/sweep2_partial_through_D4R1.out
26
```

```
$ rg -o 'PY_S11_XKIN_ANISO_D3_ROOT_COINCIDENCE_R2_R3_K_SOLUTION' /home/trevnorris/.s11_build/sweep2_partial_through_D4R1.out | head -1
PY_S11_XKIN_ANISO_D3_ROOT_COINCIDENCE_R2_R3_K_SOLUTION
```

(The round-4 Codex leg parsed these tags: up to 4 distinct radical bases per scalar,
probe_radical_classes.py in its scratch.)

## Item 2 — the pivot inverse

```
$ sed -n '1060,1066p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    for free_col in free_cols:
        rhs = -sp.ImmutableMatrix([mat[row, free_col] for row in selected_rows])
        try:
            pivot_values = selected_minor.inv() * rhs
        except Exception as exc:
            record_issue(f"generic nullspace pivot inverse unavailable after {type(exc).__name__}; used gauss_jordan_solve")
            pivot_values = selected_minor.gauss_jordan_solve(rhs)[0]
```

Analyst timings: D3 2×2 radical pivot inv() 13.3 s (out_d3_null_whole.txt); D4 3×3 inv() did
not complete in 500 s (out_d4_inv3.txt).

## Item 3 — the determinant route (staged, orchestrator re-run of the Codex leg's probes)

```
$ sed -n '1002,1011p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
def exact_determinant(matrix: sp.MatrixBase) -> sp.Expr:
    simplified = reduced_matrix(matrix)
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
$ timeout 590 python3 -u /tmp/s11f4_leg_qNaO/probe_exact_det_stages.py 2>&1 | tail -7
STAGE reduced_matrix seconds=0.367
STAGE exact_domain_matrix_already_reduced seconds=0.012
DOMAIN EX
STAGE DomainMatrix.det.as_expr seconds=36.874
DOMAIN_DET ops=209247 srepr_len=5044921 sqrt_occurrences=663
STAGE reduced_expr_domain_det seconds=30.419
REDUCED_DOMAIN_DET ops=41342 srepr_len=990239 sqrt_occurrences=140
```

```
$ timeout 590 python3 -u /tmp/s11f4_leg_qNaO/probe_same_minor_reduction.py 2>&1 | tail -6
MINOR (0, 1, 2) RAW seconds=1.097 ops=3825 srepr_len=92062 sqrt_occurrences=12
MINOR (0, 1, 2) REDUCED seconds=1.053 ops=1206 srepr_len=29284 sqrt_occurrences=4
MINOR (0, 1, 2) VALUE_RESIDUAL_STATUS True
MINOR (1, 2, 3) RAW seconds=2.365 ops=10755 srepr_len=258336 sqrt_occurrences=36
MINOR (1, 2, 3) REDUCED seconds=2.815 ops=10755 srepr_len=258336 sqrt_occurrences=36
MINOR (1, 2, 3) VALUE_RESIDUAL_STATUS True
```

## Item 4 — the simplifier sites, both

```
$ sed -n '1392,1403p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    for idx, root in enumerate(roots, start=1):
        m_r = sp.ImmutableMatrix(route.matrix.subs(omegaSquared, root))
        rank = matrix_rank(m_r)
        nullity = n - rank
        stacked = sp.ImmutableMatrix.vstack(m_r, sp.ImmutableMatrix([list(k)]))
        stacked_rank = matrix_rank(stacked)
        transverse_nullity = n - stacked_rank
        m_dot_k = sp.simplify(m_r * k_vector)
        basis = generic_nullspace_vectors(m_r, rank, root=idx)
        dot_k = sp.Tuple(*[sp.simplify((vec.T * k_vector)[0]) for vec in basis])
        residuals = sp.Tuple(*[sp.simplify(k_sq * sp.ImmutableMatrix(vec) - ((vec.T * k_vector)[0]) * k_vector) for vec in basis])
        basis_payload = sp.Tuple(*[sp.ImmutableMatrix(vec) for vec in basis])
```

```
$ sed -n '1918,1923p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
        m_dot_k = sp.simplify(m_r * k_vector)
        basis = generic_nullspace_vectors(m_r, rank, root=idx)
        dot_k = sp.Tuple(*[sp.simplify((vec.T * k_vector)[0]) for vec in basis])
        residuals = sp.Tuple(*[
            sp.simplify(k_sq * sp.ImmutableMatrix(vec) - ((vec.T * k_vector)[0]) * k_vector)
            for vec in basis
```

The scoped site is unexercised in the preserved cells (no admitted strata), hence the driven
demonstration in acceptance 4:

```
$ rg '^PY_S11_XKIN_ANISO_D[23]_STRATUM_ORDERING:' /home/trevnorris/.s11_build/sweep2_partial_through_D4R1.out
PY_S11_XKIN_ANISO_D2_STRATUM_ORDERING: Tuple()
PY_S11_XKIN_ANISO_D3_STRATUM_ORDERING: Tuple()
```

## Spec cover

```
$ sed -n '404,406p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md

⭐⭐ **WHY N3 EXISTS.** `null(M_r) ∩ k^⊥` is exactly the null space of `M_r` stacked on `kᵀ`, so `nu_T`
counts the null directions having **no component along `k`** — and it is a **rank**, so it is
```

```
$ sed -n '409,410p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
⛔⛔ **CLASSIFYING THE VECTORS RETURNED BY `NullSpace` IS NOT A SUBSTITUTE AND MUST NOT BE THE PRIMARY
RESULT.** `NullSpace` returns an **arbitrary** basis. If a root's null space contains directions that are
```

```
$ sed -n '418,421p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
stiffness functional; nothing in this file says how any package's null directions lie relative to `k`.** ⛔
Any code path that assumes a two-way parallel/perpendicular split will report a wrong count without failing.

⇒ ⭐ **`N6` is retained for display only.** Emit the **residual objects**, ⛔ not a classification token as
```

The spec's only simplify-related lines (neither mandates nor forbids a Q4 simplifier; §Q5's
fully-simplified ratio object is out of this brief's scope):

```
$ rg -n 'simplify' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
426:routine's implied count — which is what lets it fail. ⛔ Do not "simplify" it to `rank + nu − D`: with `nu`
492:⭐ Build this inventory from the package's **DECLARED additive terms** as §7 gives them, before simplifying
```

## The run this round responds to

```
$ sed -n '276p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
        print(f"{tag}: {render(payload)}", flush=True)
```

```
$ ls -la /home/trevnorris/.s11_build/sweep2_partial_through_D4R1.out /home/trevnorris/.s11_build/partial_out_preserved_20260813.out
-rw-rw-r-- 1 trevnorris trevnorris  3209481 Aug 13 22:34 /home/trevnorris/.s11_build/partial_out_preserved_20260813.out
-rw-rw-r-- 1 trevnorris trevnorris 24414381 Aug 14 08:27 /home/trevnorris/.s11_build/sweep2_partial_through_D4R1.out
```

```
$ grep -o 'PY_S11_XKIN_ANISO_D4_[A-Z0-9_]*' /home/trevnorris/.s11_build/sweep2_partial_through_D4R1.out | tail -1
PY_S11_XKIN_ANISO_D4_ROOT1_N7_RESIDUAL
```

The sweep started 22:35:25 Aug 13, completed everything except XKIN_ANISO D4-D5 in ~46 min,
ground CPU-pegged (etime≈cputime, flat ~500 MB RSS) inside the D4 second-root chain from
00:53:22, and was killed at ~09:50 by the user's decision (ps observations in the session
record; the process is gone, so they are historical).

## Not measurements — decisions, marked as such

- Item 3 was REVERSED from this brief's pre-fold draft after the Codex review leg's staged
  measurement showed the swell is the DomainMatrix.det() stage, not the reduction; the
  orchestrator re-ran both probes (outputs above) before adopting. The agent analyst's
  contrary attribution is preserved with a supersession note.
- The pre-fix baseline stream is deliberately NOT referenced in the brief the builder
  receives; the post-build leg prompt carries it (round-1 lesson: a visible stream must not
  be a builder target).
- Acceptance shapes (real-minor operands; pinned D4 blocks; identity-class residual checks;
  driven scoped demonstration): orchestrator decisions incorporating both legs' defective-fix
  constructions.
- The D4 completion target is the user's (2026-08-14): "Let's set our target to D=4, and if
  we can just get that out then I'll be happy with it." D5 declared not necessary.
