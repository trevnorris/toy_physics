# Measurements behind `S11_engine_fix_round1_brief.md`

Commands and their literal output. CLAUDE.md rule 2. Regenerate; do not transcribe.

## 1. The isolated D5 reproduction and its verdict
```
$ ulimit -v 3000000; python3 ~/.s11_build/repro_d5.py
Traceback (most recent call last):
  File "/tmp/repro_d5.py", line 10, in <module>
    r = m.run_cell('MAIN', 5, q9)
  File "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py", line 1380, in run_cell
    q4_data = emit_q4(package, n, route, roots, assumptions)
  File "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py", line 1239, in emit_q4
    rank = matrix_rank(m_r)
  File "/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py", line 870, in matrix_rank
    return int(simplified.rank(iszerofunc=lambda entry: entry == 0, simplify=False))
  File "/home/trevnorris/.local/lib/python3.10/site-packages/sympy/matrices/matrixbase.py", line 3115, in rank
    return _rank(self, iszerofunc=iszerofunc, simplify=simplify)
  File "/home/trevnorris/.local/lib/python3.10/site-packages/sympy/matrices/reductions.py", line 242, in _rank
    _, pivots, _ = _row_reduce(mat, iszerofunc, simpfunc, normalize_last=True,
  File "/home/trevnorris/.local/lib/python3.10/site-packages/sympy/matrices/reductions.py", line 127, in _row_reduce
    mat, pivot_cols, swaps = _row_reduce_list(list(M), M.rows, M.cols, M.one,
  File "/home/trevnorris/.local/lib/python3.10/site-packages/sympy/matrices/reductions.py", line 109, in _row_reduce_list
    cross_cancel(pivot_val, row, val, piv_row)
MemoryError
EXC_TYPE MemoryError
EXC_REPR MemoryError()
TAGS_EMITTED 62
LAST_TAG PY_S11_MAIN_D5_ROOT_COINCIDENCE_R1_R2_COEFF_REAL_STATUS_OPERANDS
ISSUES []
```

## 2. The defective zero test, both call sites
```
$ sed -n '868,870p;920p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
def matrix_rank(matrix: sp.MatrixBase) -> int:
    simplified = sp.Matrix(matrix).applyfunc(lambda entry: sp.factor(sp.cancel(entry)))
    return int(simplified.rank(iszerofunc=lambda entry: entry == 0, simplify=False))
    basis = simplified.nullspace(iszerofunc=lambda entry: entry == 0, simplify=False)

$ grep -n 'matrix_rank(\|iszerofunc' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
868:def matrix_rank(matrix: sp.MatrixBase) -> int:
870:    return int(simplified.rank(iszerofunc=lambda entry: entry == 0, simplify=False))
920:    basis = simplified.nullspace(iszerofunc=lambda entry: entry == 0, simplify=False)
1239:        rank = matrix_rank(m_r)
1242:        stacked_rank = matrix_rank(stacked)
```

## 3. Diagnostics are only reachable after the whole loop
```
$ grep -n 'except Exception' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
756:    except Exception as exc:  # exact Groebner can be unavailable for some coefficient domains
775:        except Exception as exc:
955:        except Exception:
1317:    except Exception as exc:
1450:    except Exception:
1668:            except Exception as exc:

$ sed -n '1666,1670p;1678,1688p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
                if EMITTER.count == before_count:
                    ISSUES.append(f"{package} D{n}: cell completed without emitted tags")
            except Exception as exc:
                ISSUES.append(f"{package} D{n}: cell skipped after {type(exc).__name__}: {exc}")
    declared_pairs = [(package, n) for package in PACKAGE_ORDER for n in PACKAGE_DIMS[package]]
    main_completed = {(package, n) for package, n in completed_pairs if package == PRIMARY_PACKAGE}
    main_declared = {(PRIMARY_PACKAGE, n) for n in PACKAGE_DIMS[PRIMARY_PACKAGE]}
    if main_completed == main_declared:
        ledger = merged_export(main_dim_data, run_pairs_payload, skipped_pairs_payload)
        write_exports(ledger)
    else:
        stale = SCRIPT_DIR / "S11_exports.py"
        if stale.exists():
            stale.unlink()
        ISSUES.append("S11_exports.py not published because a declared MAIN cell did not complete")

```

## 4. MAIN D5 never completed — 188 of D4 249 suffixes absent
```
$ comm -23 <(D4 suffixes) <(D5 suffixes) | wc -l ; distinct counts ; last D5 tag
188
D4 distinct: 249
D5 distinct: 61
last D5 tag: PY_S11_MAIN_D5_ROOT_COINCIDENCE_R1_R2_COEFF_REAL_STATUS_OPERANDS
absent at D5, incl. the two that fill main_dim_data:
COEFFICIENT_ORDERING
DIM_COEFFICIENTS
```

## 5. The script has no out/ writing
```
$ grep -c "out/\|OUT_DIR" research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
0
committed out/ file: 825 lines (truncated casualty)
```
