# Grok code review — GENUINE-REBUILD of the native-P Dirac constraint gate (READING-ONLY, no scripts)

⚠️ Do NOT write or run any scripts (script execution fails in this environment and kills your run). Review by READING the code and reasoning.

A PRIOR build of this Dirac constraint-analysis gate was caught FAKING the decisive result: hardcoded `additional_G_exists=False`, no Hamiltonian carrying the `{g_a}` couplings (they were only display strings), a faked Maxwell control (`cmat=sp.zeros(2)`), and a tautological `solve([t=-1,t=+1])=∅` as its "search." The REBUILD claims to fix all of that and reports `NATIVE_P_NO_EMERGENT_GAUSS`. Your job: by reading the code, judge whether the rebuild is a GENUINE computation or still fakes/shortcuts the decisive step. A negative verdict matching the builders' prior is the seductive-trap case.

Read: software/em_charge_attribute/native_p_gate_sympy.py (main target), native_p_gate_dual.wl, native_p_gate_compare.py, reports/native_p_constraint_gate.md (REBUILD NOTE claims Q1-Q6 closed), directive_native_p_gate_rebuild.md.

Judge, with file:line evidence:
1. Does build_H2 genuinely differentiate the Lagrangian -> Hessian -> Legendre -> constraints with {g_a} as free symbols that reach the computed PB matrix? Or hardcoded?
2. Does coupling_guard genuinely fail if a coupling drops?
3. Do all six controls run the IDENTICAL build_H2->dirac_search->search_G path as native-P (Maxwell's zero-PB an OUTPUT, not a literal)?
4. Is search_G a real enumeration that CAN return nonzero (GUARD-SEARCH-CAPABLE real)?
5. The tuned-strata FC=2 directions dismissed as "shift-symmetries not Gauss chains" — genuinely computed or a dodge to bury an EM candidate? Scrutinize hardest.
6. Is the Mathematica engine genuinely independent (own differentiation/Dirac/kernel), or reading SymPy's answer?
7. Any residual hardcoded verdict / tautology / X-equiv-X ablation / dropped-term or sign bug that would spuriously zero the first-class count?

Return VERDICT: VERIFICATION_CONFIRMS or VERIFICATION_ISSUES, then per-item findings with file:line, and a bottom line. Concrete over polite.
