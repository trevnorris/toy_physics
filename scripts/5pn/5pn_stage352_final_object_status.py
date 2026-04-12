from __future__ import annotations
import re
from pathlib import Path

BASE = Path('/mnt/data')

stage346 = (BASE / '5pn_stage346_actual_coherent_tracking_branch_values_output.txt').read_text()
stage349 = (BASE / '5pn_stage349_actual_four_condition_extractor_output.txt').read_text()
stage350 = (BASE / '5pn_stage350_family1_exact_branch_locator_output.txt').read_text()


def grab(text: str, label: str) -> str | None:
    m = re.search(rf"{re.escape(label)}\s*=\s*([^\n]+)", text)
    return m.group(1).strip() if m else None

results = {
    'Lambda_ell': grab(stage350, 'Lambda_ell = L/ell'),
    'kappa': grab(stage350, 'kappa'),
    'zeta_max': grab(stage350, 'zeta_max'),
    'Pe_star_chi': grab(stage350, 'Pe_*'),
    'zeta_phys_chi': grab(stage350, 'zeta_phys(Pe_*)'),
    'rho_alpha_max_chi': grab(stage350, 'rho_alpha,max'),
    'R_tr_formula': grab(stage346, 'R_tr'),
    'epsilon_formula': grab(stage346, 'epsilon'),
    'R_target_formula': grab(stage346, 'R_target'),
    'M_mix_formula': grab(stage346, 'M_mix'),
    'N_Q_formula': grab(stage346, 'N_Q on the natural source-map branch'),
}

print('STAGE 352 — FINAL NUMERICAL OBJECT STATUS')
print('=' * 72)
print()
print('Support/source branch is numerically located:')
for key in ['Lambda_ell', 'kappa', 'zeta_max', 'Pe_star_chi', 'zeta_phys_chi', 'rho_alpha_max_chi']:
    print(f'{key}: {results[key]}')
print()
print('Actual coherent branch finish packet remains symbolic in current stack:')
for key in ['R_tr_formula', 'epsilon_formula', 'R_target_formula', 'M_mix_formula', 'N_Q_formula']:
    print(f'{key}: {results[key]}')
print()
if 'No numerical PDE-selected point is present yet in the notes/scripts.' in stage349:
    print('Verdict: no numerical PDE-selected orbit-lock point is present yet in the current notes/scripts.')
else:
    print('Verdict: numerical PDE-selected orbit-lock point may be present; manual inspection still required.')
