dry_run: true
non_binding: true

# Phase A Numeric-Literal Scan

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `003, 104, 105`
- Source files: `- stage: '003'
  role: paper_stage_tex
  path: paper/stages/stage_003.tex
- stage: '003'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage003_bdg_coupling.md
- stage: '003'
  role: sympy_script
  path: scripts/moving_throat_pde_stage003_bdg_sympy_audit.py
- stage: '003'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl
- stage: '104'
  role: paper_stage_tex
  path: paper/stages/stage_104.tex
- stage: '104'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
- stage: '104'
  role: sympy_script
  path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
- stage: '104'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl
- stage: '105'
  role: paper_stage_tex
  path: paper/stages/stage_105.tex
- stage: '105'
  role: notes_stage
  path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
- stage: '105'
  role: sympy_script
  path: scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py
- stage: '105'
  role: mathematica_script
  path: mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl`

Task:
Find numerical literals, rational coefficients, closed-form constants, or parameter assignments that could be fit-insertion-point candidates. Exclude pure identities such as harmonic normalizations, dimension labels, residual tolerances, line numbers, and pass/fail counters unless the source itself claims the value is matched, fixed, canonical, forced, or derived.

Emit only YAML:

```yaml
modality: numeric_literal
candidates:
  - candidate_key:
    anchor_stage:
    parameter_names: []
    citation:
      path:
      line:
      excerpt:
    reason:
```
