# Moving-Throat PDE Framework Paper Support Manifest

## Purpose

This directory holds the publication-scale framework paper for the moving-throat
PDE program:

- `paper/pde.tex`
- `paper/pde.pdf`

It does **not** duplicate the full derivation archive.

## Policy

The authoritative derivation, audit, and provenance surface for this paper lives
in:

- `research/pde_ledger/`

That archive is the place to look for:

- theorem-block citation anchors,
- stage-level provenance,
- SymPy audits,
- Mathematica audits,
- review notes,
- and the longer derivation narrative.

This `research/pde/` directory intentionally has **no** local `scripts/` or
`mathematica/` subtree. That is by design, not an omission.

## Primary citation object

The companion citation target for this paper is:

- *Moving-Throat PDE Derivation Companion Archive Ledger*
- DOI: `10.5281/zenodo.19699523`

Repo entry points:

- `research/pde_ledger/paper/pde_ledger.tex`
- `research/pde_ledger/paper/pde_ledger_reader.tex`
- `research/pde_ledger/notes/CITATION_MAP.md`

## Framework-facing support surface

For the framework paper, the main stable ledger anchors are:

- `MTDC-T1`, `MTDC-T2` — geometry lift and breathing reduction
- `MTDC-T3` — BdG pole-shift packet
- `MTDC-T4`, `MTDC-T5` — normalization bridge, grouped bundle, isotropy
- `MTDC-T8.2` — exact `chi_Q = 1` compact-DtN matching

See:

- `research/pde_ledger/notes/CITATION_MAP.md`

## Representative supporting artifacts

Ledger-backed audits used by the framework paper:

- `research/pde_ledger/scripts/moving_throat_pde_master_sympy_audit.py`
- `research/pde_ledger/scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- `research/pde_ledger/scripts/moving_throat_pde_stage174_minimal_pde_data_packet_sympy_audit.py`
- `research/pde_ledger/scripts/moving_throat_pde_stage175_orbit_quotient_projectors_sympy_audit.py`
- `research/pde_ledger/scripts/moving_throat_pde_stage176_isotropic_grouped_p2_target_surface_sympy_audit.py`
- `research/pde_ledger/scripts/moving_throat_pde_stage206_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py`
- `research/pde_ledger/scripts/moving_throat_pde_stage207_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py`

Auxiliary prototype scripts referenced by the framework paper appendix:

- `scripts/5pn/5pn_stage3_isotropic_overlap_model.py`
- `scripts/5pn/5pn_stage4_axisymmetric_transport.py`
- `scripts/5pn/5pn_stage5_primitive_deformation_compensation.py`

These `scripts/5pn/` files are useful executable examples, but they are **not**
the primary citation surface. The ledger remains the authoritative provenance
wrapper.

## Packaging rule

If the framework paper changes its verification appendix or citation surface:

1. update `research/pde/paper/pde.tex`,
2. update this manifest,
3. do not copy ledger audits into `research/pde/`,
4. and prefer pointing to the ledger rather than creating parallel local trees.
