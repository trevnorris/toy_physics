# Moving Throat PDE — Verification Review Plan

## Purpose

This document defines the systematic verification of the moving-throat
restricted PDE derivation.  The derivation is **ongoing** — new stages will
be added as work continues.  Multiple AI agents independently review every
derivation step (notes) and its companion SymPy audit script, then record
findings in shared per-stage review files inside this folder.

The goal is **correctness**: catch sign errors, dropped terms, unjustified
assumptions, script/notes mismatches, and logical gaps before the derivation
is treated as settled.

---

## Scope

| Item | Path | Count |
|------|------|-------|
| Foundation docs | `notes/moving_throat/moving_throat_pde_{roadmap,phase0_spec,phase1_linearized_scaffold,unit_test_dn_branch}.md` | 4 |
| Stage notes | `notes/moving_throat/moving_throat_pde_stageNNN_*.md` | 166 (stages 001-103, 108-170; gap at 104-107) |
| SymPy audit scripts | `scripts/moving_throat/moving_throat_pde_stageNNN_*_sympy_audit.py` | 147 |
| Master audit script | `scripts/moving_throat/moving_throat_pde_master_sympy_audit.py` | 1 |
| Reference summaries | `notes/summaries/*.md` | 9 |

### Stage numbering note
Stages 104-107 do not exist in the repository.  The derivation jumps from
stage 103 (core parameter status) directly to stage 108 (positive source
theorem).

### Stages without a dedicated SymPy script
The following stages are **status / consolidation** notes with no independent
computation to audit.  They should still be reviewed for logical consistency
and correct summarization of prior results.

Foundations (no script expected): roadmap, phase0_spec, phase1_scaffold, unit_test  
Setup stages: 001, 002  
Status/consolidation stages: 032, 067, 073, 076, 079, 086, 096, 100, 103,
111, 115, 119, 124, 128, 132, 136, 140

---

## Review Protocol

Every reviewing agent **must** address each applicable item below.

### GROUND RULES — agent scripts

Agents **may** write and run short SymPy scripts to spot-check specific
claims.  LLMs are unreliable at mental arithmetic and symbolic algebra,
so using SymPy for independent verification is encouraged.

However, every script an agent runs **must actually compute something**.
The test is simple: does the script call `simplify()`, `expand()`,
`solve()`, `diff()`, `integrate()`, or otherwise perform real symbolic
computation?  Or does it just `print()` the agent's own reasoning
wrapped in Python comments?

**GOOD** — script that independently checks a specific equation:
```python
import sympy as sp
x, A, DK = sp.symbols("x A DK", positive=True)
# Derive the discriminant from scratch and compare to notes' claim
disc = (DK + (8/sp.pi**2 - 16/(9*sp.pi**2))*x)**2 + 4*x**2*(...)
notes_claim = ...
print("diff =", sp.simplify(disc - notes_claim))
```

**BAD** — script that just prints commentary:
```python
print("Symbol assumption check:")
print("  rho declared positive=True -- limits scope")
print("  No upper bound enforced -- acceptable")
# This computes NOTHING and verifies NOTHING
```

Rules:
1. Scripts must be **small and focused** — check one specific thing.
2. Scripts must **actually call SymPy** to compute, simplify, or solve.
3. Do NOT write scripts that are mostly comments and print statements
   restating your own analysis.  That is useless theater.
4. The pre-generated audit script outputs are in
   `scripts/moving_throat/output/` — read those rather than re-running
   the audit scripts themselves.

### A. Notes Derivation Review
1. **Equation-level correctness** — Check every numbered/displayed equation.
   Look for sign errors, missing factors, incorrect index contractions,
   wrong limits.
2. **Logical flow** — Does each claim follow from what precedes it (in this
   stage and in prior stages)?  Are intermediate steps skipped that are not
   obviously true?
3. **Assumptions** — Are all assumptions explicitly stated?  Are they
   justified by the parent equations (phase0_spec) or by prior stages?
4. **Completeness** — Are there edge cases, degenerate limits, or branches
   of the argument that are not addressed but should be?
5. **Notation consistency** — Are symbols used in the same way as in prior
   stages?  Flag any re-definitions or ambiguities.
6. **Physical interpretation** — Does the physical narrative make sense?
   Do limits (thin wall, large L, etc.) give expected behavior?

### B. SymPy Script Review (when script exists)
Read the **script source code** (not just the output).  The script is not
trusted — it is part of the derivation and must be verified independently.

1. **Faithful implementation** — Does the script set up exactly the
   equations / objects described in the notes?
2. **Correctness of code** — Are there coding bugs (wrong variable, typo in
   expression, incorrect simplify assumptions)?
3. **Hardcoded values** — Flag any numeric literals, magic numbers, or
   pre-computed expressions that should have been derived symbolically.
   A script that plugs in the final answer and checks `answer - answer == 0`
   proves nothing.  The script must *derive* results from inputs, not just
   *verify* a pre-stated identity.
4. **Tautological checks** — Could any assertion pass trivially?  For
   example: defining `X = expr` then checking `X - expr == 0`, or
   building both sides of an identity from the same intermediate so
   cancellation is guaranteed by construction.  Identify what each check
   actually tests vs. what it *appears* to test.
5. **Symbol assumptions** — Are positivity, reality, and commutativity
   assumptions on SymPy symbols correct?  Wrong assumptions can make
   `simplify()` return 0 when the expression is not actually zero.
6. **Output agreement** — Do the printed/asserted results match the claims
   in the notes?  Check the pre-generated output file in
   `scripts/moving_throat/output/`.
7. **Coverage** — Does the script test all the key results of the stage,
   or does it skip some?  Are there important equations in the notes that
   have no corresponding script check?

### C. Cross-stage consistency (for checkpoint stages)
At each checkpoint stage, also verify:
1. Do the accumulated results match what the checkpoint claims to summarize?
2. Are any earlier issues silently papered over?

---

## Batch Grouping & Review Order

Stages are grouped into 20 review batches.  Batches should be reviewed **in
order** because later batches depend on earlier results.  Within a batch,
stages can be reviewed in parallel.

Checkpoint stages (marked **CP**) deserve extra scrutiny — they are the
load-bearing joints where errors would propagate furthest.

### Batch 0 — Foundations
| Stage | File | Script | Notes |
|-------|------|--------|-------|
| — | `moving_throat_pde_roadmap.md` | — | Goal & scope |
| — | `moving_throat_pde_phase0_spec.md` | — | Frozen parent equations |
| — | `moving_throat_pde_phase1_linearized_scaffold.md` | — | Linearized problem setup |
| — | `moving_throat_pde_unit_test_dn_branch.md` | — | Finite-throat benchmark |

### Batch 1 — Geometry Lift & Coupling (001-009)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 001 | Geometry lift and linearized PDE skeleton | No |
| 002 | Breathing reduction back to (a,L) closure | No |
| 003 | Minimal BdG-wall coupling and first pole shifts | Yes |
| 004 | Localized Maxwell + mixed-sector reduction | Yes |
| 005 | Grouped real P2 bundle and outgoing-quadrupole bridge | Yes |
| 006 | Full grouped real P2 bundle, exact projectors | Yes |
| 007 | Explicit overlap integrals, O(3) isotropy theorem | Yes |
| 008 | Minimal isotropic normalization | Yes |
| 009 | Concrete axial overlaps | Yes |

### Batch 2 — Wall Profiles & Loading (010-019)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 010 | First nonconstant finite-throat wall/brane family | Yes |
| 011 | Loaded profile selection | Yes |
| 012 | Dynamic loading | Yes |
| 013 | Selected mode normalization | Yes |
| 014 | Selected branch reachability | Yes |
| 015 | Source map from mode integrals | Yes |
| 016 | Microscopic normalization equation | Yes |
| 017 | Softening depth normal form | Yes |
| 018 | Dimensionless normalization locus | Yes |
| 019 | Support feasibility frontier | Yes |

### Batch 3 — Continuum Kernel Extraction **[CP]** (020)
| Stage | Topic | Has Script |
|-------|-------|------------|
| **020** | **Exact continuum-kernel extraction of A, DeltaK_ax, beta_0, M_mix** | **Yes** |

### Batch 4 — Kernel Continuation (021-029)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 021 | Dimensionless continuum placement | Yes |
| 022 | Split u sector | Yes |
| 023 | Generalized selected branch | Yes |
| 024 | Rank-2 support completion | Yes |
| 025 | Rank-2 selected mode normalization | Yes |
| 026 | Support direction extraction | Yes |
| 027 | Continuum selected rank-2 closure | Yes |
| 028 | Coherent local tracking | Yes |
| 029 | Tracking branch bounds | Yes |

### Batch 5 — Coherent Kernel Map **[CP]** (030)
| Stage | Topic | Has Script |
|-------|-------|------------|
| **030** | **Coherent-kernel dimensionless map and support-enhancement factor** | **Yes** |

### Batch 6 — Support & Threshold Analysis (031-039)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 031 | Support compensation theorem | Yes |
| 032 | DtN overlap zeta | No (status) |
| 033 | Zeta threshold comparison | Yes |
| 034 | Lowest twin criterion | Yes |
| 035 | Non-twin asymmetry threshold | Yes |
| 036 | Overlap boost window | Yes |
| 037 | Robin softening support lane | Yes |
| 038 | Explicit lowest lane reachability | Yes |
| 039 | Transport source asymmetry | Yes |

### Batch 7 — Physical Parameter Placement **[CP]** (040)
| Stage | Topic | Has Script |
|-------|-------|------------|
| **040** | **Physical (Pe, kappa, eta) placement map for lowest support lane** | **Yes** |

### Batch 8 — Operator & Gain Analysis (041-052)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 041 | Coupled support-source operator | Yes |
| 042 | Operator branch residual bounds | Yes |
| 043 | Entropic microclosure | Yes |
| 044 | Microscopic gain thresholds | Yes |
| 045 | Parent action gain | Yes |
| 046 | Parent thresholds | Yes |
| 047 | Equilibrium alignment | Yes |
| 048 | Thin wall confinement | Yes |
| 049 | Wall figure of merit | Yes |
| 050 | Sech-Gaussian resonance | Yes |
| 051 | Resonance thresholds | Yes |
| 052 | Final reduced verdict | Yes |

### Batch 9 — Wall Branch & Family-1 Geometry (053-066)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 053 | GNLS wall shell | Yes |
| 054 | Tanh wall branch | Yes |
| 055 | Explicit branch thresholds | Yes |
| 056 | Family-1 geometry map | Yes (refresh) |
| 057 | Family-1 healing lock | Yes |
| 058 | Family-1 threshold window | Yes |
| 059 | N=5 wall depth lock | Yes |
| 060 | Family-1 theta extraction | Yes |
| 061 | Family-1 branch verdict | Yes |
| 062 | Family-1 quadrupole PE map | Yes |
| 063 | Family-1 zeta thresholds | Yes |
| 064 | Family-1 pi thresholds | Yes |
| 065 | Master quadrupole residual | Yes |
| 066 | Family-1 direct operator window | Yes |

### Batch 10 — Full Reduced PDE Write-Up **[CP]** (067)
| Stage | Topic | Has Script |
|-------|-------|------------|
| **067** | **Full reduced moving-throat PDE write-up skeleton** | **No** (consolidation) |

### Batch 11 — Loading Ratio & Isotropic Verdict (068-072)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 068 | Quadrupole demand cancellation | Yes |
| 069 | Family-1 loading ratio window | Yes |
| 070 | Outgoing branch loading ratio finish | Yes |
| 071 | Loading ratio from minimal module | Yes |
| 072 | Family-1 minimal isotropic verdict | Yes |

### Batch 12 — Geometry Lane (073-082)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 073 | Updated reduced status | No (status) |
| 074 | Grouped P2 static geometry derivation | Yes |
| 075 | Dynamic geometry obstruction | Yes |
| 076 | Grouped P2 status update | No (status) |
| 077 | Isotropic geometry decoupling | Yes |
| 078 | Second-order geometry contamination | Yes |
| 079 | Geometry lane check verdict | No (status) |
| 080 | Single normalization defect | Yes |
| 081 | Family-1 support is automatic | Yes |
| **082** | **Reduced finish line** | **Yes** **[CP]** |

### Batch 13 — Outgoing DtN (083-089)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 083 | Outgoing normalization factorization | Yes |
| 084 | Natural source map reduction | Yes |
| 085 | Higher odd irrelevance | Yes |
| 086 | Reduced 2.5PN conditional closure | No (status) |
| 087 | Outgoing DtN fingerprint | Yes |
| 088 | chiQ fix from outgoing DtN | Yes |
| **089** | **Canonical outgoing reduced closure** | **Yes** **[CP]** |

### Batch 14 — General DtN & Core Extraction (090-103)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 090 | General DtN deformation | Yes |
| 091 | Robustness classes | Yes |
| 092 | Linearized branch selection | Yes |
| 093 | Robin outlet model | Yes |
| 094 | Mixed sidechannel pole | Yes |
| 095 | Hybrid Robin-mixed compensation | Yes |
| 096 | Outlet model status | No (status) |
| 097 | Concrete core Schur | Yes |
| 098 | Core balance compensation | Yes |
| 099 | DtN mixed tube realization | Yes |
| **100** | **Outlet core status** | **No** (status) **[CP]** |
| 101 | Parent core extraction | Yes |
| 102 | Parent balance family | Yes |
| 103 | Core parameter status | No (status) |

### Batch 15 — Positive Source & Mouth Dynamics (108-119)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 108 | Positive source theorem | Yes |
| 109 | Positive source families | Yes |
| 110 | Penetration families | Yes |
| 111 | Mouth source bias status | No (status) |
| 112 | Mouth boundary layer | Yes |
| 113 | Mouth bias map | Yes |
| 114 | Parent mouth threshold | Yes |
| 115 | Mouth boundary layer status | No (status) |
| 116 | Coupled mouth fixedpoint | Yes |
| 117 | Family-1 mouth fixedpoint | Yes |
| 118 | Outlet consistent mouth closure | Yes |
| 119 | Coupled mouth status | No (status) |

### Batch 16 — Core-to-Mouth Gain (120-128)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 120 | Core-to-mouth gain map | Yes |
| 121 | Normalized mouth gain family | Yes |
| 122 | Family-1 actual mouth gains | Yes |
| 123 | Self-matched mouth susceptibility | Yes |
| 124 | Mouth gain status | No (status) |
| 125 | Self-consistent mouth branch | Yes |
| 126 | Equal normalized singular limit | Yes |
| 127 | Unique regular canonical branch | Yes |
| 128 | Mouth branch selection status | No (status) |

### Batch 17 — Rigidity & Corrections (129-140)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 129 | Positive deformation expansion | Yes |
| **130** | **First-order rigidity kernel at canonical family-1 point** | **Yes** **[CP]** |
| 131 | Representative positive families | Yes |
| 132 | Mouth rigidity status | No (status) |
| 133 | Full profile residual | Yes |
| 134 | First-order selected correction | Yes |
| 135 | Family-1 actual correction | Yes |
| 136 | Full mouth correction status | No (status) |
| 137 | Coevolving core-mouth map | Yes |
| 138 | Frozen traction fixedpoint | Yes |
| 139 | Renormalized canonical branch | Yes |
| 140 | Core-mouth coevolution status | No (status) |

### Batch 18 — Linear Defect Transport & Final Result (141-152)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 141 | Linear defect transport | Yes |
| 142 | Hybrid outlet projection | Yes |
| 143 | Bare mixed port slippage | Yes |
| 144 | DtN similarity slippage | Yes |
| 145 | Parent compensation rigidity | Yes |
| 146 | Off-family normal coordinate | Yes |
| 147 | Microscopic log channels | Yes |
| 148 | Exact branch drifts | Yes |
| 149 | Bundle inversion four drifts | Yes |
| **150** | **Bundle transport tangent-compensation theorem** | **Yes** **[CP]** |
| 151 | Off-bundle slippage | Yes |
| 152 | No linear grouped-P2 feed-down into scalar off-bundle slippages | Yes |

### Batch 19 — Grouped Outlet & Similarity Closure (153-169)
| Stage | Topic | Has Script |
|-------|-------|------------|
| 153 | Linear grouped outlet map | Yes |
| 154 | Microscopic grouped obstructions | Yes |
| 155 | Physical slope collapse | Yes |
| 156 | Axisymmetric loading mismatch | Yes |
| 157 | Static self-similarity | Yes |
| 158 | Wall-normalized load shape | Yes |
| 159 | Outgoing load factorization | Yes |
| 160 | Weak axisymmetric outgoing slippage | Yes |
| 161 | Outgoing port co-loading theorem | Yes |
| 162 | Transfer shape theorem | Yes |
| 163 | Effective transfer shape collapse | Yes |
| 164 | Coherent tracking defect | Yes |
| 165 | Microscopic coherent slippage | Yes |
| 166 | Triangular normal form | Yes |
| 167 | Branch invariant coordinates | Yes |
| 168 | Microscopic monomials | Yes |
| 169 | Similarity orbit closure | Yes |
| **170** | **Orbit quotient closure** | **Yes** **[FINAL]** |

---

## Review File Format

Each stage has a review file in this folder:

```
review/stage_NNN_review.md       (e.g. stage_020_review.md)
review/foundations_review.md     (for the 4 foundation docs)
```

All reviewing agents write into the **same file** for a given stage.  The
file has a standard template (pre-populated in each stub).  Each agent
appends a review section under `## Agent Reviews`.

### Verdict codes
| Code | Meaning |
|------|---------|
| **PASS** | No issues found; derivation and script are correct |
| **MINOR** | Cosmetic or notation issues only; no mathematical errors |
| **ISSUE** | One or more substantive errors or gaps found (detail required) |
| **BLOCK** | Critical error that invalidates downstream stages |

A stage is considered **verified** when at least 2 independent agents give
PASS or MINOR with no unresolved ISSUE/BLOCK findings.

---

## Process for Reviewing Agents

1. **Read the stage notes** in full.  Also read the immediately prior stage
   to understand what is being carried forward.
2. **Read the SymPy script source code** (if it exists) line by line.
3. **Read the pre-generated script output** from
   `scripts/moving_throat/output/`.  Do NOT re-run the audit script.
4. **Work through the derivation** step by step.  Check every equation
   against protocol items A.1–A.6.
5. **Compare script to notes** per protocol items B.1–B.7.
6. **Write short focused SymPy checks** if you need to independently
   verify a specific equation, sign, or identity.  See ground rules
   above — the script must actually compute something.
7. **Write your review** into the stage's review file, following the
   template.  Use your model name and the date.
8. **Do not alter another agent's review.**  Append yours below the
   separator line.

### For checkpoint stages (CP)
Also verify cross-stage consistency (protocol C).  Read the preceding batch
of stages and confirm the checkpoint accurately consolidates them.

### For status stages (no script)
Focus on whether the status summary faithfully reflects the results of the
stages it references.  Check for any claims that were not actually
established.

### SymPy script outputs
Pre-generated script outputs are stored in `scripts/moving_throat/output/`.
Each file is named `moving_throat_pde_stageNNN_*_sympy_audit.txt` and contains
a metadata header (date, elapsed time, exit code, status) followed by the full
stdout/stderr of the script.

**Agents should read the output file rather than re-running the script.**
If the output file shows FAIL or TIMEOUT, that is itself a finding to report.

To regenerate outputs (e.g. after script edits):
```bash
bash scripts/moving_throat/run_all_audits.sh           # incremental (only changed)
bash scripts/moving_throat/run_all_audits.sh --force    # re-run all
```

### Reference materials
Agents may consult:
- `notes/summaries/*.md` — paper summaries for the parent theory
- `notes/moving_throat/moving_throat_pde_phase0_spec.md` — frozen parent equations
- `scripts/moving_throat/moving_throat_pde_master_sympy_audit.py` — master audit
- `scripts/moving_throat/output/_summary.txt` — pass/fail summary for all scripts

---

## Tracking

Progress is tracked by the state of the review files themselves.  A stage
is in one of three states:

- **Pending** — review file exists but has no agent reviews yet
- **In Review** — one or more agents have written reviews, but < 2 PASS/MINOR
- **Verified** — at least 2 independent PASS/MINOR verdicts with no open ISSUE/BLOCK

### Quick-check commands
To see which stages still need review:
```bash
grep -rL "### Agent:" notes/moving_throat/review/stage_*_review.md
```
To see which stages have unresolved issues:
```bash
grep -l "ISSUE\|BLOCK" notes/moving_throat/review/stage_*_review.md
```

---

## Extending This Plan

The derivation is ongoing.  When new stages are added:

### Adding new stage review files

Run the helper script to generate stubs for any new stages that don't yet
have review files:

```bash
scripts/moving_throat/generate_review_stubs.sh
```

This script scans `notes/moving_throat/` for stage files, checks which ones
already have a review stub, and creates missing ones.

### Adding new batches

If new stages extend beyond the current batch 18 range (stages 141-152),
or if a new logical phase begins:

1. Add a new `### Batch NN` section to the table above.
2. Update the batch lookup in `generate_review_stubs.sh`.
3. Update the "Scope" table counts at the top of this file.

### Renumbering / inserting stages

If stages are inserted between existing ones (e.g., stage 152a), use a
suffix convention for the review file: `stage_152a_review.md`.  The helper
script handles this.

### Key principle

The review folder is the single source of truth for verification status.
Keep the plan in sync with the actual stage files, and the tracking commands
above will continue to work.
