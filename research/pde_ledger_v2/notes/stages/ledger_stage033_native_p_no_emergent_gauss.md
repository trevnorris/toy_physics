# ledger_stage033_native_p_no_emergent_gauss

## Status

**Part IV — Charge. IV-4 (build-order 033; the LAST stage of the 4-stage Part-IV split,
user decision 2026-07-22).** The self-contained Dirac–Bergmann constraint-class analysis of the
native-`P` medium's own field content. Where 030–032 built and landed the *throat*-charge mechanism
(the static ±w substrate, the puncture-deflection response, the BC-ensemble sign R1), stage 033 asks the
complementary question of the *candidate self-hosting electric theory* — can the medium's own `P`-field
Lagrangian carry an **exact U(1)/Maxwell Gauss law** on its own? — and records the honest **NO**:

- **`NATIVE_P_NO_EMERGENT_GAUSS`** — exact U(1)/Maxwell Gauss is **non-native** to the `P`-field medium. At
  **quadratic order** the two native families **THEORY-A** and **THEORY-C** are Dirac–Bergmann
  **all-second-class** (every constraint `SECOND_CLASS_COMPONENT`; weak PB matrix full rank; **FC = 0**; **no
  first-class Gauss chain**), and — within the argued boundary partition + fixed-seed scanned representatives —
  the tuned rank-drop strata that *do* carry first-class directions have a **zero Gauss descendant**
  (`DESCENDANT_ZERO`): the FC direction cannot descend to a `∝k` Gauss law. Six able-to-fail controls anchor
  the machinery. Verdict token **`NATIVE_P_NO_EMERGENT_GAUSS`** (both A and C, both engines, both exit 0).

**⚠ CRITICAL TYPING (do not overstate).** This is a **CHARACTERIZED-DEPARTURE, recorded first-class** — a
load-bearing no-go, NEVER softened, NEVER "rescued." It is **NOT** a knob, **NOT** a reduction, **NOT** an R1
deferral. It discharges **no** knob and **does not shrink the irreducible count** (`docs/model_map.md` §4, the
honest departure ledger). Reads out as: **"EM is NOT exact Maxwell."** The named able-to-fail alternates that
the verdict would have re-derived to had a first-class Gauss chain been found are
**`FIRST_CLASS_GENERIC_EM_CANDIDATE`** (from the symbolic open stratum) / **`FIRST_CLASS_TUNED_INVERSE_DESIGN`**
(from the tuned locus).

Ledger-local landing label (surviving-solution framing): the native-`P` medium does not **generically** host an
emergent Maxwell Gauss law; a hypothetical missed measure-zero first-class stratum could only be a **TUNED /
inverse-design** artifact, so the **generic** no-go stays decisive.

**⭐ This is the LAST stage of Part IV.** After 033, charge is complete: 030 paid for the static substrate, 031
dressed it with the target-blind puncture-deflection mechanism, 032 landed the honest sign R1, and 033 records
the first-class departure that the whole EM sector inherits. This is the **Stage-1 constraint gate only** —
compactness, charge quantization, deconfinement, and native ±w current supply remain **downstream** (do NOT
claim them).

## Purpose

Record, as a first-class characterized departure, that the medium's own `P`-field content
`{p = P_∥, u, s, b, …}` cannot support exact Maxwell electromagnetism as a native gauge structure. The decisive
object is the Dirac–Bergmann constraint algebra of each native family at quadratic order: a first-class Gauss
generator would show up as (i) a nonzero first-class count in the weak Poisson-bracket matrix, and (ii) a
first-class primary whose Hamiltonian descendant is a `∝k` (Gauss-shaped) spatial action. Both engines run each
family — and six named controls — through the **identical** pipeline

```text
input Lagrangian  →  Legendre transform (primaries only from Hessian nulls)
                  →  Dirac closure (iterate constraints to closure)
                  →  PB-kernel first-class / Gauss search.
```

and find, for both native families, FC = 0 with `gauss_candidates = 0` on the regular open kinetic stratum, and
`DESCENDANT_ZERO` on every first-class direction of the tuned rank-drop strata. The six controls prove the
machinery is not blind: it finds a genuine first-class Gauss chain exactly when one exists (Maxwell), and
returns the correct distinct class token in five other decidable cases.

Consumes **nothing numeric** from 030/031/032 — it is an independent constraint analysis (see "Consumes / cites"
below). It is the charge-sector twin of the light-sector departure **`FAIL_CAUCHY_STRAY_LONGITUDINAL`** (Part
III / stage 003 / pathA_36): "a second-class pair, NOT Maxwell's first-class Gauss."

## 1. The two native families and the shared pipeline

Each native family is supplied as a quadratic Lagrangian in the medium's own fields and reduced through the
single pipeline above; nothing is pre-diagonalized or hand-fed a constraint. The pipeline builds the momentum
map `Π_i = ∂L/∂v_i`, the velocity Hessian `H_ij = ∂²L/∂v_i∂v_j`, the primary constraints from the Hessian left
null vectors, the representative Hamiltonian by inverting the regular Hessian minor, then iterates Dirac
consistency (`{Φ, H}` on the primary-bracket kernel) to closure, and finally searches the weak Poisson-bracket
matrix kernel for a first-class primary that descends to a `∝k` Gauss action.

**THEORY-A** carries the field set `{p₁₋₃, u₁₋₃, s, λ, σ, b}` with retained couplings
`{g_tA, g_sA, g_dA, g_bA}` (a scalar multiplier `λ` and an axial `σ` enforce `s` and `K·u`, and `b` couples via
`g_bA b K·p`). **THEORY-C** carries `{p₁₋₃, u₁₋₃, s, λ, σ, b, ξ}` with retained couplings `{g_tC, g_sC, g_dC}`
(an extra multiplier `ξ` enforces `b`, so `g_bC` is absent). In both, the kinetic core is
`½(ṗ·ṗ + ṡ² + ρ_u u̇·u̇) + g_t ṗ·u̇` (A adds `½ḃ²`), and the potential carries the transverse `u` operator
`K²I − KKᵀ`, the mixing `g_s K² p·u`, and the longitudinal `½ g_d (K·p)²`.

## 2. THEORY-A signature (`PASS_THEORY_A_SIGNATURE`, `PASS_GUARD_COUPLINGS_ENTER_A`)

On the regular open kinetic stratum (all retained couplings symbolic, `g_tA² ≠ ρ_uA`):

```text
Hessian rank / nullity        = 8 / 2
primary + secondary constraints = 8   with stages [0, 0, 1, 1, 2, 2, 3, 3]
constraint classes            = all SECOND_CLASS_COMPONENT (×8)
weak PB matrix                : rank 8,  FIRST_CLASS 0,  SECOND_CLASS 8
gauss_candidates              = 0,  additional_G_exists = False
```

The full signature tuple asserted (both engines) is
`(8, 2, 8, [0,0,1,1,2,2,3,3], SECOND_CLASS×8, 8, 0, 8, 0, False, True)`. The Dirac tower closes in four stages;
every constraint is second-class; the weak PB matrix is nonsingular (rank = constraint count 8 ⇒ FC = 0). No
first-class primary exists, so the Gauss search returns zero candidates.

**`GUARD-COUPLINGS-ENTER` (A) PASS.** Every native coupling `{g_tA, g_sA, g_dA, g_bA}` provably survives into the
computed momentum map / Hessian / constraints / PB matrix (none silently drops); dropping any one at input (e.g.
`g_bA → 0`) fires the guard.

## 3. THEORY-C signature (`PASS_THEORY_C_SIGNATURE`, `PASS_GUARD_COUPLINGS_ENTER_C`)

On the regular open kinetic stratum:

```text
Hessian rank / nullity        = 7 / 4
primary + secondary constraints = 12  with stages [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3]
constraint classes            = all SECOND_CLASS_COMPONENT (×12)
weak PB matrix                : rank 12,  FIRST_CLASS 0,  SECOND_CLASS 12
gauss_candidates              = 0,  additional_G_exists = False
```

Signature tuple `(7, 4, 12, [0,0,0,0,1,1,1,1,2,2,3,3], SECOND_CLASS×12, 12, 0, 12, 0, False, True)`. Same
conclusion: full-rank weak PB matrix, FC = 0, no Gauss candidate. **`GUARD-COUPLINGS-ENTER` (C) PASS** over
`{g_tC, g_sC, g_dC}`.

The `FC = n_constraints − rank` entry of each signature tuple is a **true-by-construction rank–nullity
identity** (it always holds); the load-bearing signature elements — mutated by `PASS_THEORY_{A,C}_SIGNATURE`,
which give a multiplier a kinetic term and recompute — are the Hessian rank/nullity, the constraint count and
stages, the all-second-class classification, and the `FC 0 / SC` counts.

## 4. Kinetic-Hessian determinant and the sole tuned degeneracy (`PASS_KINETIC_HESSIAN_DETERMINANT`)

The physical (6×6) kinetic-Hessian determinant factors, in **both** families, as

```text
det H_phys = (ρ_u − g_t²)³        (per (p_i, u_i) component:  ρ_u − g_t²).
```

So the **only** additional kinetic degeneracy — the only surface on which the open-stratum all-second-class
result can change — is the tuned locus **`g_t² = ρ_u`**. Off it, the pipeline is on the regular open stratum
(§§2–3, FC = 0). Perturbing the `u`-kinetic input so `det ≠ (ρ_u − g_t²)³` fires the tooth.

## 5. The `DESCENDANT_ZERO` hardening guard on the tuned locus (`PASS_HARDENING_DESCENDANT_{A,C}`, `PASS_PRIMARY_FC_ACCOUNTING_{A,C}`)

On the tuned degeneracy surface `g_t² = ρ_u` the Hessian nullity increases (A: 2 → **5**; C: 4 → **7**) and
first-class directions can appear. The potential-null residual there is solved **symbolically** (not hardcoded)
for the two `g_t = ±1` signs; within the argued boundary partition + scanned representatives, the **common
Hessian/potential-null locus** is the observed stratum that develops first-class directions (**FC = 2** at that
locus). For every first-class primary that the search rejects, both engines compute the Hamiltonian descendant
`{primary, H₂}` and its field action `{q, {primary, H₂}}`, and require it be certified either
**`DESCENDANT_ZERO`** (the descendant is identically `0`) **or** non-gradient (no `∝k` spatial action). The
guard **fails the run** if any rejected direction hides a nonzero `∝k` Gauss descendant — so a tuned FC
direction cannot masquerade as an emergent Gauss law.

Reproduced audit (both engines, both families):

```text
FC-bearing strata checked      = 2   per family
rejected directions checked    = 4   per family
every rejected direction       : descendant_zero = True,
                                 secondary_action_non_gradient = True,
                                 descendant_rejection_certified = True,
                                 reason = DESCENDANT_ZERO
accounting:  primary_fc == rejected + candidates   (candidates = 0)
```

The accounting tooth (`PASS_PRIMARY_FC_ACCOUNTING_{A,C}`) enforces `primary_fc == rejected + candidates`; the
descendant tooth injects a genuine Maxwell FC primary (found live by the same engine) with a nonzero `∝k`
descendant and confirms it fires at its own assert.

## 6. The six able-to-fail controls (`PASS_CONTROL_*`, `PASS_GUARD_SEARCH_CAPABLE`)

Each control is an input Lagrangian through the **same** pipeline; the class token and `(FC, SC, Gauss cand.,
Hessian nullity)` are computed at runtime and compared against the expected tuple, and each has its **own**
ablation that fires `FIRED_AT_OWN_ASSERT`:

| # | Control | Expected class token | FC | SC | Gauss cand. | Hess. nullity |
|---:|---|---|---:|---:|---:|---:|
| 1 | maxwell | `FIRST_CLASS_GAUSS` | 2 | 0 | 1 | 1 |
| 2 | gauged hard unit | `MIXED` | 2 | 4 | 1 | 2 |
| 3 | bare O(4) hard sigma | `SECOND_CLASS_RADIAL_NO_GAUSS` | 0 | 4 | 0 | 1 |
| 4 | nonconserved current | `INCONSISTENT_PRESERVATION` | 2 | 0 | 1 | 1 |
| 5 | Coulomb-gauge Maxwell | `SECOND_CLASS_NO_LOCAL_GAUGE` | 0 | 8 | 0 | 3 |
| 6 | global U(1) complex scalar | `GLOBAL_CHARGE_NO_LOCAL_GAUSS` | 0 | 0 | 0 | 0 |

**`GUARD-SEARCH-CAPABLE` PASS.** Maxwell and the gauged-hard-unit control must yield **nonzero** Gauss
candidates — so the search demonstrably finds a real first-class Gauss chain when one exists, and its absence
for native A/C is meaningful, not vacuous. Control 4's Gauss preservation carries the explicit continuity defect

```text
{Gauss, H}  =  −nc_j1 − 2·nc_j2 − 3·nc_j3 + nc_rho_dot        (≠ 0 without a conservation rule)
```

which is nonzero absent a continuity rule ⇒ `INCONSISTENT_PRESERVATION`; imposing conservation (the ablation)
flips it to `CURRENT_CONSISTENT` and fires the tooth.

## 7. Honest scope (MANDATORY — it does NOT weaken the verdict) (`PASS_HONEST_TUNED_SCOPE`, `PASS_DECISION_ORDER_BRANCH2`, `PASS_VERDICT_TOTALITY`)

Separate the two coverage classes:

- **(a) DECISIVE — the fully-symbolic open kinetic stratum.** FC = 0 for A and C for **all** retained
  couplings (the whole regular locus `g_t² ≠ ρ_u`). This is a symbolic result, not a scan.
- **(b) ARGUED + FIXED-SEED SCANNED — the tuned `g_t² = ρ_u` locus.** NOT an exhaustive symbolic
  stratification. The scanned coverage is two parts:
  - a **6-point semidefinite-boundary scan per family** — 2 signs `g_t = ±1` × {generic semidefinite, first
    rank-drop non-common representative, common Hessian/potential-null locus}; the common-null representative is
    the FC = 2 one, all with `gauss_candidates = 0` (⇒ `DESCENDANT_ZERO`); nullity 5 (A) / 7 (C) throughout;
  - a **fixed-seed randomized rank-drop sweep of 6 points per family, 12 total across A + C** (seeds
    **A = `260713`, C = `260715`**; per-sign split **3 + 3**), all `FC = 0, gauss_candidates = 0,
    additional_G_exists = False`, nullity 5 (A) / 7 (C).

Pinned in the audit are the seed, the sign split, the cardinality (6 boundary / family, 12 randomized total),
and the **decisive per-point signatures** (`FC = 0, gauss = 0, additional_G = False`) — NOT the specific sampled
coordinates (a different RNG may sample different points). A hypothetical missed measure-zero first-class Gauss
stratum would be a **TUNED / inverse-design** artifact (`missed_measure_zero_class = TUNED_INVERSE_DESIGN`), so
the physical no-go — "native `P` does not **generically** host emergent EM" — stays decisive
(`generic_no_go_decisive = True`). The scan is **not** promoted into an exhaustive classification; there is no
unqualified-exhaustiveness "only."

**Decision order (branch 2).** With FC = 0 for both families **on the open stratum** at quadratic order, the
audit records `BRANCH_2_QUADRATIC_ABSENCE`: a regular nonlinear gauge identity would have a nontrivial leading
linearization, so the quadratic absence is the decisive branch. 033 records the **quadratic** verdict + this
decision-order argument; it does **not** execute an all-orders proof. The verdict is re-derived from the
computed searches (`PASS_VERDICT_TOTALITY`): the criterion is **no Gauss candidate**
(`additional_G_exists = False`) — with no first-class primary descending to a `∝k` Gauss action in either the
open or the tuned searches it lands `NATIVE_P_NO_EMERGENT_GAUSS`. This is the crux of the no-go: the tuned
common-null points **DO develop FC = 2**, but those first-class directions have a **zero Gauss descendant**
(`DESCENDANT_ZERO`), so no emergent Gauss law appears — the verdict follows from the absent Gauss candidate, NOT
from an absence of first-class directions. Feeding a genuine first-class Maxwell search (whose FC direction DOES
descend to a `∝k` Gauss action) through the open path re-derives
`FIRST_CLASS_GENERIC_EM_CANDIDATE` instead (the able-to-fail witness — the tooth reads the computed result, not
a duplicated literal, per the 030 `X≡X` lesson).

## 8. Departure typing (recorded honestly)

- The verdict is a **CHARACTERIZED-DEPARTURE, first-class** (`docs/model_map.md` §4). It is a load-bearing
  no-go, kept clean and simple, NEVER softened or rescued.
- It is **not a knob and discharges no knob**; it is **not a reduction / codimension edge**; it **does not
  shrink the irreducible count**. Registered as edge **R66** (a departure edge that discharges nothing) plus a
  **structural-constant / departure-support** row (the constraint-class signature + the `g_t² = ρ_u` locus +
  the six-control decision table — dimensionless structural facts, no `[L, T, M]` to reduce).
- It is the **quadratic Stage-1 constraint gate** only; compactness, quantization, deconfinement, and native
  ±w current supply are downstream and are **not** claimed here.

## 9. Source-to-stage predicate manifest

Completeness certificate: **no silently-dropped source claim.** Every source tooth / Q-claim in scope from the
source build `software/em_charge_attribute/reports/native_p_constraint_gate.md`
(+ `native_p_gate_sympy.py` / `native_p_gate_dual.wl`) lands as **PRESERVED** (folded as-is), **REPLACED_BY_STRONGER**
(a stronger reconstruction tooth), or **SCOPED_OUT** (harness / file-I/O contract, with reason). The partition
is disjoint + exhaustive, computed at runtime in **both** engines from the same canonical
`(id, disposition, owner)` triples:

```text
partition = { PRESERVED: 33, REPLACED_BY_STRONGER: 15, SCOPED_OUT: 3 }   (51 total),
manifest digest (SHA-256) = 6b191e77fefe24c9000445f01e4e2c6154ab1bb9b15bb40a6d1515dc463a7e9d.
```

Both engines assert the identifier set (51 unique), the three-way disposition set, the `STAGE033_*` owner
prefix on every row, the exact counts, and the committed digest at runtime and agree (`PASS_SOURCE_PREDICATE_MANIFEST`;
dropping one row fires it). Representative dispositions (matching the engine `SOURCE_MANIFEST` exactly):
**PRESERVED (33)** — the `Q1`/`Q2`/`Q5`/`Q6` closed-rebuild claims (Legendre, Maxwell-computed, Gauss-search,
Dirac-closure), the two `GUARD-COUPLINGS-ENTER`, `GUARD-SEARCH-CAPABLE`, the two `HARDENING_DESCENDANT`, the
per-family Hessian-rank/nullity + constraint-count + stages + all-second-class + Gauss-absence claims, the six
control classifications (`CONTROL_*`), the boundary-scan / randomized-sweep / fixed-seed claims, and the
source-first / shear-duplicate / decision-order-branch2 bookkeeping; **REPLACED_BY_STRONGER (15)** — the `Q3`
six-controls claim (→ native ablations), the `Q4` Wolfram-independence rebuild, the two per-family PB
rank/FC/SC strengthenings, the determinant `(ρ_u − g_t²)³` upgrade, the six per-tooth control ablations
(`ABLATION_*`, owner `STAGE033_ENV_OWN_ASSERT`), the two symbolic common-null solves (`COMMON_NULL_A`/`_C`), the
argued+scanned honest-scope reframe, and the computed-verdict totality; **SCOPED_OUT (3)** — the
`argparse`/`--out-dir` harness, all JSON artifact writing, and the cross-engine file comparator (the print-only /
zero-file-I/O / independent-tokens contract).

## 10. Consumes / cites

- **Consumes NOTHING numeric from 030/031/032.** 033 does not touch the puncture-deflection response matrix `m`,
  `m_gg`, the four `A_X`, the kernel determinant `D`, or the G0 zero-mode kernel — it is an independent
  Dirac–Bergmann constraint-class analysis of the native-`P` medium Lagrangian families A and C.
- **Cites (provenance, not consumed):** the Part IV charge context — native `P` = the medium's own field content
  `{p = P_∥, u, s, b, …}`, the candidate self-hosting electric theory, complementary to the throat-charge
  mechanism of 030–032 — and the sibling light-sector departure **`FAIL_CAUCHY_STRAY_LONGITUDINAL`** (Part III /
  stage 003 / pathA_36: "a second-class pair, NOT Maxwell's first-class Gauss"). Both are the model's honest
  "EM is not exact Maxwell" departures in `docs/model_map.md` §4.

## Verification

- **Dual-engine, both exit 0, genuinely independent routes.**
  `scripts/ledger_stage033_native_p_no_emergent_gauss_sympy_audit.py` — **SymPy 33 teeth**.
  `mathematica/ledger_stage033_native_p_no_emergent_gauss_mathematica_audit.wl` — **Mathematica 33 teeth**.
  Standalone, print-only, assert-zero (`raise SystemExit(1)` / `Exit[1]`), no argparse harness, no JSON/YAML
  payload, **zero file-I/O between engines**; float-/machine-real-free payload throughout. Each engine reaches
  `NATIVE_P_NO_EMERGENT_GAUSS` (A and C) on its own and prints its own tokens — cross-engine agreement is that
  they independently produce the same tokens, not a compare pass.
- **The `.wl` is a genuinely INDEPENDENT route** — re-authored, not a mirror of the `.py`: it derives the
  common-null locus by native `Solve` / `Reduce[ForAll]` (rather than hardcoding the tuning), uses native
  primitives (`Solve`, `Reduce`, `NullSpace`, `Inverse`, `Det`, `MatrixRank`), and carries its **own** six
  per-tooth control ablations.
- **Per-tooth ablation** (env switch `LEDGER_STAGE033_MUTATION`): every folded predicate has a mutation that
  fires at its **own** assert, including the six control ablations (`FIRED_AT_OWN_ASSERT`, both engines), the
  determinant / Legendre / Dirac-closure teeth, the `GUARD-COUPLINGS-ENTER` / `GUARD-SEARCH-CAPABLE` /
  `HARDENING_DESCENDANT` guards, and the verdict tooth (which re-derives from a mutated **computed** search
  input to a NAMED different verdict, never the `NATIVE_P_NO_EMERGENT_GAUSS` literal).
- **Tri-review outcome (falsification-first — recorded transparently, not hidden).** Fidelity **CLEAN**;
  adversarial **CONCERNS(2)** + a fidelity **NIT** — all three VERIFIED non-blocking, leaving **no real coverage
  gap**:
  1. **`PASS_HONEST_TUNED_SCOPE` is a FRAMING/summary tooth** — its mutation flips a self-set flag
     (`exhaustive_tuned_stratification`). Its substantive content (the 6-point boundary signature, the 12-point
     FC = 0 sweep, the symbolic-open FC = 0) is genuinely enforced by `PASS_BOUNDARY_SCAN_{A,C}`,
     `PASS_RANDOMIZED_SWEEP_{A,C}`, and `PASS_THEORY_{A,C}_SIGNATURE` (computed-data teeth). So there is **no
     coverage gap**, but the framing tooth itself is **redundant** (a known robustness item, not remediated).
  2. **The SymPy engine uses a separate `focused_mutation_probe` path under mutation** (vs the `.wl`'s inline
     injection). Adversarially traced: the probe re-derives from the **same** pipeline objects (`native_model`,
     `build_hamiltonian`, `dirac_closure`, `gauss_search`, `control_result`, `common_null_derivation`) — no
     hidden vacuity; a **verified-safe structural smell**, not a defect.
  3. **The signature's `FC = n_constraints − rank` element is a true-by-construction rank–nullity identity**
     (it cannot fail on its own). The real FC/SC/count guards are the **other** signature elements, mutated by
     `PASS_THEORY_{A,C}_SIGNATURE` — a **fidelity NIT**, not a coverage hole.

  Arbiter re-runs of both engines reproduce the tokens and the manifest digest; the tri-review leaves the
  departure decisive.

## Downstream consumers

- **Part IV is COMPLETE** after 033 (030 substrate → 031 mechanism → 032 sign R1 → 033 departure).
- **Parameter register:** edge **R66** (the `NATIVE_P_NO_EMERGENT_GAUSS` departure — discharges NO knob, NOT a
  reduction, does NOT shrink the irreducible count) + a structural-constant / departure-support row (the
  constraint-class signature `{A: rank 8, FC 0, SC 8; C: rank 12, FC 0, SC 12}`, the `g_t² = ρ_u` locus, the
  six-control decision table — "structural, dimensionless"). Cross-referenced to the sibling departure
  `FAIL_CAUCHY_STRAY_LONGITUDINAL` (Part III / pathA_36).
- **`docs/model_map.md` §4:** the `NATIVE_P_NO_EMERGENT_GAUSS` departure bullet (already present) — kept
  consistent: exact U(1)/Maxwell Gauss proven non-native (Dirac–Bergmann second-class at quadratic order; no
  first-class chain).
- **Part VII:** enters the honest departure ledger alongside the light-sector twin — a first-class recorded
  no-go, not a reduction credit.

## Provenance

- **Physics source:** `software/em_charge_attribute/reports/native_p_constraint_gate.md` (the rebuilt native-`P`
  constraint gate) + `software/em_charge_attribute/native_p_gate_sympy.py` +
  `software/em_charge_attribute/native_p_gate_dual.wl` (the source `.wl` mirrored the `.py`; the stage `.wl` is
  RE-AUTHORED as a materially distinct Wolfram route). The file-I/O comparator `native_p_gate_compare.py` is
  STRIPPED (zero file-I/O contract).
- **Consumes:** nothing numeric from 030/031/032 (independent constraint analysis).
- **Governing:** `notes/ledger_v2_blueprint.md` §5 (reshape spec) + §6 (per-tooth ablation);
  `notes/part4_charge_atomic_split.md` (IV-4 = the native-`P` departure); `docs/model_map.md` §3.4 (charge) +
  §4 (honest departure ledger). Reshape directive + review trail:
  `research/pde_ledger_v2/_scratch/stage033_reshape_directive.md`. ⛔ **Not retained** — it lived in gitignored `_scratch/` and no copy survives; this line records that a directive existed, it is not an auditable citation.
