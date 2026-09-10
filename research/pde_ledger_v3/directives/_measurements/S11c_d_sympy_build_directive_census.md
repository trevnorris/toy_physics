# S11c-d SymPy build-directive census — MECHANICAL FACTS

Scope: fact-lookup ONLY (existence, verbatim retrieval, literal-match counting, shape/cardinality of
named stored objects). No adjudication, no algebraic/symbolic residual, no derived predicate.
Gathered 2026-09-10. Working dir: `/var/projects/toy_physics/research/pde_ledger_v3`.

Inputs (sizes as on disk):
- `scripts/S11c_b_exports.py` — 61,608,243 B (base)
- `scripts/S11c_c1_exports.py` — 906,467 B (c1 delta)
- `scripts/S11c_c2_exports.py` — 22,441,522 B (c2 delta)
- `scripts/ledger_fold.py` — loader; public API `load_model(base_path, *delta_paths)` (line 102),
  `check_consumer(fold, import_keys)` (line 199). Each `*_exports.py` exposes a public
  `LEDGER = MappingProxyType({...})` accessor; rows are restored from srepr strings via an internal
  `_restore`/`eval`. `load_model` returns `(merged_dict, audit_dict)`; `audit` is a plain dict with keys
  `overwrites` and `source_row_counts`.

`load_model` over all three parents ran fine (peak RSS ~1.6 GB); no fallback to raw-text was needed for
Part A. The giant closed-row census (Part C) was done by BOTH raw-text literal-match parsing (no eval) AND
sympy `.atoms()` cross-check; the two agree.

---

## PART A — 3-parent fold mechanics

### Command (A.1 + A.2) — `part_a.py`
```python
import sys, os
sys.path.insert(0, "scripts")
from ledger_fold import load_model, check_consumer, AmbiguousSymbolError, ClosureError, ManifestError
fold, audit = load_model("scripts/S11c_b_exports.py",
                         "scripts/S11c_c1_exports.py",
                         "scripts/S11c_c2_exports.py")
print("audit keys =", sorted(audit.keys()))
print("source_row_counts =", audit["source_row_counts"])
print("final fold row count =", len(fold))
print("overwrites =", audit["overwrites"])
# disjointness of exact keys per file via module.LEDGER.keys()
# check_consumer over the two closed rows as roots:
res = check_consumer(fold, ["s11cc2ClosedSlabOperator", "s11cc2ClosedCouplingKernel"])
print("closure size =", len(res["closure"]))
print("symbol_edges =", len(res["symbol_edges"]))
print("dimension_edges =", len(res["dimension_edges"]))
```
### Literal stdout
```
=== PART A.1: load_model 3-parent fold ===
SUCCESS: load_model returned
type(fold) = dict
type(audit) = dict
audit keys = ['overwrites', 'source_row_counts']
source_row_counts = [('scripts/S11c_b_exports.py', 2441), ('scripts/S11c_c1_exports.py', 44), ('scripts/S11c_c2_exports.py', 70)]
final fold row count = 2555
overwrites (count) = 0
overwrites (list)   = []

base key count = 2441
c1 key count   = 44
c2 key count   = 70
c1 INTERSECT c2 (exact keys) = []
base INTERSECT c1 = []
base INTERSECT c2 = []

=== PART A.2: check_consumer over the 3-parent fold ===
check_consumer(closed roots) SUCCESS
import_keys       = ['s11cc2ClosedCouplingKernel', 's11cc2ClosedSlabOperator']
closure size      = 213
symbol_edges (n)  = 273
dimension_edges(n)= 45
```

### Facts (A.1)
- `load_model(base, c1, c2)` **succeeds**. Return = `(dict fold, dict audit)`.
- Base row count = **2441**; c1 delta = **44**; c2 delta = **70**.
- Final fold row count = **2555** (= 2441 + 44 + 70; strictly additive).
- `overwrites` = **[] (0)**. `audit` exposes exactly two fields: `overwrites`, `source_row_counts`.
- Exact-key intersections are **all empty**: base∩c1 = base∩c2 = c1∩c2 = ∅. The three parents write
  disjoint key sets; the fold is a pure union with no exact-key replacement.

### Facts (A.2 — the fold guard)
- Guard functions **exist**: `check_consumer(fold, import_keys)` (line 199); helpers
  `_identity_rows`, `_referenced_atoms`, `_producer_identities`; raises `AmbiguousSymbolError` /
  `ClosureError` / `ManifestError`. Also present: `assert_lookups_equal_manifest`,
  `assert_delta_is_minimal`.
- Called with the two closed rows as roots, `check_consumer` **resolves with no ambiguity**:
  closure size **213**, symbol_edges **273**, dimension_edges **45**. No `AmbiguousSymbolError` raised.

---

## PART B — c2 delta OWN write-keys (raw text)

### Command (B.3) — structural key extraction
```
sed -n '11,80p' scripts/S11c_c2_exports.py | grep -oP "^'[^']+'"
```
`_LEDGER` spans lines 10-81 of `scripts/S11c_c2_exports.py`; each top-level key begins at column 0.

### B.3 — COMPLETE c2 write-key set (70 keys, verbatim)
```
s11cc2ClosedCouplingKernel
s11cc2ClosedSlabOperator
s11cc2Coefficientm1Profile
s11cc2Coefficientm1ProfileDimension
s11cc2Coefficientw1Profile
s11cc2Coefficientw1ProfileDimension
s11cc2FieldeW
s11cc2FieldeWDimension
s11cc2Fieldtheta
s11cc2FieldthetaDimension
s11cc2Fieldu1
s11cc2Fieldu1Dimension
s11cc2Fieldu2
s11cc2Fieldu2Dimension
s11cc2Fieldu3
s11cc2Fieldu3Dimension
s11cc2FourierW1ProfileHatTransfer
s11cc2FourierW1ProfileHatTransferDimension
s11cc2FourierW1ProfileJetHat1
s11cc2FourierW1ProfileJetHat1Dimension
s11cc2FourierW1ProfileJetHat2
s11cc2FourierW1ProfileJetHat2Dimension
s11cc2FourierW1ProfileJetHat3
s11cc2FourierW1ProfileJetHat3Dimension
s11cc2MiddleMomentum1
s11cc2MiddleMomentum1Dimension
s11cc2MiddleMomentum2
s11cc2MiddleMomentum2Dimension
s11cc2MiddleMomentum3
s11cc2MiddleMomentum3Dimension
s11cc2OutgoingNormalMomentum
s11cc2OutgoingNormalMomentumDimension
s11cc2TestA0
s11cc2TestA0Dimension
s11cc2TestA1
s11cc2TestA1Dimension
s11cc2TestA2
s11cc2TestA2Dimension
s11cc2TestE
s11cc2TestEDimension
s11cc2TestPhi
s11cc2TestPhiDimension
s11cc2TestTheta
s11cc2TestThetaDimension
s11cc2Time
s11cc2TimeDimension
s11cc2TrialA0
s11cc2TrialA0Dimension
s11cc2TrialA1
s11cc2TrialA1Dimension
s11cc2TrialA2
s11cc2TrialA2Dimension
s11cc2TrialE
s11cc2TrialEDimension
s11cc2TrialPhi
s11cc2TrialPhiDimension
s11cc2TrialTheta
s11cc2TrialThetaDimension
s11cc2X1
s11cc2X1Dimension
s11cc2X2
s11cc2X2Dimension
s11cc2X3
s11cc2X3Dimension
s11cc2Y1
s11cc2Y1Dimension
s11cc2Y2
s11cc2Y2Dimension
s11cc2Y3
s11cc2Y3Dimension
```
(70 keys total. Note: `s11cc2OutgoingNormalMomentum` + its `Dimension` exist as write-keys, though
the closed rows use `s11cc2MiddleMomentum*` for internal legs — see Part C.7.)

### B.4 — casing confirmation of the candidate keys
| Candidate (as asked) | Status | Exact spelling found | `*Dimension` companion |
|---|---|---|---|
| `s11cc2ClosedSlabOperator` | EXISTS | `s11cc2ClosedSlabOperator` | none (no Dimension row) |
| `s11cc2ClosedCouplingKernel` | EXISTS | `s11cc2ClosedCouplingKernel` | none (no Dimension row) |
| `s11cc2Fieldtheta` | EXISTS | `s11cc2Fieldtheta` (lowercase `theta`) | `s11cc2FieldthetaDimension` |
| `s11cc2FieldeW` | EXISTS | `s11cc2FieldeW` | `s11cc2FieldeWDimension` |
| `s11cc2Fieldu1` | EXISTS | `s11cc2Fieldu1` | `s11cc2Fieldu1Dimension` |
| `s11cc2Fieldu2` | EXISTS | `s11cc2Fieldu2` | `s11cc2Fieldu2Dimension` |
| `s11cc2Fieldu3` | EXISTS | `s11cc2Fieldu3` | `s11cc2Fieldu3Dimension` |
| `s11cc2Coefficientw1Profile` | EXISTS | `s11cc2Coefficientw1Profile` | `s11cc2Coefficientw1ProfileDimension` |
| `s11cc2Coefficientm1Profile` | EXISTS | `s11cc2Coefficientm1Profile` | `s11cc2Coefficientm1ProfileDimension` |
| `s11cc2FourierW1ProfileHatTransfer` | EXISTS | `s11cc2FourierW1ProfileHatTransfer` | `s11cc2FourierW1ProfileHatTransferDimension` |
| `s11cc2FourierW1ProfileJetHat1` | EXISTS | `s11cc2FourierW1ProfileJetHat1` | `s11cc2FourierW1ProfileJetHat1Dimension` |
| `s11cc2FourierW1ProfileJetHat2` | EXISTS | `s11cc2FourierW1ProfileJetHat2` | `s11cc2FourierW1ProfileJetHat2Dimension` |
| `s11cc2FourierW1ProfileJetHat3` | EXISTS | `s11cc2FourierW1ProfileJetHat3` | `s11cc2FourierW1ProfileJetHat3Dimension` |

All 13 candidates EXIST with the exact casing asked. The two Closed* rows carry NO `*Dimension`
companion (they have no `dimension_key`); every carrier/hat/field/coefficient/momentum key does.

### B.5 — per-(anchoring,density) case STRUCTURE of the two Closed rows
Command (first 700 chars of each value srepr) and all-caps `Str('...')` structural token census.
```
=== s11cc2ClosedSlabOperator (head) ===
Tuple(Tuple(Tuple(Str('LAB_HELD'), Str('RHO4_CONSTANT')), Tuple(Tuple(Str('VALUE'),
  Tuple(Tuple(Str('U'), Tuple(Mul(Symbol('epsilon_shape'), Add(...
=== s11cc2ClosedCouplingKernel (head) ===
Tuple(Tuple(Tuple(Str('LAB_HELD'), Str('RHO4_CONSTANT')), Tuple(Tuple(Str('VALUE'),
  Tuple(Tuple(Str('TRANSVERSE_TO_THICKNESS'), Tuple(Tuple(Str('THETA'), Mul(...

ALL-CAPS Str structural tokens (token: count)
s11cc2ClosedSlabOperator:  LAB_HELD:2  MATERIAL_ADVECTED:2  RHO4_CONSTANT:2  RHOBR_CONSTANT:2
   VALUE:4  U:8  THETA:8  E_W:8  COMPUTED_BRANCH_BINDINGS:4  FOURIER_PROFILE_BINDINGS:4
   MULTIGRADE:4  DIMENSION_L_T_M:4
s11cc2ClosedCouplingKernel: LAB_HELD:2  MATERIAL_ADVECTED:2  RHO4_CONSTANT:2  RHOBR_CONSTANT:2
   VALUE:4  THETA:16  E_W:16  DIV_U:16  TRANSVERSE_TO_THICKNESS:8  THICKNESS_TO_TRANSVERSE:8
   COMPUTED_BRANCH_BINDINGS:4  FOURIER_PROFILE_BINDINGS:4  MULTIGRADE:4  DIMENSION_L_T_M:4

row metadata tail (both): 'value_kind':'COMPUTED_OBJECT', 'class':'DERIVED', 'step':'S11c-c2',
   'route':'F9A_ABSENT'
```
Structure (as observed — NOT a Python dict; a nested `Tuple`-of-pairs association list):
- Each Closed row value = `Tuple` of **4 case-entries**. Each entry =
  `Tuple( Tuple(Str(α), Str(ρ)), <payload> )`.
- Case key = `Tuple(Str(α), Str(ρ))` with **α ∈ {`LAB_HELD`, `MATERIAL_ADVECTED`}**,
  **ρ ∈ {`RHO4_CONSTANT`, `RHOBR_CONSTANT`}**. Each token appears exactly twice ⇒ the full 2×2 = 4-case
  grid, keyed exactly as those verbatim strings.
- Payload per case is keyed under `Str('VALUE')` (`VALUE`:4 = once per case), holding an inner
  association list:
  - Slab operator components: `Str('U')`, `Str('THETA')`, `Str('E_W')` (8 each = 2 per case ×4).
  - Coupling kernel: an outer sector layer `Str('TRANSVERSE_TO_THICKNESS')` /
    `Str('THICKNESS_TO_TRANSVERSE')` (8 each), then components `Str('THETA')`, `Str('E_W')`,
    `Str('DIV_U')` (16 each).
  - Both rows also carry per-case metadata sub-blocks keyed `COMPUTED_BRANCH_BINDINGS`,
    `FOURIER_PROFILE_BINDINGS`, `MULTIGRADE`, `DIMENSION_L_T_M` (4 each).
- Row metadata: `class='DERIVED'`, `step='S11c-c2'`, `value_kind='COMPUTED_OBJECT'`,
  `route='F9A_ABSENT'`; no `dimension_key`.

---

## PART C — 3-D Fourier-element census in the two closed rows

Two independent methods, agreeing:
(a) raw-text literal-match parsing of the stored value srepr (no eval), and
(b) sympy `.atoms(...)` on the restored row (peak RSS ~1.6 GB).

### C.6 — Fourier-of-profile HAT symbols: PRESENT/ABSENT + occurrence counts + literal arguments

Literal `Function('NAME')(` occurrence counts in each row:

| Hat name | in ClosedSlabOperator | in ClosedCouplingKernel |
|---|---|---|
| `s11cc2FourierW1ProfileHatTransfer` | PRESENT, 88 | PRESENT, 550 |
| `s11cc2FourierW1ProfileJetHat1` | PRESENT, 100 | PRESENT, 136 |
| `s11cc2FourierW1ProfileJetHat2` | PRESENT, 100 | PRESENT, 136 |
| `s11cc2FourierW1ProfileJetHat3` | PRESENT, 100 | PRESENT, 136 |
| `s11cc1_w1_profile_hat_transfer` (c1-level) | ABSENT, 0 | ABSENT, 0 |
| `s11cc1_w1_profile_jet_hat_1/2/3` (c1-level) | ABSENT, 0 | ABSENT, 0 |

sympy cross-check — distinct c2 hat functions applied in each row:
`['s11cc2FourierW1ProfileHatTransfer','s11cc2FourierW1ProfileJetHat1','...JetHat2','...JetHat3']`
in BOTH rows. The c1-level snake_case hats are NOT applied functions in these rows (0).

**Literal arguments** (each hat is a **3-argument** applied function — one arg per spatial component;
verbatim srepr per distinct argument-list shape, with occurrence counts). Two distinct argument shapes
occur, identical in both rows:

Shape 1 — "transfer argument" `(k_out − k_in)` per component (the dominant shape):
```
Add(Mul(Integer(-1), Symbol('s11cc1_k_input_1', real=True)), Symbol('s11cc1_k_output_1', real=True)),
Add(Mul(Integer(-1), Symbol('s11cc1_k_input_2', real=True)), Symbol('s11cc1_k_output_2', real=True)),
Add(Mul(Integer(-1), Symbol('s11cc1_k_input_3', real=True)), Symbol('s11cc1_k_output_3', real=True))
```
(HatTransfer: x502 coupling / x40 slab; each JetHat: x88 coupling / x40–x52 slab.)

Shape 2 — "middle-leg / intermediate-momentum argument", two sub-forms, each ×24 in every hat & row:
```
# (k_out − k_mid)
Add(Symbol('s11cc1_k_output_1', real=True), Mul(Integer(-1), Symbol('s11cc2MiddleMomentum1', real=True))), ... (_2, _3)
# (k_mid − k_in)   [srepr: MiddleMomentum − k_input]
Add(Mul(Integer(-1), Symbol('s11cc1_k_input_1', real=True)), Symbol('s11cc2MiddleMomentum1', real=True)), ... (_2, _3)
```
So every hat appears BOTH at the transfer difference `(k_out−k_in)` AND at the two middle-leg
differences `(k_out−k_mid)` and `(k_mid−k_in)`, where `k_mid = s11cc2MiddleMomentum{1,2,3}`.
(Physics labels not asserted — distinction is by the verbatim argument text only.)

### C.7 — momentum symbols: PRESENT/ABSENT + literal name-mention counts

| Momentum name | ClosedSlabOperator | ClosedCouplingKernel |
|---|---|---|
| `s11cc1_k_output_1` | 1188 | 3820 |
| `s11cc1_k_output_2` | 1188 | 3856 |
| `s11cc1_k_output_3` | 1188 | 3820 |
| `s11cc1_k_input_1` | 1044 | 2904 |
| `s11cc1_k_input_2` | 1044 | 2904 |
| `s11cc1_k_input_3` | 1044 | 2904 |
| `s11cc1_q_out_input` | ABSENT, 0 | ABSENT, 0 |
| `s11cc1_q_out_output` | ABSENT, 0 | ABSENT, 0 |
| `s11cc2MiddleMomentum1` | 340 | 484 |
| `s11cc2MiddleMomentum2` | 340 | 484 |
| `s11cc2MiddleMomentum3` | 340 | 484 |
| `s11cc2OutgoingNormalMomentum` | PRESENT, 952 | PRESENT, 3422 |

Middle-momentum symbol names found (exact): `s11cc2MiddleMomentum1`, `s11cc2MiddleMomentum2`,
`s11cc2MiddleMomentum3`. `s11cc1_q_out_input` / `s11cc1_q_out_output` do NOT appear in either closed row.

### C.8 — Integral (3-D measures) in the closed rows: PRESENT

| Metric | ClosedSlabOperator | ClosedCouplingKernel |
|---|---|---|
| literal `Integral(` count | 218 | 526 |
| sympy `.atoms(Integral)` distinct | 26 | 90 |

(The two differ because `.atoms` returns a set of structurally-distinct Integral subexpressions;
the literal count counts every textual occurrence. Both are non-zero ⇒ Integrals are present.)

Distinct integration variables (verbatim), identical set in BOTH rows:
```
s11cc2Y1, s11cc2Y2, s11cc2Y3            # real-space d^3y  (limits -oo, oo)
s11cc1_k_output_1, _2, _3               # momentum d^3k_out (limits -oo, oo)
s11cc1_k_input_1, _2, _3                # momentum d^3k_in  (limits -oo, oo)
s11cc2MiddleMomentum1, _2, _3           # momentum d^3k_mid (limits -oo, oo)
```
Every limit spec is `Tuple(Symbol('<var>', real=True), -oo, oo)`. Per-variable limit-spec counts:
slab — Y{1,2,3}:218, k_output:202, k_input:84, MiddleMomentum:12; coupling — Y{1,2,3}:526,
k_output:510, k_input:168, MiddleMomentum:12. So momentum convolutions ARE represented as explicit
`Integral` nodes over both position (d³y) and momenta (d³k_out, d³k_in, d³k_mid) — NOT as bare
hats-at-a-momentum-difference and NOT via DiracDelta (see C.9).

### C.9 — DiracDelta (3-D momentum deltas) in the closed rows: ABSENT

| Metric | ClosedSlabOperator | ClosedCouplingKernel |
|---|---|---|
| literal `DiracDelta(` count | 0 | 0 |
| sympy `.atoms(DiracDelta)` | 0 | 0 |

Both methods agree: **zero** DiracDelta in either closed row. The directive's claim that c2 stripped the
3-D deltas from the closed rows is verified literally.

### C.10 — contrast: `dtn_kernel` (a c1 import row) DOES carry the 3-D deltas
`dtn_kernel` is a c1 write-key (`scripts/S11c_c1_exports.py`, line 85); value type `Tuple`.
Command: `c1.LEDGER["dtn_kernel"]["value"]` then `.atoms(...)`.
```
dtn_kernel atoms(DiracDelta) count: 3
   DiracDelta(Add(Mul(Integer(-1), Symbol('s11cc1_k_input_2')), Symbol('s11cc1_k_output_2')))
   DiracDelta(Add(Mul(Integer(-1), Symbol('s11cc1_k_input_1')), Symbol('s11cc1_k_output_1')))
   DiracDelta(Add(Mul(Integer(-1), Symbol('s11cc1_k_input_3')), Symbol('s11cc1_k_output_3')))
dtn_kernel atoms(Integral) count: 0
dtn_kernel momentum-like free symbols: ['q_out','s11cc1_k_input_1','_2','_3',
   's11cc1_k_output_1','_2','_3','s11cc1_q_out_input','s11cc1_q_out_output']
```
So:
- `dtn_kernel` has **3 DiracDelta** nodes, one per component, each at the transfer difference
  `(k_output_i − k_input_i)`, and **0 Integral** nodes.
- Its momentum legs: `s11cc1_k_output_1/2/3`, `s11cc1_k_input_1/2/3`, `q_out`,
  `s11cc1_q_out_input`, `s11cc1_q_out_output`. (It also carries the c1 hats as bare SYMBOLS —
  `s11cc1_w1_profile_hat_transfer`, `s11cc1_w1_profile_jet_hat_1/2/3` — not as applied functions.)
- Confirms: the 3-D momentum deltas live on `dtn_kernel` (c1 import), NOT on the c2 closed rows.

---

## Compact summary table

| # | Fact | Result |
|---|---|---|
| A.1 | load_model(base,c1,c2) | SUCCESS; rows 2441 + 44 + 70 = **2555**; overwrites=0; keys disjoint (all pairwise ∩ = ∅); audit fields = {overwrites, source_row_counts} |
| A.2 | check_consumer over 2 closed roots | resolves, NO ambiguity; closure=213, symbol_edges=273, dimension_edges=45 |
| B.3 | c2 own write-keys | **70** keys (full list above) |
| B.4 | 13 candidate keys | ALL EXIST, casing as asked; `Fieldtheta` lowercase `theta`; Closed* rows have NO Dimension companion; all carrier/hat/field/coeff/momentum keys have a `*Dimension` |
| B.5 | case structure | Tuple-of-pairs assoc-list; case key = `Tuple(Str(α),Str(ρ))`, α∈{LAB_HELD,MATERIAL_ADVECTED}, ρ∈{RHO4_CONSTANT,RHOBR_CONSTANT}; 4 cases; payload under `Str('VALUE')`; slab comps U/THETA/E_W; coupling sectors TRANSVERSE_TO_THICKNESS/THICKNESS_TO_TRANSVERSE × comps THETA/E_W/DIV_U |
| C.6 | c2 hats | all 4 (HatTransfer, JetHat1/2/3) PRESENT in both rows; c1 snake_case hats ABSENT; each hat = 3-arg fn; args = transfer `(k_out−k_in)` AND middle-leg `(k_out−k_mid)`,`(k_mid−k_in)` with k_mid=MiddleMomentum |
| C.7 | momenta | k_output_{1,2,3}, k_input_{1,2,3}, MiddleMomentum_{1,2,3}, OutgoingNormalMomentum PRESENT; q_out_input, q_out_output ABSENT |
| C.8 | Integrals | PRESENT (slab 218/26, coupling 526/90 literal/atoms); vars = Y1/Y2/Y3 (d³y) + k_out/k_in/MiddleMomentum (d³k), all limits −oo..oo |
| C.9 | DiracDelta | **0 in both** closed rows (literal + atoms) |
| C.10 | dtn_kernel (c1) | 3 DiracDelta at (k_out−k_in), 0 Integral; legs k_out/k_in + q_out/q_out_input/q_out_output; deltas live here, not on closed rows |
