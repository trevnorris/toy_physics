# S9b flow-only oracle comparison: the river and acoustic-flow check

**Status:** this record interprets the accepted check (`8d91f448`). The orchestrator wrote it; its own review
is pending. It changes nothing in the S9b record (`steps/S9b_light_bending.md`, accepted `ff3bac63`), which
applies no oracle.

**Evidence.**
- **Check:** Codex wrote `scripts/S9b_flow_oracle_check.py`, its Wolfram script and the data in
  `scripts/S9b_flow_oracle_check_data/`, and the check was reviewed under G1. Disposition:
  `directives/_measurements/S9b_flow_oracle_check_review_r2_disposition.md`.
- **Transcript:** `scripts/out/S9b_flow_oracle_check.out`, sha256 `e66ff635…`.
- **Lookups:** two files in `directives/_measurements/`. "O:" cites a section heading of
  `S9b_flow_oracle_check_outcome_lookups.md`; "C-n" cites a section of `light_sector_closeout_lookups.md`.

## Why the check exists

S9b's record review left a question open: how does S9b's flow-only condition compare with the published
river-model and acoustic-flow profiles? The user chose "Queue it after record". The check would be new
computation, outside the record and after its acceptance, written by Codex and reviewed under G1
(`directives/_measurements/S9b_record_r2_review_disposition.md`, lines 19–21; C1). This record states what
the check printed.

## What was compared

S9b's kept flow-only restriction (`R_FLOW_ONLY` Part B, deflection and radar, both engines) was compared with
equations quoted from two papers. The check's quote file stores each equation with its citation and the
conditions the paper attaches (C4):
- **Hamilton–Lisle (2008), the river model:**
  - (1), the line element;
  - (2), `β = (2GM/r)^(1/2)`;
  - (4), `β(r) = [2M(r)/r]^(1/2)`.
- **Visser (1998), acoustic black holes:**
  - (4), the acoustic metric;
  - (5), continuity;
  - (6)–(7), the force balance;
  - the §7 background, `v = √(2GM/r)` with `ρ ∝ r^(−3/2)`, and (53);
  - §8, `v = c r₀²/r²` (54), and (55);
  - §10, (85)–(87).

Each oracle's own mass and speed parameters stay separate from S9b's orbital `GM` and `c₀`. A substitution
that sets them equal is printed separately and labelled conditional. The check runs three comparisons: the
flow profile, the exchange of brane material, and the optical cone.

## What it printed

**1. Flow profile.**
- **S9b's velocity.** In every case the check reads, S9b's flow-only velocity is `V = ±√2 √GM/√r`, or `0` in
  the zero-velocity cases (O: "S9b's flow-only velocities").
- **Against the oracles with their own parameters,** the squared-profile zero test reduces to an
  identification of parameters, in both engines (O: "zero tests of the squared-profile residuals"):
  - against Hamilton–Lisle (2), to `GM == gmH`;
  - against Hamilton–Lisle (4), to `GM == massH[r]`, which holds only if the oracle's mass function is the
    constant `GM`;
  - against Visser §7, to `GM == gmV`.

  Variants that carry a clock or speed factor print the matching identification with that factor.
- **Under the labelled common-parameter substitution,** both operands are `2GM/r` and the test prints `True`.
  Each engine's case-1 deflection domain has 10 such tests, from Hamilton–Lisle (2), Hamilton–Lisle (4) and
  Visser §7 (C2).
- **On the general branch,** S9b's velocity integration constant is live. There the 10 tests without oracle
  parameters reduce to `constant == 0`, in each general domain (O). S9b's kept branch sets that constant to zero by its supplied
  velocity counting (R-S21-01).
- **Against Visser §8,** the test reduces to `2*py`GM*py`s9b$r^3 == radius0O^4*soundO^2`, which reads
  `2GM r³ = r₀⁴ c²`. That equation holds at one radius, so it does not identify two profiles (O: "against
  Visser section 8").

**2. Exchange of brane material.**
- **Visser's continuity reproduces S9b's exchange.** Evaluated on S9b's flow-only fields, Visser's continuity
  (5) gives S9b's own implied exchange `j_n`. The balance-operator zero test prints `True` in all 18
  comparisons (O: "Comparison 2").
- **The zero-exchange density.** In cases 1–2 of each family, `j_n` vanishes only if `3ρ_br + 2rρ_br′ = 0`;
  in the remaining cases, whose velocity is zero, the test is `True`. Both engines solve that equation as `ρ_br = C₁/r^(3/2)`
  (O: `DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT`). Visser's §7 background quotes the same law,
  `ρ ∝ r^(−3/2)`. The check imposes no oracle density law: the power law comes out as the zero-exchange
  solution for S9b's live density.
- **Hamilton–Lisle** states no continuity operand. That comparison prints `UNMATCHED`, 18 times (O).
- **Light speed on this background.** WL's attached domain prints `c0 == Sqrt[muPerp[r]/rhoBr[r]]` (O). On a
  zero-exchange background the density varies with radius. So whether light keeps the speed `c₀` there
  depends on how the brane's shear response varies with its density. That is the object of R-S1-02, and
  nothing here answers it.

**3. Optical cone.**
- **The cone is a supplied input.** The optics compared are supplied, not an S9b result: the transcript prints
  `"SUPPLIED leading-eikonal input; not an S9b result"` (O).
- **The match with Visser (4).** Of 256 optical residuals, three print `0`, all against Visser (4), under
  conditional speed and binding identifications (C3).
- **The other residuals.** Every optical residual is computed on the general flow-only domains, where S9b's
  velocity is a live function (O: `WORKING_ON` counts). The other 253 zero tests print equations, none
  evaluated at S9b's kept velocity (C3):
  - 240 contain the live velocity, for example the test against Hamilton–Lisle (2) at transcript lines
    9810–9820 (O);
  - 13 contain no velocity symbol: they are conditions on a speed identification or a metric coefficient.
- **No FORM control.** The accepted version has no FORM control on the optical comparison (C5).

**Controls** (C5; check controls record):
- **Visser FORM control:** it moved 8 of the 18 balance values, namely cases 1 and 2 of each engine's
  deflection and radar families. Each moved delta vanishes only where `ρ_br′ = 0`. The other 10 deltas print
  `0`.
- **Hamilton–Lisle FORM control:** dropped by the user's choice, after it did not complete within 50 GiB
  (`run_023`, `run_025`).
- **Review legs' ablations:** the round-2 legs ablated the zero test three ways (an operand shift, removal of
  an assumption, and a FORM change), and each moved the tests it affects (review disposition).

## Conditions with no S9b counterpart

The identification notes print 18 rows as `UNMATCHED` (C5):
- frame and time;
- speed units;
- mass calibration, for each paper;
- Visser's background continuity, which has no source;
- Visser's force balance (6)–(7), which includes an "externally imposed driving potential" `Φ`;
- a barotropic equation of state; an inviscid flow; an irrotational scalar potential;
- §7's and §8's equation of state and density;
- scalar sound versus transverse light;
- the oracles have no shear-bearing sheet and no off-sheet exchange;
- two rows on how each engine's dispersion is bound.

In addition, S9b prints no line element, and Hamilton–Lisle has no continuity operand.

## What this supports and what it does not

**It supports a conditional correspondence.** Suppose S9b's orbital `GM` is identified with each oracle's
mass parameter; the check does not supply that identification. Then S9b's kept flow-only profile is the
river-model profile and Visser's §7 profile. With zero exchange, S9b's live brane density then follows the
`r^(−3/2)` law of Visser §7.

**It does not support:**
- that the brane is the river model's medium or Visser's acoustic fluid (the equation-of-state,
  barotropic, inviscid and irrotational rows are `UNMATCHED`);
- any result of general relativity;
- any statement about S9b's optics, which are a supplied input here;
- that either paper's conditions hold at this model point.

Whether each oracle applies here remains open.

## Second route: what the analog leaves out, and who owns it here

The register's second route asks why a prior result was limited, and whether that objection is answered
here. The check names the following objections. It answers none of them.

1. **What drives the flow.** Visser's background obeys the force balance (6)–(7), which includes an
   externally imposed driving potential `Φ`. The check computes no force comparison
   (`force_accounting_Visser`). So it is not known whether this model's drain sustains the flow without an
   imposed potential. The drain's source and boundary data belong to S12 (R-S12-02), and the force balance
   to O2's carried handoffs.
2. **Where the material goes.** Visser's continuity has no source. In that setting the §7 velocity forces the
   `r^(−3/2)` density. Here the brane can exchange material with the bulk through `j_n`, which the oracles
   lack (`no_off_sheet_exchange`). So this model can carry the same flow in two ways: with that density and
   zero exchange, or with another density and nonzero exchange. The check selects neither. Which one the
   drain produces is again S12's object (R-S12-02), and whether the knit can meet it belongs to R-S21-01.
3. **Scalar sound versus transverse light.** The analog carries a scalar wave, while light here is transverse
   with two polarizations. A shared cone does not identify the two wave equations
   (`scalar_acoustics_vs_transverse_optics`). The polarization and direction exclusions stay as S9b scopes
   them; no owner is added.
4. **The kinematics of analog gravity are not its dynamics.** The survey's limit on analog gravity applies
   unchanged: a matching flow profile does not derive `GM` or a law of gravity. It stays where the S9b
   register pass lists it.

These questions are carried as handoffs, not as new requirement entries. The register's light-closing pass
adds this comparison to R-S21-01 as evidence.
