# stage012 — PREDICTION, written BEFORE the rewrite and before any comparison

Written 2026-07-26 at HEAD `6f5c4424`, from the *committed* `.out` and the *pre-rewrite* `.py`.
Purpose: a mismatch predicted beforehand is evidence; one explained afterward is a rationalisation.

Stage: `ledger_stage012_dtn_pole_ladder_robin`. Basis **`(L,M,T)`**, declared at `.py:193`.

## 1. What the committed `.out` already renders (`out:109-113`, bare compound tuples)

| `.out` label | value |
|---|---|
| `[energy]` | `{2, 1, -2}` |
| `[four-volume]` | `{4, 0, 0}` |
| `[P]` | `{-2, 1, -2}` |
| `[rho]` | `{-4, 0, 0}` |
| `[K]=[P]-5[rho]` | `{18, 1, -2}` |
| `[alpha]` | `{-1, 0, 0}` |
| `[c_S^2]` | `{2, 0, -2}` |
| `[c_S]` | `{1, 0, -1}` |
| `[k]` | `{-1, 0, 0}` |
| `[tan_argument=k*L0]` | `{0, 0, 0}` |
| `[Z00_prefactor=-k]` | `{-1, 0, 0}` |
| `[Z00]` | `{-1, 0, 0}` |
| corrupt `[K]+(1,0,0)` | `{19, 1, -2}` |

## 2. PREDICTIONS

**P1 — every rendered value above is UNCHANGED by the rewrite.** The rewrite is behaviour-preserving;
a changed dimension value is a regression, not a refactor. Any movement here stops the stage.

**P2 — `.py` stdout stays byte-identical** (live stdout vs `tail -n +7` of the committed transcript
diffs to exactly 6 wrapper lines), as it did for 004 and 011 across a 114-line rewrite.

**P3 — `omega_dim` and `mass_dim` will be UNREACHABLE in the `.wl`, and will need waiving.**
Reason: at `.py:556-575` they are *inputs* to `base_dims`, not outputs of `walk()`, and the `.out`
renders neither. This is the **same shape** as stage011, where `OmegaDim`/`MassDim` were exactly the
two unreachable quantities for exactly this reason. If they turn out reachable, that is a *better*
outcome than predicted — note it.

**P4 — ~5 cosmetic name mismatches need a join mapping**, `.py` symbol → `.out` label:
`energy_dim`→`[energy]`, `four_volume_dim`→`[four-volume]`, `pressure_dim`→`[P]`,
`cs_squared_dim`→`[c_S^2]`, `cs_dim`→`[c_S]`, `z00_prefactor_dim`→`[Z00_prefactor]`.
Names are join keys, not values — aligning them is not a fidelity violation. Aligning *values* is.

**P5 — a fractional exponent (3/2) appears in the corrupt walk.** The module stores exact
`sympy.Rational`, so this should pass; if a float appears anywhere, that is a defect.

**P6 — 8 of 12 compared quantities are blind to at least one axis transposition** (per the coverage
census), incl. the fully dimensionless `[tan_argument]={0,0,0}`. So the transposition ablation must be
checked for *how many* quantities catch it, not merely that it fails. Expect roughly 4 detectors.

## 3. ⭐ THE CROSS-STAGE CHECK — the first real one available

stage011's committed `DIM|` lines and stage012's `.out` should agree on **six shared quantities**,
both in basis `(L,M,T)`:

| quantity | stage011 (committed `.out`) | stage012 (committed `.out`) | predicted |
|---|---|---|---|
| energy | `{2, 1, -2}` | `{2, 1, -2}` | **AGREE** |
| four-volume | `{4, 0, 0}` | `{4, 0, 0}` | **AGREE** |
| pressure | `{-2, 1, -2}` | `{-2, 1, -2}` | **AGREE** |
| rho | `{-4, 0, 0}` | `{-4, 0, 0}` | **AGREE** |
| K | `{18, 1, -2}` | `{18, 1, -2}` | **AGREE** |
| c_S² | `{2, 0, -2}` | `{2, 0, -2}` | **AGREE** |

**Prediction: all six agree.** This is a genuine cross-stage consistency result, not a tautology —
the two stages compute them by different routes in different scripts, and nothing has ever compared
them before.

⭐ Independent corroboration against `docs/model_map.md` §2, which was written from a separate
six-agent fan-out and not from these scripts:
- `[K] = M L¹⁸T⁻²` ⇒ `(L,M,T) = {18, 1, -2}` ✓ matches both stages
- `ρ0 = [L⁻⁴]` ⇒ `{-4, 0, 0}` ✓ matches both stages

If any of the six disagree, that is a **real cross-stage conflict** and the stage stops for
adjudication — it is a model question, not a code question.

---

# OUTCOME — `.wl` emission step, recorded after the fact

**P1 — HELD.** All 13 previously-rendered values are unchanged. 17 records emitted (more than the ~13
predicted): the corrupt walk exposed four additional quantities.

**P5 — HELD, and this is the exact-rational check.** The corrupt walk carries
`corrupt_walk.cs_dim={3/2, 0, -1}`, `corrupt_walk.k_dim={-3/2, 0, 0}`,
`corrupt_walk.tan_argument_dim={-1/2, 0, 0}`, `corrupt_walk.z00_prefactor_dim={-3/2, 0, 0}`.
Exact rationals throughout — **no floats**.

**⭐ §3 CROSS-STAGE CHECK — ALL SIX AGREE.** Predicted in writing before the comparison existed:

| quantity | stage011 | stage012 | |
|---|---|---|---|
| energy | `{2, 1, -2}` | `{2, 1, -2}` | AGREE |
| four-volume | `{4, 0, 0}` | `{4, 0, 0}` | AGREE |
| pressure | `{-2, 1, -2}` | `{-2, 1, -2}` | AGREE |
| rho | `{-4, 0, 0}` | `{-4, 0, 0}` | AGREE |
| K | `{18, 1, -2}` | `{18, 1, -2}` | AGREE |
| c_S² | `{2, 0, -2}` | `{2, 0, -2}` | AGREE |

Two stages, two scripts, independent routes, never compared before this. Both corroborate
`docs/model_map.md` §2 (`[K] = M L¹⁸T⁻²`, `ρ0 = [L⁻⁴]`), which was written from a separate fan-out.
**This is the first positive cross-stage consistency result of the rewrite** — the prior findings were
all conflicts or gaps.

**One UNREACHABLE, and it is physically meaningful, not a gap.** `corrupt_walk.z00_dim`: the corrupt
path assigns the sentinel `"dimensionful_tan_argument"` at `.wl:268` because the corrupted tangent
argument genuinely *is* dimensionful (`{-1/2, 0, 0} ≠ 0`). That is the corrupt probe working as
designed. Turning it into an exponent vector would require new computation, which print-only forbids.

**P3 — STILL OPEN.** `omega_dim` and `mass_dim` were neither emitted nor declared unreachable. They
will surface as `py_only` when the `.py` is rewritten; resolve them then.

**P4/P6 — pending the `.py` rewrite.** Note the `.wl` adopted the `.py`'s names directly (including
dotted `clean_walk.cs_dim` forms), so P4's cosmetic mismatches should not materialise.

⚠ **Disclosure to weigh at the `.py` step:** Codex reported that while harvesting names, an `rg`
diagnostic displayed one Python line containing the corrupt-delta constructor, and states it was not
used to choose any value or axis. Independent corroboration that the values are live, not copied:
(a) every added line is `ToString[InputForm[<lookup>]]` with **no literal exponent anywhere** in the
diff; (b) clean and corrupt walks share one Print pattern but yield different values
(`{1,0,-1}` vs `{3/2,0,-1}`), which a hardcode would require two separate literals to fake;
(c) PASS tally unchanged 91 → 91.

---

## 4. What would falsify the approach

- A rendered value changes ⇒ the rewrite is not behaviour-preserving (P1).
- stdout moves beyond the 6 wrapper lines ⇒ the emission leaked into the transcript (P2).
- The transposition ablation is caught by 0 quantities ⇒ stage012 has no working gate, like stage021.
- Any exponent arrives as a float rather than an exact rational (P5).
