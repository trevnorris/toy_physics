# Independent physics review — S11c-c2 N6 cross-engine comparator (SCRIPT; astra-written)

## Artifact
`research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py` + its test
`research/pde_ledger_v3/scripts/test_S11c_c2_N6_cross_engine_comparator.py` (astra-written).

## What to check
This is a cross-engine MEASUREMENT INSTRUMENT: it joins the two blind N6 engines' committed tag streams by object
name, prints `operand_A` (SymPy), `operand_B` (WL) and the typed three-valued `A − B`, and DECIDES NOTHING. Verify
the WORKING instrument does what the governing physics requires and contains no false-agreement / false-disagreement
path. The governing contract is `directives/S11c_c2_SHARED_PHYSICS.md` §N8/T7 (~477) + §5c (~303) + §1b (~76); the
build directive is `directives/S11c_c2_N6_comparator_build_directive.md` (⚠ read it to know what was ASKED, but ⛔ an
artifact can satisfy its directive and still be wrong — form your OWN view of what the N6 cross-engine test can
establish from the engines' actual outputs FIRST).

## What you are handed
The script + test; the two engines' committed outputs (WL `mathematica/out/S11c_c2_N6_mathematica_audit.out` ⚠ 400 MB,
targeted reads only, NEVER cat; SymPy `scripts/out/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.out`); the engine
sources `scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py` + `mathematica/S11c_c2_N6_mathematica_audit.wl`;
the governing spec + build directive; the reuse base `scripts/S11c_{c1,b,a}_cross_engine_comparator.py`.

## Required method (SCRIPT branch — derive independently; ⛔ code-reading alone has missed real defects)
⭐ **Write your OWN derivation/ablation scripts BEFORE trusting the comparator, and save BOTH each script AND its
literal stdout to named absolute paths — ⛔ without these, your claims are discarded.** A prose "I checked X" is worth
nothing.
- ⛔⛔ **FORM ABLATION IS MANDATORY.** Change the STRUCTURE of a load-bearing object in a /tmp COPY of the comparator
  (flip a sign AND an off-diagonal; collapse two distinct axes into one; repoint a paired NAME to a different object's
  payload) and report the LITERAL diff. A COEFFICIENT rescale tests arithmetic; only a FORM change tests physics.
- **One-sided corruption** for every claimed-independent route: corrupt ONE engine's operand for one family and
  confirm ONLY that family's residual moves (if the other side moves too, the join is not what it claims).
- Probe for: a value verified using the predicate that produced it; a conclusion emitted as an unconditional literal;
  an `assert` that precedes the value it guards (report any assert on a measured payload — the script must PRINT, not
  assert); a per-case VERDICT token (`PASS`/`AGREE`/`COVARIANT`/…) the script decides (BANNED).

## SPECIFIC scrutiny (the directive made these decisions — verify the SCRIPT honors them AND that they are correct)
1. ⛔⛔ **The PIT SEAL.** The two engines use DIFFERENT prime sets + sample points (blindness). Confirm the comparator
   NEVER differences `PROBE_NUMERATORS`/`numerator_denominator` across engines and NEVER re-samples to "align". Ablate:
   force a cross-engine PIT subtraction in a /tmp copy — the real comparator must not have that path.
2. ⭐ **The SCHEMA BRIDGE (the load-bearing design call — SCRUTINIZE HARD).** The directive pre-registers ONLY
   mechanical name/CAS identities (grade `{1,η,σ}↔(η,σ)`, face `SUM↔0`, axis-order, `_NODES` literal decode) and emits
   a `{matched, sympy_only_column, wl_only_slot}` pairing table; it ⛔ FORBIDS pre-registering block-name /
   kernel-family / component-vs-formal-basis equality (a physics identification = the forbidden blanket collapse).
   VERIFY the script does exactly this: does it silently join columns whose keys are NOT mechanically identical (a
   manufactured join / false agreement)? Does an unmatched block/kernel key correctly become a residual/accounting
   row, not a count? Construct a synthetic block-name mismatch and confirm it is NOT auto-joined.
3. **The carrier-rep residual.** WL blind `C_E` vs SymPy imported-slab `C_E` must produce a COMPUTED three-valued
   residual (DO-NOT-FOLD), ⛔ not be skipped as a "seal". Confirm the residual is printed.
4. **One-sided support.** Nonzero-support must be typed `NONZERO_WITNESSED` vs `NO_NONZERO_FOUND` (the latter =
   UNDECIDED, ⛔ not "structurally zero"); ⛔ native support booleans are never subtracted across engines.
5. **`_NODES` reconstruction.** Confirm reconstruction follows `ARITHMETIC_DAG.nodes[ref]` and that a
   `(srepr,"formal_jet")` literal ROUND-TRIPS (⛔ is not treated as an opaque name → a manufactured disagreement); a
   missing `ARITHMETIC_DAG` is `parse_failed`, not a zero residual.
6. **Join-set completeness + no leak.** All shared reconcile + covariance families (incl. `R_COV_CONTROL_DELTA`,
   `SOURCE_CONTROL_DELTA`, `PHI_DOMAIN_CENSUS`, `ACTUAL_CONTROL_PARAMETERS`) are joined; SymPy-only families are
   `sympy_only` accounting (⛔ never a disagreement, ⛔ no manufactured WL operand); NO residual target / expected
   value / prime literal is baked in (rule 5).
7. **DoD teeth.** Confirm zero-extract from two NONEMPTY declared-shared containers trips an OPERATIONAL FAILURE (not
   read as agreement/disagreement); `join>0` is NOT the exit criterion; an all-`deferred_oversize` run fails.

## Physics filter
Report a finding only if it catches a way the physics/measurement could be WRONG (a manufactured join, a missed
channel, a differenced PIT, a leaked value, a false-zero symbolic channel, an assert-before-emit). ⛔ Not style.

## Ablation sandbox + ops
⛔ Copy the comparator to /tmp and ablate the COPY; ⛔ never modify the working tree. ⚠ The comparator runs on a 400 MB
WL stream — keep the lazy materialize/release discipline; if a run OOMs, check `free -h` (a leftover process, not the
comparator). Save every ablation script + its literal stdout to named absolute paths and report those paths. This is
pure Python/SymPy (⛔ no Mathematica kernel — no seat/timeout-kernel rules apply). End with SOUND / NOT-SOUND + ranked
findings, each with the file/line + the ablation stdout that proves it.
