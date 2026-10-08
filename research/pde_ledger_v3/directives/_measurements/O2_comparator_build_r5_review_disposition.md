# O2 comparator build r5: review dispositions (orchestrator)

**Artifact:** the r5 comparator, its tests and its ablation harness, as frozen in
`_scratch/s9b_build/o2_comparator_build_review_baseline_r5.sha256`:

| File | sha |
|---|---|
| comparator | `16a15825…` |
| tests | `325120e9…` |
| harness | `21f70e9a…` |

The fresh Codex author (session `01a118a9…`, third revision) made it to amended directive `44f78617`, from brief
`_scratch/s9b_build/o2_comparator_repair5_prompt.md`. The user chose this repair route on 2026-10-07: the same author,
one fix, with W3 only stated as a limit. Builder run: 83 tests OK; 47 harness ablations, each failing its designated
test; full guarded run exit 0. Sources were last edited before every run. The artifact is Codex-written, so it
gets fresh Claude + Grok (G1), reviewed until clear (G4).

**Legs.** Both used the identical prompt `_scratch/s9b_build/o2_comparator_build_review_prompt_r5.md` (review round
6), and both reported before adjudication.
- **Grok:** CLEAR, no findings (`_scratch/s9b_build/o2_comparator_build_review_r5_grok.txt`; evidence copied to
  `…_r5_grok_evidence/`).
- **Fresh Claude (opus):** CLEAR, no findings, four non-blocking notes (`_scratch/s9b_build/o2_comparator_build_review_r5_claude.md`;
  evidence copied to `…_r5_claude_evidence/`).

**What the legs established independently:**
- **Coverage:** 160 joined and 204 unjoined, with 0 unaccounted. Each leg made its own enumeration first.
- **Joins:** all 160 resolve in both engines and the joins are injective.
- **Algebraic residuals:** Claude recomputed 15 rows in Mathematica, including the transpose layout and a relation;
  Grok recomputed a set of SymPy rows. Every residual is 0.
- **Balances:** the closed parts match on recomputation, including `hold_normal` in both engines.
- **Unformable residuals:** `not_formed` with true reasons, and an equal name is never a zero (a FORM ablation of
  that guard fails its test).
- **Printed limit:** its counts match independent recounts.
- **Controls:** all 30 item-6 bullets are controlled, and each harness knife fails its test.
- **One-sided corruptions:** WL `VR`→`VR^2`, a WL `XiW'`→`XiW''` inside an OPEN occurrence, and a renamed PY profile
  are all printed.
- **No target:** no verdict token, and exit is 0 on a nonzero residual.

Each check of mine below is a mechanical lookup: sed, grep and sha256sum on the artifact, plus verbatim retrieval of
the legs' filed stdout. Commands and literal output are in `O2_comparator_build_r5_review_disposition_lookups.md`.

## Round-5 findings in the r5 artifact

**W1 resolved.** In the production output, `"wl::6"` and `"wl::3,8"` now occur on 0 lines (r4: 6 each).
`"wl::D"` and `"wl::Map"` remain on 14 lines each, down from 20. Retrieval shows where they still sit:
- raw operand serialization (`"Symbol","wl::D"`, 14 lines);
- the raw symbol census `structure()` → `named_operands_and_heads` (comparator L838–859).

The printed `raw_syntax_policy` labels that census as diagnostics, not an OPEN comparison. No
`named_OPEN_operands` inventory contains them, and they appear in no difference pair (0 lines each; lookups W1
detail, helper `_scratch/s9b_build/gen/o2_r5_w1_detail.sh`).

**W2 resolved.** The harness now names `term_orientation` on 4 lines (r4: 0). Both legs ran the harness, and every
knife fails its test.

**W3 resolved.** The printed `not_compared` (L992–993) now names:
- held derivative or variation variables;
- held aggregate sets and binders;
- scalar coefficients beyond their orientation sign;
- provenance and constructor options.

## Notes, adjudicated

| # | Note (leg) | Disposition |
|---|---|---|
| N1 | A negative non-integer SymPy rational coefficient (`Rational(-1,2)`) is not read as a sign flip by `term_orientation` or `signed_actions`. Both flip only for a `Number` that starts with `-` (L793, L1151). On a synthetic fixture this gives a spurious difference when the signs agree and an empty one when they disagree. (Claude N1) | **Note; carried to the record.** Under the user's finite-contract filter, this is a path the contract does not list, and the output is faithful on the measured streams. Two legs measured that independently: Claude reports `terms_with_Rational_on_sign_path 0` on all 8 balance sides, with 0 differing rows under a variant that counts negative rationals; Grok reports `sign_mismatches=0 negative_rational_terms=0`. The record states that the comparator's orientation reading has not been validated for negative rational coefficients, so any reuse on other streams must repair or check this first. |
| N2 | `entry_trace_NATIVE_HOLD_LOAD` pairs a PY trace string with WL `B_HOLD_LIVE/Entries/3/Origin`, the MechanicalFaceSupport entry. It is provenance only, and the strings are identical. (Claude N2) | **Note.** It changes no physics row. The record claims nothing from trace rows. |
| N3 | The control for "a nested sibling removed" (O7) removes a record sibling, not an OPEN-argument sibling. (Claude N3) | **Note.** The bullet has a control that fails, and the leg found the output faithful. |
| N4 | Inventories are unioned across same-role occurrences within a row. Component-indexed roles and per-component balances limit the effect. (Claude N4) | **Note; carried to the record.** This is a form of the stated multiplicity limit. The record does not read an empty difference on a row with repeated same-role occurrences as per-occurrence agreement. |

**Carried from round 5 (note N1 there):** the generalized-rates pairing differs in scope. PY's is normal-only and
WL's is rotational plus normal. PY's rotational work is unjoined, with a true reason.

## A measurement for the record (not a comparator defect)

The comparator prints WL-only `xi_w''(r)` live objects on the in-plane `MomentumCurrent` and `MaterialEnergyCurrent`
balance entries. The Claude leg confirmed independently that this is faithful:
- WL's `OpenFirstVariation` entries carry `XiW` derivative orders `{{1, 43}, {2, 12}}`;
- PY's `OPEN_MomentumFlux_*` occurrences carry only `{1: 105}`.

Its FORM ablation that forgets derivative order makes this difference disappear (135 leaf paths) and fails two
tests. This is a cross-engine difference in OPEN content. Interpreting it belongs to the O2 record (sub-step 7),
under M1: it is preserved, never designed away.

## Verdict

**ACCEPT r5.** Both legs are CLEAR. Nothing outstanding changes what the comparator computes or what may be
claimed from its output, provided the record carries the limits above (N1, N4, the generalized-rates scope) together
with the printed `not_compared` list (G4). This closes O2 sub-step 6.

The history ran six review rounds, two authors and 17 findings. The closed results never moved. The final scope is
the user's finite contract (amendment 1, `44f78617`), and the final repair route was also the user's choice.
