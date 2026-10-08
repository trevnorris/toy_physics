# S9b: repair and linked-brane decision list (orchestrator)

**Author:** Claude (orchestrator), 2026-10-08. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
both NEEDS REVISION; dispositions in `directives/_measurements/S9b_repair_dl_review_disposition.md`).

This list sets what the next S9b spec version (v10) must contain, and routes the build findings still owed. It names
objects and sources. It states no expected value, sign or relation between symbols.

## D1. Base
- **v10 starts from v8** (`directives/S9b_SHARED_PHYSICS.md` at `c2f1cf2b`, cleared for build). All of v8's
  identifications, references, order counting, Parts A–C and limits stand unless an item below changes them.
- **v9's Part D is not carried** (`05b1a5d5`, preserved, not accepted). Its review found two things:
  - with no steady brane balance, its conditions reduced to Parts B/C;
  - its embedding source was the charge sector.
- **The steady balance now exists** as O2, accepted at `72866fcf` (`steps/O2_steady_brane_balance.md`, §8 handoff).

## D2. Adopted premises
Each premise below is user-selected and is an **adopted substrate input to a conditional model**. None is derived.
Every result that depends on one is flagged with it.

**From O2 (user-selected 2026-10-06):** P1–P4, `directives/O2_premise_decision_list.md` (`77d2c39a`).

**New for S9b:**

| Premise | Content | Qualification carried with it |
|---|---|---|
| **P5** (2026-10-07) | The brane's momentum density is `ρ_br V`. | If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |
| **P6** (2026-10-08) | In steady flow the brane's in-plane stress is an isotropic pressure `p_br(ρ_br)` that depends on `ρ_br` only, because shear relaxes under steady load (P1). `p_br` stays a general function. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. In the optical regime, light still sees the elastic transverse stiffness `μ_⊥`. | Scoped to Part D. It supplies the steady in-plane part of O2's stress `𝒯_br^live` only. |
| **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥ ∝ ρ_br^α` and `α` one live symbol. This is the user's selected live exponent. | Scoped to Part D. There it supplies `μ_⊥` as a function of `ρ_br` alone. Its other O1 `ℳ_⊥` dependences are set aside by this premise, and results are flagged. Parts A–C keep `μ_⊥(x)` and `ρ_br(x)` independent, as in v8. |
| **w-parity** (2026-10-07) | For an electrically neutral mass, the far field is symmetric under `w → −w`. So `ξ_w` and the material `w`-velocity vanish there. | This is a **labelled neutral-sector restriction**. Parts A–C keep `ξ_w` live as in v8, which covers the charged case. Part D is stated for the neutral sector, and prints the label. |

## D3. Part D: the linked brane (new)
**The object.** Part D is the brane's far-field steady in-plane momentum balance. Its ingredients:
- O2's conditional hold balance, from the record's §8 handoff and its sources;
- P3, P5, P6 and w-parity, as supplied;
- v8's mass balance `∇·(ρ_br V) = −j_n`;
- the density link.

**The deliverable.** Part D re-expresses the Part B conditions with that balance in force. It prints:
- which of `δ`, `V`, `ρ_br` and `j_n` it determines relative to `GM`;
- which of them stay free;
- for each condition, the `j_n` it implies.

**What must be true:**
- **Supplied pieces.** The spec supplies O2's balance pieces and the premises as equations, from the O2 sources. It
  does not compose them in advance; the engines compose them.
- **OPEN operands.** Every O2 OPEN operand that P3–P6 and w-parity do not supply stays a general live unknown, not
  an engine-chosen family (v9 review). This includes its admissible gradient and history dependence (O2 record §8).
- **Live symbols.** `j_n`, `p_br(ρ_br)` (and with it `c_comp`), `α` and the asymptotic `ρ_br⁰` stay live. No order or
  value is assigned to `c₀/c_comp`.
- **Bulk-density route.** v8's Part C is unchanged.

## D4. Induced metric and order
- **v8's sentence.** v8 says the induced-metric mass balance differs from the flat form by "a relative `O(ε)`
  correction to the implied `j_n`". That holds only if `∂_r[(∂ξ_w)²] = O(ε/r)` (v9 review, Grok).
- **One rule for v10.** The supplied mass balance is on the coordinate `d³x` measure (O2 record §8). Any reading on
  the induced measure carries the scale `∂_r[(∂ξ_w)²]` live, or states the condition above. No order is attached
  without one or the other.
- **Neutral Part D.** With `ξ_w = 0` (w-parity), the supplied `g_ij` reduces to `δ_ij`, and the engines print that
  reduction. This limits the claim; it adds no term.

## D5. The O2 record's returned obligation
The O2 record returned one interpretation to the orchestrator: the WL-only `ξ_w''` dependence, together with the
seven WL-only first-derivative keys at the same locations.

**Where they sit.** The record places all eight in the momentum-flux and energy-flux roles (its §4 location table).

**Routing:**
- **The obligation stays the orchestrator's** (M1; record §8). S9b does not discharge it. The S9b record carries it
  as undischarged. It is distinct from any step's ownership of a physical input.
- **Retained keys.** Part D retains all eight keys, at every §4 location, for every piece of O2 content that no
  adopted premise supplies (record §8).
- **Supplied content.** Where P5, P6 or the density link supplies an object, the premise states that object's
  dependence, and every result that uses it is flagged.
  - The spec author determines, from the O2 sources, which content each premise supplies.
  - The author does not compose the balance (D3).
  - The user chose P6 over keeping the stress OPEN (2026-10-08).
- **Register.** A register entry follows only if a sourced requirement is established. An unresolved difference
  alone does not create one (record's rest-on criterion).

## D6. Owed build findings: routed to the repair build directive
These come from the build review preserved at `bb94b885`, as listed in that commit's message. They concern implementation, not spec physics. They go to the
repair build directive, which gets its own G2 pass, and the build legs verify them (E2):
- **B1.** The exports bind `w`, `q` and `L` onto unrelated upstream rows (bulk normal coordinate, wave-norm
  coordinate, half-interval size) by bare-symbol equality (`F9B_EQUAL`). The rows are retagged as corroborated S9b
  KNOBs. Both the binding and the corroboration/status promotion are unsupported.
- **B2.** The every-`b` requirement is restated, not reduced. No condition is solved relative to `GM` per stratum,
  and no implied `j_n` is printed (with the `V` sign symbolic).
- **B3.** The radar `ln(1/b²)` coefficient comes from a hand rule, not from Part A's round trip.
- **B4.** Path independence is computed on a placeholder.
- **B5.** The demo refusal is a typed literal.
- **B6.** The SymPy Part C rows conjoin domain predicates built from the pre-response amplitude `d` (seven
  occurrences in the `δ = 0` row).
- **B7.** The tag sets are not parallel across engines.
- **Exclusion.** `S9b_exports.py` is not added to any fold until a repaired delta is accepted (`bb94b885`).

## D7. Sequence and roles
1. **This list.** One Codex + Grok pass, folded once.
2. **Spec v10.** Codex (gpt-6.1-sol) authors it from v8, this list and the O2 sources. A fresh Claude agent and
   Grok review it until clear.
3. **Build directive.** A repair-and-extension build directive (decision list, G2) leads to one build round. It
   covers D6 and Part D, in both engines.
4. **Afterwards.** Build legs review until clear, with a FORM ablation each. Then the comparator, then the record.

## Not decided here
- Any expected value, sign or limit. This includes the relation of `c_comp` to `c₀`, and of `j_n` to zero.
- The laws for S8 inertia and stress, and O1/O3–O7. The exception is Part D, where P5, P6 and the density link
  supply content as scoped above.
