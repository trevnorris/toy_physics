# S15 reading: how the momentum crossing the region around drain 2 depends on the separation (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-26. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The lookups behind this list are in
`directives/_measurements/S15_reading_setup_lookups.md` ("L-n").

**Gate.** This list gets one Codex + Grok pass and is folded once (G2). It also names the objects the build computes,
so that content is reviewed until clear before any builder (spec row).

**Revision 2** folds pass 1 once (`directives/_measurements/S15_reading_decision_list_review_p1_disposition.md`,
F1–F4). The reviewed revision 1 is preserved at `1c6695df`.

## Why this build
- The user's S15 choice **"Keep it general"** reads, verbatim: "Any radial inflow and density profile around each
  drain. S15 then says which profiles give an inverse-square pull, and that becomes gravity's flow requirement,
  compared with light's in S21." (L1).
- S15's list leaves that reading to the record: "Which profiles give which dependence on the separation is the
  record's reading of the printed objects. The build prints the objects." (`directives/S15_decision_list.md:45–46`,
  L2). The accepted build prints the surface and region integrals (L4). Their dependence on the separation is not
  printed, and stating it needs computation.
- **User decision (2026-10-26).** The user chose "Reading build (Recommended)", whose option text reads: "Both engines
  independently compute how the named pieces depend on the separation, for general profiles. A comparator then diffs
  them, and the record reads the compared output. This keeps the two-engine check."

## Supplied (unfalsifiable in this build)
- Everything S15's list supplies, unchanged: the sections "Supplied equations" and "The drains, their fields and the
  region", and the domain and pressure conditions (`:120–121`). This is revision 3, unchanged since `91248b67` (L2).
- The objects of that list's "What the build computes", items 1 and 2, as it names them. The accepted S15 build
  printed them (`3192123b`), and the accepted comparator compared them (`5a019394`) (L3).
- Each engine constructs them again from the supplied equations. The `.wl` imports nothing; it may reuse the text of
  `mathematica/S15_two_drain_region.wl`. Neither engine reads the other engine's code or output.

## The separation and the family (supplied definitions)
- **The separation:** `D = X_1 − X_2`, `d = |D|`. The domain requires `d > R` (L2).
- **The family** (a labelled choice). At one instant, configurations that differ in drain 1's centre, `X_1`.
  - **Held** (labelled): the region and its motion (`X_2`, `Ẋ_2`, `r_in`, `R` and their rates), and the form of
    `p_br`.
  - **Changing:** every other quantity: `u_k`, `ρ_k`, `W`, `ϱ`, `η`, `ζ`, `Π_o`, `f` and `Ẋ_1`. Each one's change
    along the family is a general live object, with the same arguments as the quantity.
- **The separation derivative** of an object `Q`:
  `𝒟Q = d/dε Q(X_1 + εD, F + εF′, …)|_{ε=0}`, where each changing quantity `F` goes to `F + εF′` and `F′` is its
  change.
- **The departure from the inverse square:** `E_Q = 𝒟Q + 2Q`. `Q` falls off as the inverse square of the separation
  along the family exactly where `E_Q = 0` for every `D` on the domain.
- **Its two parts:** `E_Q = E_Q^geo + E_Q^resp`, where
  - `E_Q^geo = (𝒟Q)|_{F′=0} + 2Q` is the departure when no quantity other than `X_1` changes along the family;
  - `E_Q^resp = E_Q − E_Q^geo` is what the changes add.

  The parts are these objects, whatever form an engine writes a term in.
- **Where they exist.** Each departure object is stated where the separation derivative of its object exists, and it
  carries the condition under which it does, labelled where it is printed. For instance, the separation derivative
  of the whole pressure, and of the whole minus the expansion, needs `p_br` differentiable at the total density. The
  response part of the expansion needs `p_br′` differentiable at `ρ_2 + ϱ`.

## The profile family (a labelled choice, items 3 and 4 only)
`u_k(r, t) = A_k(t) r^{a_k}` and `ρ_k(r, t) = B_k(t) r^{b_k}`, for `k = 1, 2`. The amplitudes `A_k`, `B_k` are general
live functions of `t`. The exponents `a_k`, `b_k` are real constants, kept as symbols. Everything else stays general,
the changes `F′` included.

## What the build computes
Both engines construct these independently, with the `.wl` blind.

**The objects `Q`:**
- every piece and group S15's item 1 prints, at each surface, for each Cartesian component (L4);
- the region's momentum and its rate, each whole and by group, and the integrated force density (item 2; L4);
- at each surface, three named sums: the rest part (carried, whole pressure, viscous, other), the transport (centre
  and radius), and the crossing (rest plus transport); and the same three sums over both surfaces.

**For each `Q`:**
1. `Q` itself.
2. `E_Q^geo` and `E_Q^resp`, for general profiles, exact wherever they exist on the domain, with every live quantity
   kept.
3. `E_Q^geo` on the profile family, with its dependence on `d` explicit and exact. Where the integrand holds no general
   field besides the family's amplitudes, the region's quantities and the velocities, the integrals are evaluated in
   closed form. Otherwise the object is printed with the family substituted and left unevaluated.
4. For each `Q` whose item-3 object is in closed form:
   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
   - `Q` on the family under each condition.
5. For each drain `k`, the part of `Q` that depends on that drain's velocity: `Q` minus `Q` at `Ẋ_k = 0`.
6. Dimensional homogeneity of items 2 and 5, each change `F′` carrying the dimension of `F`.

Every object is exact on S15's domain wherever it exists, and it carries the conditions under which it exists: S15's
pressure conditions, and for each departure object the conditions above. No limit or expansion is taken
in place of an exact object. Where an engine does not obtain an object exactly, it prints a label saying so, beside
what it did print.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value. The build computes no force law, coefficient,
normalization or mass, and identifies no printed object with a pull or with the force on a throat.

## Model point
- **Regime:** none chosen; as S15's list states it (`:129–141`). The family compares configurations at one instant,
  and every object keeps its time dependence.
- **Freezes:** none intended. Every quantity that is not held changes along the family, and its change stays live.
- **Labelled choices:** S15's, plus the family, its held quantities, and the profile family (items 3 and 4).
- **Transfer:** as S15's list states it.

## For the comparator and the record (not in the build directive)
- **The comparator** compares items 1–6 across engines. It also compares item 1 with the accepted S15 transcripts.
- **The record.** Its decision list is written after this build's comparator is accepted. It reads S15's objects and
  this build's, as S15's list asks (`:145–156`):
  - for each `Q`, the printed conditions under which `E_Q^geo` vanishes. These are the conditions under which `Q`
    falls off as the inverse square when no quantity other than `X_1` changes along the family. Beside them it reads
    `E_Q^resp`, what any change adds. Neither is read as a condition for `E_Q = 0` while the changes are live;
  - gravity's flow requirement for S21, stated as those conditions, with `E_Q^resp` and the conditions of existence
    beside them;
  - the velocity-dependent parts.

  It checks the prior-art rows as S15's list says (`:160–163`).

## Retained OPEN
- Everything S15's list retains as OPEN (`:174–186`).
- For general profiles, the conditions under which `E_Q^geo` vanishes. Item 2 prints it; this build does not solve it.
- Whether, and how, any quantity changes with the separation. `E_Q^resp` prints what such a change would add. No
  change is modelled, so no condition for `E_Q = 0` with the changes live is computed.

## Intended claim
For each object S15 prints: its departure from the inverse square of the separation, exact, for general profiles,
with every quantity's change along the family live, wherever it exists. On the power-law profile family, the
conditions under which each closed-form geometric departure vanishes: the conditions under which the object falls off
as the inverse square when no quantity other than `X_1` changes along the family. The claim is conditional on P5, P6,
`R-S14a-01` and the labelled choices. It gives no force law and no pull, and no printed object is the force on a
throat.

## Next
Pass 1 is folded once, here. The spec content changed, so it gets a second Codex + Grok pass, and is reviewed until
clear. Then the build (Codex) and its G1 review (fresh Claude + Grok) until clear; the comparator; and the S15 record
list.
