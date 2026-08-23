# S11b SymPy engine — X-1 repair directive (energy-basis independence omits §5's total-divergence equivalence)

## Authority and boundary
Repair `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py` IN PLACE (committed
baseline `864d6f41`). `CLAUDE.md` binds; `directives/S11b_SHARED_PHYSICS.md` is the sole physics authority —
implement it as written; ⛔ do not modify or relitigate the spec's physics. ⛔ Add no expected value, count,
or acceptance criterion (rule 5): the orchestrator
holds the withheld references and diffs on our side. This is a NARROW fix confined to the **energy-basis
construction** (Task B0-energy) and the coefficient assignment that flows from it; every other object —
the face response, the affinity law, the B2c static map, the dispersion route, roots, stability, the
transverse decoupling result, the breathing-mode form, and the entire S11 carry-forward — must stay
**physically** unchanged.

⚠⚠ **This repair edits the SHARED stored-energy density, so many DOWNSTREAM symbolic objects will legitimately
change their COORDINATES** (any expression carrying a folded energy coefficient re-expresses). ⛔ Do **not**
demand byte-identity of the equations of motion, `U_LONG`, `sigma_eta`, `mu_theta`, `p_W`, the dispersion
matrix, or the transverse stiffness. ⭐ The invariant to preserve is the **physics**, not the spelling: the
stored-energy functional **modulo a total in-plane divergence** (equivalently, the Euler–Lagrange equations
it generates), the dispersion **locus/roots**, the stability **classification**, transverse **decoupling**,
the B2c/G13 **slice relation**, and the breathing **form** — each as a function of whatever free coefficients
survive. Confine any "leave unchanged" claim to those PHYSICS invariants and to the NAMED non-target objects
below; ⛔ do not bless the current basis count or the current coefficient list as correct (that both leaks a
reference and can shield the defect).

## The defect — X-1
`construct_energy_basis()` (L463) judges basis independence with `independent_columns` (L406), which tests
**pointwise polynomial independence** of the invariants over the raw field-component symbols (the nine
`G_i_j`, `grad_theta_i`, `grad_e_W_i`). That test implements part of §5's symmetry group but **omits one
explicit bullet of it**:

> §5 (`directives/S11b_SHARED_PHYSICS.md:286-287`), inside "THE SYMMETRY GROUP, STATED IN FULL":
> *"Equivalence modulo total divergences — two densities differing by a total in-plane divergence are the
> same term; ⛔ do not count both."*

Two constructed invariants can be pointwise-distinct polynomials yet equal **as energy terms** (integrated
densities) because they differ by a total in-plane divergence `∂_i V_i`. The current independence test cannot
see that, so it can carry more invariants than §5's group admits, over-counting the basis and assigning a
free coefficient to a term that is not independent under §5.

## The fix — implement §5's full symmetry group in the basis judgment; re-express, ⛔ do not delete
1. **Independence modulo total divergence.** Make the energy-basis independence judgment honor §5's
   total-divergence-equivalence bullet, so `ENERGY_BASIS`, `ENERGY_BASIS_COUNT`,
   `ENERGY_BASIS_INDEPENDENT_TERMS`, and `ENERGY_BASIS_OMISSIONS` (emitted at L618-621) reflect the **full**
   §5 symmetry group and not merely pointwise independence. Two constructed invariants that differ by a total
   in-plane divergence `∂_i V_i` (for some bilinear `V_i` in the fields and their gradients) are **one** term.
   ⭐ Compute this equivalence — do not assert it: the standard CAS witness is the Euler–Lagrange operator,
   since a constant-coefficient quadratic density is a pure total divergence **iff** its EL derivative with
   respect to every field vanishes identically. Emit the objects that establish which constructed invariants
   collapse together; ⛔ do not type the result.
2. **Eliminate a redundant invariant by REWRITING it, folding its coefficient into the retained ones — ⛔
   never by dropping its term.** Setting a redundant invariant's coefficient to zero **deletes physics**
   (the functional loses that stiffness); the correct operation replaces the redundant invariant by its
   total-divergence-equivalent linear combination of the RETAINED invariants and adds its coefficient into
   theirs. The emitted stored-energy density (`full_energy`, L498) must remain the **same functional modulo a
   total in-plane divergence** as the pre-repair baseline — i.e. it must generate the **same Euler–Lagrange
   equations**.
3. **Prove preservation with a genuine two-route residual (no tautology).** Emit
   `ENERGY_REEXPRESSION_RESIDUAL` = (the full symmetry-enumerated density, one free coefficient per
   constructed invariant) − (the reduced re-expressed density, redundant coefficients folded in), and exhibit
   it as a total in-plane divergence by emitting its **Euler–Lagrange derivative with respect to each field**
   (⛔ emit the derivative object; ⛔ do not `assert` it zero). This residual must be able to **fail**: under a
   wrong fold (a one-sided corruption of one folded coefficient) the EL derivative must become nonzero. ⭐
   Verify independence of the two routes by one-sided corruption — corrupt only the enumeration route; the
   reduced route must not move.
4. **Scope of the independence change.** The fix belongs in the **un-substituted** basis construction
   (before any mode ansatz). ⚠ State explicitly whether the other `independent_columns` uses that operate on
   **mode-substituted** reductions — `impermeable_reduced`/`imp_pivots` (L624-625) and
   `flux_reduced`/`flux_pivots` (L632-633) — are affected. (Once a mode ansatz is substituted the density is
   algebraic in `eta,e_W`, so the divergence structure is generally gone; confirm this by computation and
   report it, ⛔ do not assume it.)

## Downstream rows the fix touches — carry them, ⛔ do not leave a stale reference
The energy coefficient assignment (`full_energy`, L498-510), the moduli tuple in `STABILITY_CONDITION`
(L1533), and the dimension map's per-coefficient rows (the energy-coefficient entries in the `dimensions`
dict, `symbol_links`, and `symbol_dimensions`, ~L1795-1830) all enumerate the energy coefficients by hand. If
the reduction removes an independent coefficient, every one of these must be brought into agreement with the
reduced basis **by construction** (⛔ never by a hand-typed edit that asserts the new list) — drive them from
the same reduced-basis object so they cannot drift, and ⛔ leave no dangling symbol for a folded-away
coefficient. State which tags legitimately change (that is expected, not a regression).

## Export-chain integrity — this engine is ON the export chain (⛔ non-negotiable)
The engine writes `S11b_exports.py`, which downstream steps import. The full inherited export contract is
G7 + the `F1`–`F9` chain rules (§11, `S11b_SHARED_PHYSICS.md:975-982`) — ⛔ implement it whole, ⛔ do not
paraphrase a selected subset (F9's own measured failure was a dropped clause). Preserve, verified by a
**real-consumer** run (`import S11b_exports` executes the module — ⛔ never read a cached `.out`):
- **Carry-forward fidelity.** Every key in `S11_exports.LEDGER` (1663 rows) appears in `S11b_exports.LEDGER`
  with an identical value, except any `F9c` collision (which prefixes `s11b_` and leaves the imported row
  untouched). ⛔ No S11 row silently dropped or mutated.
- **F9 applied WHOLE** (three-valued, TOTAL over the imported row shapes — ⛔ not a paraphrase); the `v_0`
  false-merge guard holds (the drain is `v_dr`; ⛔ the bare key `v_0` never appears in the S11b rows).
- **D1 — every emitted primary object reaches the flat export.** This repair adds NEW primary emissions
  (`ENERGY_REEXPRESSION_RESIDUAL`, its Euler–Lagrange derivatives, the collapse evidence) and changes the
  basis rows. ⚠ `import S11b_exports` and the `D3` round-trip validate only what was WRITTEN — a new row that
  is computed but never placed in the `LEDGER` passes both. So `D1` (iterate-the-emitted-collection) must
  carry **every** non-local primary emission into the flat export MECHANICALLY, with unique keys and complete
  class/dimension metadata (`F1`/`F3` in-row evidence). ⛔ No emitted primary object left out of the ledger.
- **§11 import identity (⛔ do not redeclare an import).** Every inherited object (`c_s0`, `mu_R`, `rho_br`)
  stays bound to its imported `S11_exports.LEDGER` object under the SAME name (§11 L949/L960/L961) — ⛔
  carry-forward fidelity of the upstream ROW is not enough; the repaired energy must still USE the imported
  object. If a coefficient fold would modify an imported energy coefficient, write the result as an EXPLICIT
  composite of the imported object plus the retained S11b knob(s); ⛔ never rename a folded/effective
  coefficient onto an imported symbol, and ⛔ never redeclare an import under a new identity. `rho_m`/`v_dr`
  originate here (⛔ no upstream row minted).
- **Digests** pin exactly `{own source, S11_exports.py, S11b_SHARED_PHYSICS.md}`; **D3** round-trip present;
  **`_RELATIONALS`** reviver present; export FROZEN (`MappingProxyType` outer + per-record, `del _LEDGER`);
  **F6** publish gate intact.

## The three script clauses + corollaries (verbatim, non-negotiable)
1. PRINT computed objects; ⛔ never a prose conclusion. 2. PRINT the residual; ⛔ do not `assert` it zero.
3. Interpretation → the step record. ⛔ A hand-typed CAS object is still hand-typed: every emitted object is
REACHED BY COMPUTATION; the physical symbols are combined by hand **only** in constructing the invariants and
the ansatz, and every control re-enters at that construction, ⛔ never at a result. ⛔ No tautological
residual — the two operands of `ENERGY_REEXPRESSION_RESIDUAL` come from INDEPENDENT routes (full enumeration
vs reduced re-expression), verified by one-sided corruption (break the enumeration route only; the reduced
route must not move).

## Supplied vs withheld
SUPPLIED (use freely): §5's full symmetry group including the total-divergence bullet; the field content and
the constructed invariants already in the engine; the entire S11 import. WITHHELD (rule 5): the resulting
basis count, which invariant is redundant, the fold coefficients, and any identity among the invariants — the
orchestrator computes these and diffs against the engine's emitted objects; a mismatch is a **finding**, ⛔
not a build failure. ⛔ Do not iterate the engine toward any count or coefficient value.

## Run discipline (Python)
This is a pure-Python SymPy engine — ⛔ no Mathematica, no kernel/seat limits. Run the engine yourself to a
scratch log; long-but-printing is fine, long-and-silent is the failure. The orchestrator captures and commits
the `.out`; ⛔ do not commit.

## Acceptance — executable, DECISIVE, no expected values (rule 5)
- **Independence honors §5's full group.** By one-sided corruption: an invariant that is a
  total-divergence-equivalent of retained invariants is **no longer** counted independent (adding it does not
  raise the emitted basis rank/count), while a genuinely independent invariant still **is** counted. ⛔ Do not
  state the resulting count.
- **Physics preserved (the load-bearing check).** `ENERGY_REEXPRESSION_RESIDUAL`'s Euler–Lagrange derivative
  is identically zero (the reduced density differs from the full enumeration by a total in-plane divergence),
  and this residual FAILS (EL derivative ≠ 0) under a corrupted fold. The dispersion locus/roots, stability
  classification, transverse decoupling (coupling identically zero), the B2c/G13 slice relation, and the
  breathing-mode form are behaviourally unchanged as functions of the surviving free coefficients (the
  orchestrator diffs these against the withheld references and against the sibling engine).
- **Export integrity.** `import S11b_exports` runs; the 1663 S11 rows are carried identically; F9 whole; the
  digests, D3, `_RELATIONALS`, freeze, and F6 gate are intact.
- **No value leak.** ⛔ No numeric count, coefficient relation, or "expected" appears anywhere in the engine.

## Report (§13) — under 15 lines
The edits (site:line), the computed evidence that the independence judgment now honors §5's total-divergence
equivalence, the `ENERGY_REEXPRESSION_RESIDUAL` two-route construction and the one-sided-corruption
demonstration that it can fail, which downstream tags legitimately changed coordinates (and confirmation the
physics invariants did not), and the export-integrity real-consumer result. ⛔ State no count or coefficient
value.
