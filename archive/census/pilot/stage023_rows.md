# Census rows — stage023, both engines

**Artifact class:** census rows per `notes/census/CENSUS_SCHEMA.md`. **Stage:** 023.
**Builder leg.** ⛔ No merge below is adjudicated (§8.3): every cross-artifact QID is **PROPOSED**.

## 0. The two artifacts, and the locus abbreviation

| handle | full path |
|---|---|
| **PY** | `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py` |
| **WL** | `/var/projects/toy_physics/research/pde_ledger_v2/mathematica/ledger_stage023_nullspace_underdetermination_mathematica_audit.wl` |
| **NOTE** | `/var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage023_nullspace_underdetermination.md` |

⚠ **Stated deviation from §7.2's locus rule.** §7.2 requires each locus to carry its own full path. Every
`PY:NNN` / `WL:NNN` / `NOTE:NNN` below expands to exactly one of the three paths above and nothing else, so
no line number in this file is ambiguous. Recorded in §J so it is auditable rather than silent.

**Off-stage evidence loci** (all opened and read for this pass — §9.1 rule 1; none carried from an attribution):

| handle | full path |
|---|---|
| REG | `/var/projects/toy_physics/research/pde_ledger_v2/notes/parameter_register.md` |
| PREREG | `/var/projects/toy_physics/docs/pathA_preregistration.md` |
| MAP | `/var/projects/toy_physics/docs/model_map.md` |
| S008 | `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py` |
| S009 | `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage009_flat_slab_return_residual_sympy_audit.py` |
| S017 | `/var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage017_grouped_p2_lane_isotropy.md` |
| S022 | `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage022_cross_l_fingerprints_sympy_audit.py` |
| S043 | `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage043_irreducible_count_range_sympy_audit.py` |
| SRCMAP | `/var/projects/toy_physics/research/pde_ledger_v2/notes/stage023_pathA34_nullspace_underdetermination_source_map.md` |

---

## A. The reported-result set (§5.8 output 17) — the starting set, taken as fixed

From `notes/census/REPORTED_RESULTS.md` §2.1. ⛔ Not amended (§7.1.1). Three results:

| RES id | value-bearing loci used as the walk's entry points |
|---|---|
| `RES:023:RETURN-UNDERDETERMINED` | `PY:388-393` (rank0 / native nullity / return-augmented rank / return-moving nullity), `PY:395-413` (the two `e_{Z_ℓ,ret}` witnesses, their constraint deltas and `ΔT0`/`ΔT1`), `PY:205-217` (the 11 dofs), `PY:384` (the three collected constraints); `.wl` counterpart `WL:143-151` |
| `RES:023:A0-FORM` | `PY:342` (`A0`), `PY:345` (`expected_A0`), `PY:349` (residual), `PY:351` (realized ω order), condition `PY:273-280`, `:316`; `.wl` counterpart `WL:199,201,203,640` |
| `RES:023:A1-FORM` | `PY:343`, `PY:346-348`, `PY:350`, `PY:352`; `.wl` counterpart `WL:200,202,204,641` |

⛔ **Verdict-token loci are NOT entry points.** REPORTED_RESULTS also lists `PY:862` and `PY:1036` against
`RES:023:RETURN-UNDERDETERMINED`; both are `verdict_residual(baseline["verdict"], …)` asserts. §7.1.1 forbids
a verdict token as the closure walk's entry point while allowing it as test-(a) evidence, so this census uses
them **only** as test-(a) evidence and walks from the computed loci above. **This single adjudication is what
puts the whole dimensional block, the provenance-class firewall and the decoupling detector out of the
universe** (§F) — it is the highest-leverage decision in this pass and is flagged for the review leg.

**.wl reachability adjudication.** REPORTED_RESULTS cites only `.py` loci. §7.2 states the `.py` and `.wl` of
one stage are two artifacts and that "a dual-engine value has at least two occurrences", and `NOTE:839-841`
records the two engines emitting the same `native_nullspace_dimension=8`, `return_moving_nullity=2` and
`A0/A1` residual-zero. The `.wl` rows therefore answer to the **same three results**, reached through the
`.wl`'s own value-bearing content (`WL:146-148`, `WL:203-204`). Recorded as an adjudication, not a given.

---

## B. Granularity rules applied (stated so they can be contested)

- **(G1)** An **unnamed** numeric literal inside an expression is a closure **leaf**, never its own occurrence
  — route 1 admits *named* quantities (§7.1). (So the `2` in `Keta + 2*TOmega` and the `2` in
  `2*c_s**3` are `C-SELF` leaves, not rows. Both are flagged in the owning rows' notes.)
- **(G2)** Binding site = the file-and-line of the assignment **or assertion** expression (§7.2). In a `.py`
  every assignment line is a binding site whatever its scope (§7.2's Wolfram rules are stated as *differing*
  from the `.py`). In the `.wl`, `Module`/`With` locals belong to the enclosing definition head (§7.2 rule 2)
  and a pattern definition is an occurrence at its definition line where the artifact names its result
  (§7.2 rule 3).
- **(G3)** A dict/association **entry that introduces a value** is a binding site; an entry that merely
  re-exports an already-bound local is not.
- **(G4)** A multi-target assignment (`K0c, Keta, TOmega = sp.symbols(...)`) is ONE binding site and N
  occurrences — the triple is `(artifact, binding-site, quantity)` (§7.2).
- **(G5)** `$Assumptions` (`WL:26`) is a binding site (§7.2 rule 4) carrying one occurrence per quantity it
  constrains; it is the `.wl`'s only declaration site for the free symbols.

**LIVE/RETIRED (§7.4): every row below is `LIVE`.** Checked, not assumed: the register's only two retired
rows are `λ_Pu` (`REG:139`) and `α_aniso` (`REG:159`), neither of which names any stage023 quantity; the
stage023 edge `R42` is carried live at `REG:309`; `NOTE` carries no retirement marker.

---

## C. Rows — artifact **PY** (63 occurrences)

Occurrence ID = `QID:<name>@PY#<line>`. `PF` = `PHYSICS-FED`; `SR` = `SELF-REFERENTIAL`;
`CONV` = `CONVENTION-LADEN`; `SB` = `should_be_tier`; `basis` = `should_be_basis`; `Δ` = delta flag.
Evidence codes `E1…E12` resolve in §E. Reachability witness codes `W1…W3` resolve in §B/§A:
**W1** = `RES:023:RETURN-UNDERDETERMINED`, **W2** = `RES:023:A0-FORM`, **W3** = `RES:023:A1-FORM`;
hop kind `expr` = expression walk, `premise` = `premise-dependence` hop (§3.3.1(4)).

### C.1 Free-symbol declarations — axis B `B-DECLARED-UNASSIGNED` (16 rows)

Every row: axis B = `B-DECLARED-UNASSIGNED` (evidence **E1**); closure = `{C-FREE}`; `PF = C-UNRESOLVED`;
`SR = false`; near-miss buckets: none.

| # | QID | site | axis A | CONV | is_tier | SB | basis | Δ | witness | ev |
|---|---|---|---|---|---|---|---|---|---|---|
| PY-01 | `QID:omega` | PY:163 | `A-UNADJUDICATED` | false | `no-tier:unadjudicated` | `no-tier:unadjudicated` | `none` | N | W1/W2/W3 expr | E1,E2 |
| PY-02 | `QID:a` | PY:164 | `A-REDUCIBLE-UNDERIVED` | **UNADJUDICATED** | `tier1-debt` | `no-tier:convention` | `convention-candidate` | **Y** | W2/W3 expr | E1,E3 |
| PY-03 | `QID:c_s` | PY:165 | `A-UNADJUDICATED` | false | `no-tier:unadjudicated` | `tier1-debt` | `named-route` | **Y** | W2/W3 expr | E1,E4 |
| PY-04 | `QID:M0` | PY:166 | `A-UNADJUDICATED` | false | `no-tier:unadjudicated` | `tier3-emergent` | `physical-picture-expectation` | **Y** | W2 expr | E1,E5 |
| PY-05 | `QID:D1` | PY:167 | `A-UNADJUDICATED` | false | `no-tier:unadjudicated` | `tier3-emergent` | `physical-picture-expectation` | **Y** | W3 expr | E1,E5 |
| PY-06 | `QID:D0` | PY:169 | `A-UNADJUDICATED` | false | `no-tier:unadjudicated` | `tier1-debt` | `named-route` | **Y** | W1 expr | E1,E6 |
| PY-07 | `QID:K0c` | PY:170 | `A-REDUCIBLE-UNDERIVED` | false | `tier1-debt` | `tier3-emergent` | `named-route` | **Y** | W1/W2 expr | E1,E7 |
| PY-08 | `QID:K_eta` | PY:170 | `A-REDUCIBLE-UNDERIVED` | false | `tier1-debt` | `tier3-emergent` | `named-route` | **Y** | W1/W3 expr | E1,E7 |
| PY-09 | `QID:T_Omega` | PY:170 | `A-REDUCIBLE-UNDERIVED` | false | `tier1-debt` | `tier3-emergent` | `named-route` | **Y** | W1/W3 expr | E1,E7 |
| PY-10 | `QID:Z0_ret` | PY:171 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1/W2 expr | E1,E8 |
| PY-11 | `QID:Z1_ret` | PY:171 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1/W3 expr | E1,E8 |
| PY-12 | `QID:Omega_U` | PY:175 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 expr | E1,E9 |
| PY-13 | `QID:Omega_W` | PY:175 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 expr | E1,E9 |
| PY-14 | `QID:R_mix` | PY:175 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 expr | E1,E9 |
| PY-15 | `QID:g_U` | PY:175 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 expr | E1,E9 |
| PY-16 | `QID:g_W` | PY:175 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 expr | E1,E9 |

⚠ **PY-01** is recorded **BLOCKED on `should_be_basis`** — see §J item 1. ⚠ **PY-07/08/09** carry an
**intra-occurrence conflict** (§10.3) — see §H.

### C.2 Computed symbolic expressions — axis B `B-EXECUTED` (41 rows)

Every row in this block: axis B = `B-EXECUTED` (evidence **E10**); axis A = `A-REDUCED` (evidence **E11** —
the reduction is performed *at the binding site itself* and the loci it reduces TO are the "reduces to"
column); closure carries at least one `C-FREE` leaf and **no** `C-FIELDEQ`/`C-EXTERNAL`/`C-PEER` leaf, so
`PF = C-UNRESOLVED` (evidence **E12**); `SR = false`; `CONV = false`;
`is_tier = no-tier:unadjudicated` (§5.7: `A-REDUCED` fails tier 3's `PHYSICS-FED` conjunct and has no other
tier); `SB = tier3-emergent`; `basis = physical-picture-expectation`; `Δ = Y`.
⛔ These rows are in **no** §4 near-miss bucket: `C-UNRESOLVED` may not be asserted into
`executed-but-not-physics-fed` (§3.3).

| # | QID | site | reduces to (closure, all leaves) | witness |
|---|---|---|---|---|
| PY-17 | `QID:omega_u` | PY:221 | `C-FREE{Omega_U@PY:175}` + `C-SELF{control flag @PY:202}` | W1 expr |
| PY-18 | `QID:P_port` | PY:222 | `C-FREE{Omega_U,g_W,R_mix,g_U @PY:175}` | W1 expr |
| PY-19 | `QID:Delta_port` | PY:223 | `C-FREE{Omega_U,Omega_W,R_mix @PY:175}` | W1 expr |
| PY-20 | `QID:N0_from_port` | PY:224 | via PY:222,PY:223 → `C-FREE{Omega_U,Omega_W,R_mix,g_U,g_W}` | W1 expr |
| PY-21 | `QID:P0_raw` | PY:225 | via PY:224 + `C-FREE{D0@PY:169}` | W1 expr |
| PY-22 | `QID:K0` | PY:284 | `C-FREE{K0c@PY:170}` | W1/W2 expr |
| PY-23 | `QID:K1` | PY:285 | `C-FREE{K_eta,T_Omega@PY:170}` + `C-SELF{literal 2 @PY:285}` | W1/W3 expr |
| PY-24 | `QID:Z0` | PY:287 | `C-FREE{Z0_ret@PY:171}` | W1/W2 expr |
| PY-25 | `QID:Z1` | PY:288 | `C-FREE{Z1_ret@PY:171}` | W1/W3 expr |
| PY-26 | `QID:T0` | PY:302 | via PY:284,PY:287 → `C-FREE{K0c,Z0_ret}` | W1/W2 expr |
| PY-27 | `QID:T1` | PY:303 | via PY:285,PY:288 → `C-FREE{K_eta,T_Omega,Z1_ret}` + `C-SELF{2}` | W1/W3 expr |
| PY-28 | `QID:epsilon0` | PY:304 | via PY:284,PY:287 → `C-FREE{K0c,Z0_ret}` | W2 expr |
| PY-29 | `QID:epsilon1` | PY:305 | via PY:285,PY:288 → `C-FREE{K_eta,T_Omega,Z1_ret}` + `C-SELF{2}` | W3 expr |
| PY-30 | `QID:positive_bounded` | PY:316 | via PY:273-280 → `C-SELF{positivity assumptions @PY:170,:171}` + `C-FREE{K0c,K_eta,T_Omega,Z0_ret,Z1_ret}` | W2/W3 expr |
| PY-31 | `QID:A0` | PY:342 | `C-SELF{v0@PY:337}` + `C-CONVENTION{sp.I, time convention @PY:612}` + `C-FREE{a,omega,c_s,M0}` + via PY:302 | W2 expr |
| PY-32 | `QID:A1` | PY:343 | `C-SELF{v1@PY:338, power1@PY:339}` + `C-CONVENTION{sp.I}` + `C-FREE{a,omega,c_s,D1}` + via PY:303 | W3 expr |
| PY-33 | `QID:expected_A0` | PY:345 | `C-SELF{typed pathA_29 form, uncited}` + `C-CONVENTION{sp.I}` + `C-FREE{a,omega,M0,c_s}` + via PY:304 | W2 expr |
| PY-34 | `QID:expected_A1` | PY:346 | `C-SELF{typed form + literal 2, uncited}` + `C-CONVENTION{sp.I}` + `C-FREE{a,omega,D1,c_s}` + via PY:305 | W3 expr |
| PY-35 | `QID:A0_residual` | PY:349 | via PY:342, PY:345 | W2 expr |
| PY-36 | `QID:A1_residual` | PY:350 | via PY:343, PY:346 | W3 expr |
| PY-37 | `QID:A0_order` | PY:351 | via PY:342 | W2 expr |
| PY-38 | `QID:A1_order` | PY:352 | via PY:343 | W3 expr |
| PY-39 | `QID:GENERATOR_DOFS` | PY:205 | `C-FREE{the 11 symbols @PY:169,:170,:171,:175}` | W1 expr |
| PY-40 | `QID:dofs` | PY:381 | via PY:205 | W1 expr |
| PY-41 | `QID:base_constraints` | PY:384 | via PY:225 + `C-FREE{K0c,K_eta,T_Omega}` + `C-SELF{2}` | W1 expr |
| PY-42 | `QID:constraints` | PY:386 | via PY:384 (baseline: selector list empty, PY:385) | W1 expr |
| PY-43 | `QID:constraint_jacobian` | PY:387 | via PY:386, PY:381; routine `sp.diff` PY:370 contributes no leaf (§3.3.1(1)) | W1 expr |
| PY-44 | `QID:rank0` | PY:388 | via PY:387; `Matrix.rank()` PY:375 is a generic linear-algebra routine ⇒ **no leaf** (§3.3.1(1)) | W1 expr |
| PY-45 | `QID:native_nullspace_dimension` | PY:389 | via PY:388, PY:381 | W1 expr |
| PY-46 | `QID:grad_T0` | PY:390 | via PY:302, PY:381 | W1 expr |
| PY-47 | `QID:grad_T1` | PY:391 | via PY:303, PY:381 | W1 expr |
| PY-48 | `QID:return_augmented_rank` | PY:392 | via PY:387, PY:390, PY:391 | W1 expr |
| PY-49 | `QID:return_moving_nullity` | PY:393 | via PY:392, PY:388 | W1 expr |
| PY-50 | `QID:witness_unit_vector` | PY:397 | `C-SELF{1,0}` + via PY:381 (dof identity) | W1 expr |
| PY-51 | `QID:constraint_deltas` | PY:398 | via PY:387, PY:397 | W1 expr |
| PY-52 | `QID:delta_T0` | PY:402 | via PY:390, PY:397 | W1 expr |
| PY-53 | `QID:delta_T1` | PY:403 | via PY:391, PY:397 | W1 expr |
| PY-54 | `QID:moves_return` | PY:411 | via PY:402, PY:403 | W1 expr |
| PY-55 | `QID:native_null_moves_return` | PY:422 | via PY:393 | W1 expr |

Two further `B-EXECUTED` rows recompute the witness content at the assertion site (REPORTED_RESULTS cites
`PY:826-829` as the emission):

| # | QID | site | reduces to | witness |
|---|---|---|---|---|
| PY-56 | `QID:witness_constraint_residuals` | PY:821 | via PY:387, PY:397 | W1 expr |
| PY-57 | `QID:preserves_every_constraint` | PY:825 | via PY:821 | W1 expr |

*(PY-17 … PY-55 = 39 entries, plus PY-56 and PY-57 = **41 rows** in §C.2. Tally in §C.5.)*

### C.3 Declared literals — axis B `B-DECLARED-LITERAL` (3 rows)

Every row: axis B = `B-DECLARED-LITERAL`; closure = `{C-SELF}` **fully determined** ⇒ `PF = false`
(a positive finding, not `C-UNRESOLVED`); `SR = false`; `CONV = false`; axis A = `A-REDUCED` (evidence
**E13**); `is_tier = no-tier:unclassified-nonfed` (§5.1); `SB = tier3-emergent`; `basis = named-route`;
`Δ = Y`. ⛔ In **no** §4 near-miss bucket — not `B-EXECUTED`, and not `PHYSICS-FED`, so near-miss 3
(`physics-fed-but-declared-literal`) does **not** apply. **The single thing standing between these rows and
near-miss 3 is a missing citation locus** (§3.3's `C-PEER` rule): the value is genuinely computed at
`S022:209-210`, and neither PY nor NOTE carries a file-and-line for it.

| # | QID | site | value | axis A evidence | witness |
|---|---|---|---|---|---|
| PY-58 | `QID:v0` | PY:337 | `1` | E13 | W2 expr |
| PY-59 | `QID:v1` | PY:338 | `1/2` | E13 | W3 expr |
| PY-60 | `QID:power1` | PY:339 | `3` | E13 | W3 expr |

### C.4 Constitutive propositions (§7.1 route 2) — axis B `B-POSTULATED` (3 rows)

Every row: axis B = `B-POSTULATED` (evidence **E14**; ⛔ no code locus owed — §3.2); closure terminates at
the artifact's own declaration ⇒ `{C-SELF}`, `PF = false` (§3.3.1(4)); `SR = false`; hop kind
= `premise-dependence`. Route-2 bounds evidenced at **E15**/**E16**.

| # | QID | site | axis A | CONV | is_tier | SB | basis | Δ | witness | ev |
|---|---|---|---|---|---|---|---|---|---|---|
| PY-61 | `QID:pathA29_Z_is_premise` | PY:425 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 premise | E8,E14,E15 |
| PY-62 | `QID:pathA29_boundary_dof_none` | PY:425 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 premise | E8,E14,E15 |
| PY-63 | `QID:time_convention_exp_minus_i_omega_t` | PY:612 | `A-UNADJUDICATED` | **UNADJUDICATED** | `no-tier:unadjudicated` | `no-tier:convention` | `convention-candidate` | **Y** | W2/W3 premise | E14,E16 |

### C.5 PY tally

16 (C.1) + 41 (C.2) + 3 (C.3) + 3 (C.4) = **63 occurrences**.

---

## D. Rows — artifact **WL** (56 occurrences)

### D.1 Free-symbol declarations at `$Assumptions` — `B-DECLARED-UNASSIGNED` (16 rows)

Binding site `WL:26` for all sixteen (§7.2 rule 4 + G5). Axis A, CONV, is_tier, SB, basis, Δ and evidence
are **identical to the corresponding PY row** (same quantity, same off-stage evidence); axis B evidence is
**E1′** (declaration at `WL:26-36`; no top-level `Set` of any of them anywhere in WL — verified by reading
the file: the only bindings of these names are `Module` locals inside `caseFor`, `WL:430-450`).
Closure `{C-FREE}`, `PF = C-UNRESOLVED`, `SR = false`.

| # | QID | mirrors | axis A | is_tier |
|---|---|---|---|---|
| WL-01 | `QID:omega` | PY-01 | `A-UNADJUDICATED` | `no-tier:unadjudicated` |
| WL-02 | `QID:a` | PY-02 | `A-REDUCIBLE-UNDERIVED` | `tier1-debt` |
| WL-03 | `QID:c_s` (`cs`) | PY-03 | `A-UNADJUDICATED` | `no-tier:unadjudicated` |
| WL-04 | `QID:M0` | PY-04 | `A-UNADJUDICATED` | `no-tier:unadjudicated` |
| WL-05 | `QID:D1` | PY-05 | `A-UNADJUDICATED` | `no-tier:unadjudicated` |
| WL-06 | `QID:D0` | PY-06 | `A-UNADJUDICATED` | `no-tier:unadjudicated` |
| WL-07 | `QID:K0c` | PY-07 | `A-REDUCIBLE-UNDERIVED` | `tier1-debt` |
| WL-08 | `QID:K_eta` (`Keta`) | PY-08 | `A-REDUCIBLE-UNDERIVED` | `tier1-debt` |
| WL-09 | `QID:T_Omega` (`TOmega`) | PY-09 | `A-REDUCIBLE-UNDERIVED` | `tier1-debt` |
| WL-10 | `QID:Z0_ret` (`Z0ret`) | PY-10 | `A-IRREDUCIBLE-STRUCTURAL` | `tier1-structural` |
| WL-11 | `QID:Z1_ret` (`Z1ret`) | PY-11 | `A-IRREDUCIBLE-STRUCTURAL` | `tier1-structural` |
| WL-12 | `QID:Omega_U` (`OmegaU`) | PY-12 | `A-IRREDUCIBLE-STRUCTURAL` | `tier1-structural` |
| WL-13 | `QID:Omega_W` (`OmegaW`) | PY-13 | `A-IRREDUCIBLE-STRUCTURAL` | `tier1-structural` |
| WL-14 | `QID:R_mix` (`Rmix`) | PY-14 | `A-IRREDUCIBLE-STRUCTURAL` | `tier1-structural` |
| WL-15 | `QID:g_U` (`gU`) | PY-15 | `A-IRREDUCIBLE-STRUCTURAL` | `tier1-structural` |
| WL-16 | `QID:g_W` (`gW`) | PY-16 | `A-IRREDUCIBLE-STRUCTURAL` | `tier1-structural` |

⛔ Every `mirrors` entry is a **PROPOSED** merge only (§8.3); until the physics leg adjudicates, each pair is
two QIDs. §G reports both readings.

### D.2 Computed expressions — `B-EXECUTED` (37 rows)

Same block treatment as §C.2: axis A `A-REDUCED` (E11), `PF = C-UNRESOLVED` (E12), `SR = false`,
`CONV = false`, `is_tier = no-tier:unadjudicated`, `SB = tier3-emergent`,
`basis = physical-picture-expectation`, `Δ = Y`, no §4 near-miss bucket.

| # | QID | site | reduces to (closure) | witness |
|---|---|---|---|---|
| WL-17 | `QID:omega_u` (`ou`) | WL:105 | `C-FREE{OmegaU@WL:26}` + `C-SELF{factor arg}` | W1 expr |
| WL-18 | `QID:P_port` (`PPort`) | WL:105 | `C-FREE{OmegaU,gW,Rmix,gU}` | W1 expr |
| WL-19 | `QID:Delta_port` (`DeltaPort`) | WL:105 | `C-FREE{OmegaU,OmegaW,Rmix}` | W1 expr |
| WL-20 | `QID:N0_from_port` (`N0Port`) | WL:105 | via WL-18, WL-19 | W1 expr |
| WL-21 | `QID:P0_raw` (`P0Raw`) | WL:105 | via WL-20 + `C-FREE{D0}` | W1 expr |
| WL-22 | `QID:baseline_port_packet` | WL:120 | via WL:105 (`portKernel[1]`) | W1 expr |
| WL-23 | `QID:P0_raw` (2nd site, `P0port`) | WL:121 | via WL-22 | W1 expr |
| WL-24 | `QID:K1` (`K1dc`) | WL:122 | `C-FREE{Keta,TOmega}` + `C-SELF{literal 2}` | W1/W3 expr |
| WL-25 | `QID:GENERATOR_DOFS` (`rankDofs`) | WL:123 | `C-FREE{the 11 symbols @WL:26}` | W1 expr |
| WL-26 | `QID:base_constraints` | WL:124 | via WL-23, WL-24 + `C-FREE{K0c}` | W1 expr |
| WL-27 | `QID:constraint_jacobian` (rule) | WL:126 | via its arguments; `D[…]` contributes no leaf | W1 expr |
| WL-28 | `QID:null_basis` (rule) | WL:132 | via WL-27; `NullSpace` is a generic routine ⇒ **no leaf** (§3.3.1(1)) | W1 expr |
| WL-29 | `QID:return_gradients` (rule) | WL:134 | via WL-27 | W1 expr |
| WL-30 | `QID:return_moving_nullity` (rule) | WL:136 | via WL-28, WL-29; `MatrixRank` ⇒ **no leaf** | W1 expr |
| WL-31 | `QID:T0` (`T0dc`) | WL:138 | `C-FREE{K0c,Z0ret}` | W1/W2 expr |
| WL-32 | `QID:T1` (`T1dc`) | WL:139 | via WL-24 + `C-FREE{Z1ret}` | W1/W3 expr |
| WL-33 | `QID:epsilon0` (`eps0`) | WL:140 | `C-FREE{Z0ret,K0c}` | W2 expr |
| WL-34 | `QID:epsilon1` (`eps1`) | WL:141 | via WL-24 + `C-FREE{Z1ret}` | W3 expr |
| WL-35 | `QID:constraint_jacobian` (2nd site, `Jbase`) | WL:143 | via WL-27, WL-26, WL-25 | W1 expr |
| WL-36 | `QID:null_basis` (2nd site, `basis`) | WL:144 | via WL-28, WL-26, WL-25 | W1 expr |
| WL-37 | `QID:return_gradients` (2nd site, `Greturn`) | WL:145 | via WL-29, WL-31, WL-32, WL-25 | W1 expr |
| WL-38 | `QID:rank0` | WL:146 | via WL-35 | W1 expr |
| WL-39 | `QID:native_nullspace_dimension` | WL:147 | via WL-36 | W1 expr |
| WL-40 | `QID:return_moving_nullity` (2nd site) | WL:148 | via WL-36, WL-37 | W1 expr |
| WL-41 | `QID:witness_unit_vector` (`z0Unit`) | WL:150 | `C-SELF{UnitVector index}` + via WL-25 | W1 expr |
| WL-42 | `QID:witness_unit_vector` (`z1Unit`) | WL:151 | `C-SELF{UnitVector index}` + via WL-25 | W1 expr |
| WL-43 | `QID:A0` (`A0lead`) | WL:199 | `C-SELF{v0@WL:197}` + `C-CONVENTION{I}` + `C-FREE{a,omega,cs,M0}` + via WL-31 | W2 expr |
| WL-44 | `QID:A1` (`A1lead`) | WL:200 | `C-SELF{v1@WL:198, literal 3}` + `C-CONVENTION{I}` + `C-FREE{a,omega,cs,D1}` + via WL-32 | W3 expr |
| WL-45 | `QID:expected_A0` | WL:201 | `C-SELF{typed pathA_29 form, uncited}` + `C-CONVENTION{I}` + `C-FREE{a,omega,M0,cs}` + via WL-33 | W2 expr |
| WL-46 | `QID:expected_A1` | WL:202 | `C-SELF{typed form + literal 2, uncited}` + `C-CONVENTION{I}` + `C-FREE{a,omega,D1,cs}` + via WL-34 | W3 expr |
| WL-47 | `QID:A0_residual` (`resA0`) | WL:203 | via WL-43, WL-45 | W2 expr |
| WL-48 | `QID:A1_residual` (`resA1`) | WL:204 | via WL-44, WL-46 | W3 expr |
| WL-49 | `QID:positive_bounded` (`positiveBoundedPair`) | WL:420 | `C-SELF{$Assumptions positivity @WL:26-36}` + `C-FREE{K0c,Keta,TOmega,Z0ret,Z1ret}` | W2/W3 expr |
| WL-50 | `QID:A0_order` | WL:640 | via WL-43 (assertion expression, §7.2) | W2 expr |
| WL-51 | `QID:A1_order` | WL:641 | via WL-44 (assertion expression, §7.2) | W3 expr |
| WL-57 | `QID:witness_constraint_residuals` | WL:552 (`With` local at `:573`, §7.2 rule 2) | via WL-35, WL-41/WL-42 | W1 expr |
| WL-58 | `QID:delta_return` | WL:552 (`With` local at `:574`) | via WL-37, WL-41/WL-42 | W1 expr |

*(WL-17…WL-51 = 35 entries, plus WL-57 and WL-58 = **37 rows** in §D.2. The two trailing IDs keep the
sequential numbering of §D.3/§D.4 stable; they are listed here because they belong to this block.)*

### D.3 Declared literals — `B-DECLARED-LITERAL` (2 rows)

Same treatment as §C.3: `{C-SELF}`, `PF = false`, `A-REDUCED` (E13), `is_tier = no-tier:unclassified-nonfed`,
`SB = tier3-emergent`, `basis = named-route`, `Δ = Y`.

| # | QID | site | value |
|---|---|---|---|
| WL-52 | `QID:v0` | WL:197 | `1` |
| WL-53 | `QID:v1` | WL:198 | `1/2` |

⚠ **WL has no `power1` occurrence** — the `.wl` writes the ω-power as an inline literal `3` inside
`A1lead` (`WL:200`), which is an unnamed leaf under G1. This is a genuine `.py`-only QID, recorded in §G.

### D.4 Constitutive propositions — `B-POSTULATED` (3 rows)

| # | QID | site | axis A | CONV | is_tier | SB | basis | Δ | witness | ev |
|---|---|---|---|---|---|---|---|---|---|---|
| WL-54 | `QID:pathA29_Z_is_premise` | WL:183 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 premise | E8,E14,E15 |
| WL-55 | `QID:pathA29_boundary_dof_none` | WL:183 | `A-IRREDUCIBLE-STRUCTURAL` | false | `tier1-structural` | `tier1-structural` | `none` | N | W1 premise | E8,E14,E15 |
| WL-56 | `QID:time_convention_exp_minus_i_omega_t` | WL:339 (item at `WL:360`) | `A-UNADJUDICATED` | **UNADJUDICATED** | `no-tier:unadjudicated` | `no-tier:convention` | `convention-candidate` | **Y** | W2/W3 premise | E14,E16 |

### D.5 WL tally

16 (D.1) + 37 (D.2) + 2 (D.3) + 3 (D.4) = **58 occurrences**.

**Grand total: 63 + 58 = 121 occurrences.**

---

## E. Evidence (§9), by code

- **E1** — axis B `B-DECLARED-UNASSIGNED`: declaration locus is the row's binding site; **no assignment
  exists in PY**. Verified by reading PY end-to-end: the only re-bindings of these names are new locals
  (`k0 = K0c` `PY:284`, `z0 = …` `PY:287`), never the symbol itself.
  **E1′** — the same for WL: declaration at `WL:26-36`; the only bindings are `Module` locals in `caseFor`
  (`WL:430-450`, `:447-450`), which §7.2 rule 2 excludes from global binding.
- **E2** — axis A `A-UNADJUDICATED` for `omega`. `A-INDEPENDENT-VARIABLE` **costs evidence** (§3.1): §9
  requires loci where the ledger **sweeps, integrates over, or differentiates with respect to** it. Neither
  engine does any of the three: the only ω-directed operation is power extraction
  (`as_powers_dict().get(omega)` `PY:351-352`; `Exponent[…, omega]` `WL:640-641`). ⇒ falls to
  `A-UNADJUDICATED` (§9.0), which keeps the row in tier 1's upper bound. See §J item 1.
- **E3** — axis A `A-REDUCIBLE-UNDERIVED` for `a`. Named route `a = ħ/(m_GNLS c_s0)`, recorded at
  `REG:132`; opened and read — the row reads `` | `a` (pin) | `L` | I-1 (004) | `CONV` |
  `= ħ/(m_GNLS c_s0)` | pin choice — never free ``. Executable within the framework: `ħ`, `m_GNLS` are
  action primitives (`REG:125`, `:126`; `MAP:49`) and `c_s0` is `DERIVED` at `REG:129` / edge R1 `REG:268`.
  **CONV flag `UNADJUDICATED`** (§3.4): a convention claim **is** made (`REG:132` `CONV`), clause (a)'s
  transformation group is not documented and clause (b) is untested — and §3.4's ⭐⭐ rule bites, since no
  dimensionless observable is computable from `a` in this stage (the `A0/A1` magnitude is explicitly
  DEFERRED, `NOTE:94`, `:782-784`). Buys no exclusion; reported in `convention-unadjudicated`.
- **E4** — axis A `A-UNADJUDICATED` for `c_s`. A route exists for `c_s0` (`c_s0 = √(5Kρ0⁴/m_GNLS)`,
  `REG:129`, edge R1 `REG:268`; `c_s² = 5Kρ⁴/m` at `MAP:62`) but **whether stage023's `c_s` is that
  quantity is explicitly unresolved in two independent places**: `REG:129` carries
  "⚠ **STAGE023 WORK POINTER — no adjudication** … tracks the unresolved stage023 evaluation point", and
  `NOTE:513-518` records `c_s` as "**UNDETERMINED** against 018's bulk asymptotic `c_s0`, because nothing in
  the four places that define this carrier states which `ρ` it is evaluated at". §8.1/§8.3 forbid the builder
  adjudicating that identity ⇒ the route cannot be attached to this occurrence. `should_be` carries it as
  `named-route` (`REG:268`).
- **E5** — axis A `A-UNADJUDICATED` for `M0`, `D1`. A route is named — `M0(omega)=int_brane S_leak(omega,x)
  d^3x` at `S008:315` (opened: it is a `print`, not a computation; `S008:136` declares `M0 D1 Q2 R0 R1` as
  free symbols and `S008:546` states "No dimensions are invented here for M0/D1/Q2; they remain frozen-slice
  normalized source moments"), and `REG:83` classes the moments as "constructions not parameters". **Whether
  that route is executable within this framework cannot be decided from the evidence in reach** — it runs
  through the leak-source profile `S_leak`, whose availability is not settled anywhere opened. §3.1.2 states
  no undecidable branch (see §J item 2); §3.1.1's third branch is applied by analogy ⇒ `A-UNADJUDICATED`.
- **E6** — axis A `A-UNADJUDICATED` for `D0`. **Two competing named routes, and the identity that would pick
  one is recorded as UNDETERMINED.** Route (i): `D0 = K₂ − (B̃0+Z̃0)`, opened at `S017:31`, exported to
  022/023 at `S017:113-115`, with `REG:185` carrying `D0=K̃+6T̃_Ω−(B̃0+Z̃0)`. Route (ii): `D0 := K_* − Q/Δ`
  at `PREREG:195`, applied only after the closed solve. The corpus itself records the 017 identification as
  unresolved: `REG:185` carries "⚠ **STAGE023 WORK POINTER — no adjudication** … tracks the contradictory
  `D0` closure", and `NOTE:460-476` states "**As the identification with stage017's exported D-lane it is
  `UNDETERMINED`**". ⇒ `A-UNADJUDICATED`; `should_be = tier1-debt`, basis `named-route` at `S017:31`.
- **E7** — axis A `A-REDUCIBLE-UNDERIVED` for `K0c`, `K_eta`, `T_Omega`. Named route at `REG:170`, opened
  and read: "an explicit scalar-reduction (profile+measure) of the ℓ=0 Gate-2 collective (§8.2, stage013) +
  the ℓ=1 §9.4 harmonic (stage017) sectors to the calibrated wall packet would DISCHARGE this (collapse them
  into `{μ_η,T_w,β,T_Ω}`) — `PENDING`". **Executable within this framework** (§3.1.2): every piece it needs
  is already in the ledger — 013's Gate-2 collective, 017's ℓ=1 harmonic sector, and the calibrated wall
  packet `{μ_η,T_w,β,T_Ω}` (`REG:174`, `:176`, `:182`, `:184`) — what is missing is that nobody wrote it
  ("no equation yet"). The register itself calls it "**COUNTED as reduction debt**".
- **E8** — axis A `A-IRREDUCIBLE-STRUCTURAL` for `Z0_ret`, `Z1_ret` and the pathA_29 premise. Route **is**
  named: `REG:169` — "named route: the track-3/Gate-6 **nonlinear** return would derive the transmissions
  from the medium … `PENDING`". §3.1.2 therefore requires the **framework extension** to be named, and it is:
  the **nonlinear brane↔bulk return closure (Gate 6)**, which this ledger does not model —
  `REG:309` (R42): "Gate-6's nonlinear closure (**sim-deferred**) must supply 2 independent return equations
  (or an equivalent return law)"; `REG:168` (the `d` slab row) records "the slab family is POSTULATED; the
  real return geometry = the deferred nonlinear brane↔bulk closure"; `NOTE:8` "Gate 6 (the nonlinear
  return-selector), is sim-deferred". The medium's defining properties are unchanged by that extension ⇒
  STRUCTURAL, not POSTULATE (§3.1.1).
- **E9** — axis A `A-IRREDUCIBLE-STRUCTURAL` for `Omega_U`, `Omega_W`, `R_mix`, `g_U`, `g_W`. Route named at
  `PREREG:186-195` (opened): the mixed-Schur port identities `Δ := Ω_U²Ω_W² − R_mix²`,
  `P := Ω_U²G_W + R_mix G_U`, `Λ := P/Δ`, `N0 = Λ²`, `P0 = N0/D0` — which is exactly `PY:222-225` /
  `WL:107-110`. The **framework extension** the route requires is named in the same passage: the formulas are
  "applied **AFTER the closed solve**" and "evaluated on the **derived self-consistent background**"
  (`PREREG:186-188`), the background being `R0(w)` from the static self-consistent throat balance
  `PREREG:178-184`, which the same page states "requires the **full closed matter/gauge solve** (it is
  realized numerically by the closed solve, not as a posited formula)" (`PREREG:180`). That solve is the
  ledger's standing deferred interior — `REG:300` (R33): "the same deferred nonlinear-throat interior as R30".
  ⚠ **No `REG:` row exists for any of these five** (§I).
- **E10** — axis B `B-EXECUTED`: the computation locus is the row's binding site; the input-leaf loci are the
  row's "reduces to" column. Both engines run to exit 0 (`NOTE:23-24`, dual-engine tallies); no re-run was
  performed by this census (§11: the census does not re-run a stage).
- **E11** — axis A `A-REDUCED`: the reduction is performed at the row's own binding site (the expression
  shown), and the loci of the quantities it reduces TO are the "reduces to" column.
- **E12** — axis C `PF = C-UNRESOLVED`: the walk completed and terminated on at least one `C-FREE` leaf with
  **no** `C-FIELDEQ`, `C-EXTERNAL` or `C-PEER` leaf anywhere in the closure. Per §3.3 that is
  `C-UNRESOLVED`, ⛔ not `PHYSICS-FED = false`. The two candidate physics leaves were tested and both fail:
  (i) **the operator**: §3.3.1(1) rules `rank()`/`NullSpace`/`MatrixRank` out as leaves ("a generic
  linear-algebra routine is NOT an operator"), and what *could* be `C-FIELDEQ` — the three assembled
  constraint rows — carries no cited equation locus in either engine (`PY:220-232`, `PY:384`; the `.wl`
  comment `WL:104` names "Stage017 grouped-P2 port kernel", a **stage name, not a locus**, which §3.3 sends
  to `C-SELF`); (ii) **transcription**: the pathA_29 form at `PY:344-348` / `WL:201-202` is likewise
  attributed by stage name only (`PY:344` comment; `PY:1050`; `NOTE:794`), so it is `C-SELF`. ⚠ The source
  *does* exist and was opened — `S009:194-195` (`neighbor_fraction(ε) = ε/(1+ε)`) and `S009:261`
  (`z_formula = −M0·ε0/(1+ε0)`) — so this is a **citation defect, not a fabrication**; recorded in §J item 4.
- **E13** — axis A `A-REDUCED` for `v0`, `v1`, `power1`. **The reduction is performed** at `S022:209-210`
  (opened): `radiative_power = 2*lval+1`; `radiative_coeff = compact(y_series.coeff(z, radiative_power)/sp.I)`
  — the coefficient is extracted from the outgoing-Hankel DtN series built at `S022:153-160`, and checked
  against the typed target at `S022:403-406` ("`ell={ell} derived radiative_coeff matches typed cross-ell
  target`"). It reduces TO the spherical-Bessel/Hankel kernel (`S022:153-160`) and the ℓ index. ⚠ The typed
  table `S022:136-140` is stage022's *target*, not its derivation; the derivation is `:209-210`.
- **E14** — axis B `B-POSTULATED`: where the model posits it — `PY:425` / `WL:183` for the pathA_29 premise
  (the artifact itself names it "**the keystone return premise**" at `PY:1051` / `WL:767`), `PY:612-615` /
  `WL:360` for the time convention (`computed_class: "convention"`, tag `exp_minus_i_omega_t`). ⛔ No code
  locus owed (§3.2).
- **E15** — route-2 bounds for the pathA_29 premise. **(i) depended on by a reported result**:
  `RES:023:RETURN-UNDERDETERMINED` — `NOTE:59-61` "**Why the return survives (the keystone premise):**
  pathA_29 records `Z_is_premise = True`, `boundary_dof = none` — the brane↔bulk return admittance is a
  *premise* of the flat-slab family, not something the linear theory supplies an equation to fix";
  `NOTE:796-797`. **(ii) not derived elsewhere in the ledger**: `REG:169` carries the transmissions as
  `FREE-UNREDUCED` with the derivation `PENDING`; `REG:168` records the slab family itself as POSTULATED
  ("the slab family is POSTULATED; the real return geometry = the deferred nonlinear brane↔bulk closure").
- **E16** — route-2 bounds for the time convention. **(i)** `RES:023:A0-FORM` / `A1-FORM` claim the residual's
  **sign**, and the `i` prefactor (`sp.I`, `PY:342-348`; `I`, `WL:199-202`) data-depends on the convention;
  `NOTE:409` derives `[ω]` from "Time convention `e^{−iωt}`". **(ii)** nothing in the ledger derives a time
  arrow; it is declared. CONV flag `UNADJUDICATED` (§3.4): a convention claim is made (`computed_class:
  "convention"`), but no transformation group is documented (clause a) and no invariance is demonstrated over
  the model's stated observables (clause b) — and the one observable in reach, the claimed sign, is the very
  thing the convention fixes. Buys no exclusion; reported in `convention-unadjudicated`.

---

## F. Out of scope (§5.8 output 15) — with reasons and counts

Counted at **named-object granularity** (one entry per named binding the artifact introduces), from the loci
listed. ⛔ Nothing is silently truncated; where the granularity is coarser than "every temporary", it is said.

### F.1 `reached-by-no-reported-result` (§5.8 output 16) — **211**

⭐ This is the finding-bearing exclusion and it carries its own count.

| group | loci | count |
|---|---|---|
| PY `SOURCED_DIMS` entries | PY:461-483 | 22 |
| PY `EXPECTED_DIMS` entries | PY:487-493 | 7 |
| PY walked `computed` dims | PY:525 (7 keys) | 7 |
| PY `dimension_records` | PY:542-570 | 29 |
| PY dimension-basis objects | PY:179,180,181 | 3 |
| PY dimension-check status objects | PY:520,528,534,535 | 4 |
| PY `R0`, `R1` (dimension-only, unread — the stage's own **U5**, `NOTE:638-639`) | PY:168 | 2 |
| PY `P0_physical` | PY:231 | 1 |
| PY `forward_relations_ok` | PY:310 | 1 |
| PY verdict-condition scalars (`overcancel`, `no_consistent_return`, `A0_form`, `A1_form`) | PY:317,318,362,363 | 4 |
| PY provenance items other than `time_convention` (+ `ok`, `baseline_emits_…`) | PY:576-611,636,637 | 10 |
| PY decoupling detector outputs | PY:646-659 | 7 |
| PY `conditions` / `verdict` / `verdict_read_set` | PY:703-718 | 9 |
| **PY subtotal** | | **106** |
| WL `baseDims` entries | WL:248-254 | 22 |
| WL `expectedDims` entries | WL:258-261 | 7 |
| WL `dimensionAudit` computed + status | WL:270-282 | 10 |
| WL `emitDimensionRecords` records | WL:298-332 | 29 |
| WL dimension helpers (`zeroDim`, `dimScale`, `dimensionAxisSlots`, `dimensionAxesLabel`, `dimensionComponents`, `printDimRecord`, `dimOf`) | WL:209-222 | 7 |
| WL `baselineDimAudit` | WL:286 | 1 |
| WL `P0Physical` | WL:105 (entry `:116`) | 1 |
| WL `caseFor["baseline"]` result keys (consumed only by the verdict ladder and the control battery) | WL:424 (keys `:485-492`) | 16 |
| WL provenance items other than `time_convention` (+ `Ok`, `NoConcreteEpsilonMagnitude`) | WL:342-359,396,397 | 10 |
| WL `failurePriority`, `verdictFromFlags` | WL:401,411 | 2 |
| **WL subtotal** | | **105** |

⛔ **The `reached-by-no-reported-result` set is larger than the in-universe set (211 vs 121), and it contains
the entire dimensional block of both engines.** ⭐ It contains, specifically, the very construction
`CENSUS_SCHEMA.md` §3.5 cites as the reason axis C exists (`dim_of`'s `symbol_dims[expr]` leaf lookup at
`PY:436-439` over `dims = dict(SOURCED_DIMS)` at `PY:505`). Under §7.1.1's universe rule that construction is
**out of the universe**, because no reported physical result's closure reaches it. Reported, not resolved —
see §J item 3.

### F.2 Controls and deliberate negatives — **≈ 78** (PY ≈ 52, WL ≈ 26)

Enumerated by group: PY `Mutation` fields (`PY:187-202`, 13); the counterfactual selector machinery
(`PY:236`, `:240-242`, `:246-266`, `:270` — 11 named objects); the four control symbols `q_free`, `gain0`,
`gain1`, `eta_null` (`PY:172,173,174`, 4); the mutation branches inside the builders (`PY:221`, `:289-300`,
`:306-308`, `:338`, `:339`, `:507`, `:509` — 12); the ablation harnesses (`PY:734-756`, `:759-767`,
`:770-793` — 12). WL: `selectorResiduals`/`selectorConstraints`/`Jselector`/`selectorBasis`/
`selectorReturnImage`/`selectorRank`/`selectorNativeNullity`/`selectorReturnMovingNullity`/`selectorMetadata`
(`WL:153-181`, 9); `oneSelector*` (`WL:185-187`, 3); `inject*`/`etaUnit` (`WL:189-194`, 5);
`corruptV1A1`/`corruptV1Residual`/`corruptOrderA1` (`WL:205-207`, 3); `corruptSourcedDims`/`corruptFreeDims`
+ their audits (`WL:289-296`, 4); `dynamicAblation`/`checkFailureAblation`/`portDependency` (`WL:496-550`, 3).
Backed by the substrate's own labels: `derived_from_named_pde = False`, `control_only = True`
(`PY:256-257`), `NOTE:72-75`, `NOTE:774-777`.

### F.3 Test scaffolding and display constants — **≈ 40** (PY ≈ 25, WL ≈ 15)

Verdict token strings (`PY:27-35`, 9; `WL:16-24`, 9), pass/fail counters (`PY:24-25`; `WL:12-14`), and the
assert/format helpers (`PY:46-160`, `WL:38-102`).

### F.4 `derived-proposition-not-constitutive`, `retired-or-excluded`, `paused-or-pending` — **0 each**

No assertion in either engine failed route 2's second bound; no stage023 quantity appears in a retired
register row (`REG:139`, `:159` — checked); no stage023 quantity is recorded PAUSED or PENDING-as-an-approach
(the `PENDING` at `REG:170` is a **reduction debt**, which §5 counts inside tier 1, not an exclusion).
⭐ **Provisional-exclusion sub-count = 0** (§5.8 output 15's required sub-count).

---

## G. The required output set (§5.8) — `is_tier` counts only (§6.1)

⛔ **Reported at both levels (§10.1).** Occurrence level is exact. QID level is given **twice**: under the
PROPOSED cross-engine merges, and unmerged (since §8.3 says an unadjudicated merge is not a merge).

**QID inventory (proposed):** 67 distinct quantities/propositions. Unmerged: 121.

| # | output (§5.8) | occurrence level | QID level (proposed merges) | QID level (unmerged) |
|---|---|---:|---:|---:|
| 1 | **TIER 1 — range** | **[26, 116]** | **[13, 64]** | **[26, 116]** |
| 1a | ↳ `tier1-debt` | 8 | 4 | 8 |
| 1b | ↳ `tier1-structural` | 18 | 9 | 18 |
| 1c | ↳ `tier1-postulate` | **0** | 0 | 0 |
| 2 | **TIER 2 — calibrated** | 0 | 0 | 0 |
| 3 | **TIER 3 — emergent** | 0 | 0 | 0 |
| 3a | ↳ split by axis B (§5.4) | 0 `B-EXECUTED` / 0 unexecuted | — | — |
| 3b | ↳ split by calibration-propagation (§5.5) | 0 `tier3-calibration-propagated` / 0 `tier3-held-out` | — | — |
| 4 | **`DERIVED`** (§4) | **0** | 0 | 0 |
| 5 | near-miss `executed-but-not-physics-fed` | **0** | 0 | 0 |
| 6 | near-miss `derived-in-form-but-unexecuted` | 0 | 0 | 0 |
| 7 | near-miss `physics-fed-but-declared-literal` | **0** | 0 | 0 |
| 8 | `unclassified-nonfed` | 5 | 3 | 5 |
| 9 | convention bucket (`no-tier:convention`, flag `true`) | **0** | 0 | 0 |
| 10 | `unadjudicated` | 90 | 51 | 90 |
| 11 | `convention-unadjudicated` | 4 | 2 | 4 |
| 12 | `self-referential` | **0** | 0 | 0 |
| 13 | conflict set | 6 (all **intra-occurrence**) | 3 | 6 |
| 14 | stage043 delta (three directions) | §I | §I | §I |
| 15 | out-of-scope list | §F (211 + ≈78 + ≈40; provisional sub-count **0**) | — | — |
| 16 | `reached-by-no-reported-result` | **211** | — | — |
| 17 | reported-result set | §A (3 results) | — | — |
| 18 | `unvalued-in-universe` (Route C) | n/a (QID-level, corpus-wide) | **not takeable from one stage** | — |
| 19 | `C-PEER` populations, **reported separately** | `peer-cited-in-artifact` = **0**; `peer-cited-in-stage-note` = **0** | 0 / 0 | 0 / 0 |
| 20 | `no-tier:independent-variable` | **0** | 0 | 0 |

**Arithmetic check (occurrence level):** 26 (tier 1) + 5 (`unclassified-nonfed`) + 90 (`unadjudicated`)
= 121 ✓. Tier-1 upper bound = 26 + 90 = 116 (§5.6).
**Arithmetic check (QID, proposed):** 13 + 3 + 51 = 67 ✓; upper bound 13 + 51 = 64 (§10.2.2 rule 4: no QID
already in the lower bound is re-added by an unadjudicated occurrence — none applies here, since no QID mixes
a tiered with an untiered occurrence).

⛔ **No fraction is quoted**, because §10.1's rank denominator is `OPEN-METHOD` (§13) and every fraction would
need both denominators against a rank that is not yet computable. The three row counts above are the audit
receipts.

**Why buckets 5, 7, 9, 12 and 20 are zero — stated so the zeros are legible as measurements:**
- **5** — every `B-EXECUTED` row is `C-UNRESOLVED`, and §3.3 forbids asserting `C-UNRESOLVED` into
  `executed-but-not-physics-fed` (a positive claim). The rows that *would* populate it are in §F.1.
- **7** — no row is `PHYSICS-FED`, so near-miss 3 is unreachable. The three literals in §C.3 are one
  citation locus away from it (E12, E13).
- **9** — no occurrence satisfies §3.4's clauses (a) **and** (b); the two candidates are `UNADJUDICATED`.
- **12** — no closure leaves an occurrence and returns to it through another `(artifact, binding-site,
  quantity)` triple. Checked specifically for the `T`/`ε` pair (`PY:302-305`), which is the shape most likely
  to round-trip: `T0 ← K0,Z0` and `eps0 ← K0,Z0` are siblings, not a cycle. The stage's own
  `forward_relations_ok` identity (`PY:310-315`) *is* an `X≡X`-shaped check, but it is out of universe (§F.1)
  and the stage already de-counts it (`NOTE:765`).
- **20** — see E2 and §J item 1.

---

## H. Conflicts (§10.3) — **6, all intra-occurrence**

⛔ Not resolved, not adjudicated. Both claims and both loci carried on the row.

| rows | QID | claim A (locus) | claim B (locus) |
|---|---|---|---|
| PY-07, WL-07 | `K0c` | `FREE-UNREDUCED`, **PENDING** scalar-reduction, **counted as reduction debt**, "NOT identified with the raw 013/017 densities" — `REG:170` | the identification stated **as performed**: "`K0c` ← stage013's Gate-2 collective (a,L) reduction (§8.2)" and "likely **DERIVED** manifestations, NOT new counted knobs" — `SRCMAP:250,253` |
| PY-08, WL-08 | `K_eta` | same (`REG:170`) | "`K_eta` = R29-derived `T_wβ²` (013)" — `SRCMAP:251` |
| PY-09, WL-09 | `T_Omega` | same (`REG:170`) | "`T_Omega` = 017's counted `T_Ω`" — `SRCMAP:251-252` |

Both loci opened and read. The census adjudicates its own axis values from evidence and does **not** inherit
either class (§2); the substrate's disagreement is recorded as a finding about the substrate. `NOTE:730-734`
(W4) already flags this and calls the source map "stale pre-build" — ⛔ that characterisation is the
substrate's, not a resolution, and §10.3 forbids picking.

**Inter-occurrence conflicts: 0.** No QID's occurrences disagree in `is_tier` — the two engines agree on every
shared quantity. **`mixed-adjudication` QIDs: 0.**

⚠ **Recorded but deliberately NOT entered in the conflict set:** the `D0` / `P0_raw` / `P0_physical`
contradiction (`REG:185` `M T⁻²` vs `PY:468` / `WL:250` `M L⁻¹T⁻²`; `NOTE:454-476`, W1/U8). It is a
**dimensional** CORRECT/UNDETERMINED disagreement, which §11 places on a different axis from provenance, and
the census does not rule on dimensional correctness. It is carried instead as the evidence that blocks
`D0`'s axis A (E6).

---

## I. stage043 reconciliation (§8.4) — the stage023 slice only

⛔ **Not the full 152-ID reconciliation**: §8.4 is a corpus-wide obligation and cannot be discharged from one
stage. What is takeable here is the slice of `S043`'s manifest that touches a stage023 census QID. All IDs
below were read from `S043:216-300`.

**Direction 1 — census QIDs with no `REG:` ID (the census's genuine extension): 60 of 67 proposed QIDs.**
Named among the physically substantive ones: `Omega_U`, `Omega_W`, `R_mix`, `g_U`, `g_W`, `D0`, `M0`, `D1`,
`c_s` (stage023's), `P0_raw`, `P_port`, `Delta_port`, `N0_from_port`, `v0`, `v1`, `power1`, `A0`, `A1`,
`expected_A0`, `expected_A1`, and every rank/witness object. ⭐ **All five port-kernel symbols carrying
`tier1-structural` in this census have no `REG:` ID at all** — verified by enumerating every `REG:` literal
in `S043` and searching it.

**Direction 2 — `REG:` IDs mapping to ≥ 1 census occurrence: 7** (covering 7 of the 67 proposed QIDs).

| `REG:` ID (`S043` locus) | stage043 category | census rows | census `is_tier` |
|---|---|---|---|
| `REG:b:K0c` (`S043:225`) | `FREE-UNREDUCED-DEBT` / `CAT_DEBT` | PY-07, WL-07 | `tier1-debt` |
| `REG:b:K1_eta_direction` (`S043:225`) | `FREE-UNREDUCED-DEBT` | PY-08, WL-08 | `tier1-debt` |
| `REG:b:K1_TOmega_direction` (`S043:226`) | `FREE-UNREDUCED-DEBT` | PY-09, WL-09 | `tier1-debt` |
| `REG:derived:Z0_ret_alias` (`S043:297`) | not-counted "derived" bulk | PY-10, WL-10 | `tier1-structural` |
| `REG:derived:Z1_ret_alias` (`S043:297`) | not-counted "derived" bulk | PY-11, WL-11 | `tier1-structural` |
| `REG:conv:a_pin` (`S043:335`) | `"CONV"` / `CAT_CONVENTION` (`S043:336`) | PY-02, WL-02 | `tier1-debt` (flag `UNADJUDICATED`) |
| `REG:derived:K1_sum` (`S043:298`) | not-counted "derived" bulk | PY-23, WL-24 | `no-tier:unadjudicated` |

⚠ **Two divergences, recorded not resolved.** (a) stage043 places `Z0_ret`/`Z1_ret` in the **not-counted
`derived` bulk** as aliases, with the counted content sitting in `REG:b:epsilon0`/`epsilon1` (`S043:224`);
this census places stage023's occurrences of them in `tier1-structural`. The two are not in contradiction
about the *dof count* — `REG:169` says the aliasing adds no new dof — but they are different dispositions of
the same symbol, and §2 forbids inheriting either. (b) stage043 classes `a` as `conv`; this census leaves the
convention flag `UNADJUDICATED` (E3), which buys no exclusion.

**Direction 3 — Route C `unvalued-in-universe`: not takeable.** Route C requires establishing that **no
artifact in the ledger** values the quantity — a corpus-wide finding. From this stage alone I can record that
16 quantities are `B-DECLARED-UNASSIGNED` **in both engines of stage023**; whether any is valued elsewhere is
outside this pass. ⛔ Reported as owed, not as zero.

---

## J. Blocked rows, and every place the schema did not decide

1. **BLOCKED — `should_be_basis` has no value for an independent-variable expectation.** Rows PY-01, WL-01
   (`omega`). §3.1's `A-INDEPENDENT-VARIABLE` costs evidence and §9 admits exactly three uses — swept,
   integrated over, differentiated with respect to. Neither engine does any of the three; the only
   ω-directed operation is **power/series-order extraction** (`PY:351-352`, `WL:640-641`), which the
   three-verb list does not cover. `is_tier` is forced to `no-tier:unadjudicated` (correct and conservative),
   but there is then **no admissible `should_be_basis`** (§6's four values are `named-route`,
   `physical-picture-expectation`, `convention-candidate`, `none`) by which to record that the picture reads
   ω as a coordinate. §6 forbids a departing `should_be_tier` on basis `none`, so the rows are recorded with
   `should_be_tier = is_tier`, `basis = none`, and the expectation is recorded **here** instead.
   **Missing rules:** (i) §9's independent-variable evidence list should say whether series-order extraction
   qualifies; (ii) §6 needs a basis value for an independent-variable expectation.
2. **§3.1.2 has no undecidable branch.** §3.1.1 gives three branches for STRUCTURAL-vs-POSTULATE, the third
   being `A-UNADJUDICATED`. §3.1.2 — the debt-vs-structural test, which is the split the tier-1 deliverable is
   *about* — gives only two. Rows PY-04/05, WL-04/05 (`M0`, `D1`) sit exactly there: a route is named
   (`S008:315`) and its executability cannot be established. §3.1.1's third branch was applied by analogy
   (⇒ `A-UNADJUDICATED`, keeping the rows in tier 1's upper bound), and the borrowing is recorded here rather
   than hidden. **Missing rule:** an explicit undecidable branch in §3.1.2.
3. **§3.5's motivating construction is out of the universe.** §3.5 names `PY:436-439` / `PY:505` as the
   corpus instance of the shape axis C exists to detect. §7.1.1's universe rule does not reach it, because no
   reported physical result's closure passes through the dimension block (§F.1). The consequence is
   mechanical, not editorial: `executed-but-not-physics-fed` = 0 for this stage. ⛔ Not a defect I resolved —
   the two sections are consistent and simply do not meet on this stage. **What a reviewer must decide:**
   whether a stage's own dimensional gate is intended to be censused at all under §7.1.1.
4. **`should_be_basis` has no value naming a citation defect.** All 78 `B-EXECUTED` rows are blocked from tier 3
   solely by axis C, and the blocker is that no engine and no stage note carries a **file-and-line** for the
   equation or the imported value (E12). §6's `named-route` means "actionable debt" (a reduction to run) and
   `physical-picture-expectation` means "a hunch, no route named"; neither says "the route exists, is
   already executed elsewhere, and only the locus is missing". `physical-picture-expectation` was used as the
   closest admissible value on those rows and on PY-04/05, and `named-route` on §C.3's literals (where the
   *reduction* is what the route is). **Missing rule:** a `should_be_basis` for a citation-blocked axis C.
5. **§7.1.1's ban on verdict-token entry points versus REPORTED_RESULTS' cited loci.** REPORTED_RESULTS
   §2.1 lists `PY:862` and `PY:1036` — both `verdict_residual(...)` asserts — as loci of
   `RES:023:RETURN-UNDERDETERMINED`, which §7.1.1 forbids as walk entry points. Resolved here by walking only
   the computed loci and using the verdict asserts as test-(a) evidence, exactly as §7.1.1 permits. ⚠ It is
   the single decision that most changes this stage's universe (§A), and REPORTED_RESULTS is not amendable
   (§7.1.1) — so if a reviewer reads those loci as entry points, the universe grows by most of §F.1 and this
   pass must be re-taken, not patched.
6. **`SELF-REFERENTIAL` on a `C-UNRESOLVED` closure — an ambiguity, adjudicated.** §3.3's parenthetical says
   "a closure containing an **untraceable** leaf is `C-UNRESOLVED` as a whole, and `SELF-REFERENTIAL` is then
   undefined". A `C-FREE` leaf is *not* untraceable — the walk completed and terminated at a known
   declaration — so this census reports `SELF-REFERENTIAL = false` (determined) on those rows rather than
   `undefined`. ⚠ If the reviewer reads the parenthetical as covering `C-FREE` too, then ~110 rows carry
   `SELF-REFERENTIAL = undefined` and the `self-referential` bucket becomes unmeasurable rather than 0.
7. **Locus-abbreviation deviation** from §7.2 — declared in §0.
8. **Route C (§8.4) is not takeable from one stage** — §I direction 3, recorded as owed.

**Blocked count: 2 rows** (PY-01, WL-01) blocked on one field each; **0 rows blocked on an axis value**;
**0 rows unclassifiable**. Every one of the 121 rows carries all three axes, the flag, both tiers, a basis, a
delta, LIVE/RETIRED, IDs and its §9 evidence.

**Unadjudicated count and its causes (the 90 `no-tier:unadjudicated` occurrences):**
| cause | rows |
|---:|---|
| 78 | `A-REDUCED` computed expressions failing tier 3's `PHYSICS-FED` conjunct on a `C-UNRESOLVED` closure (§5.7's last line) — no cited field-equation or peer locus anywhere in either engine (E12) |
| 4 | identity to the quantity that carries the route is recorded UNDETERMINED by the corpus itself — `c_s` ×2 (E4), `D0` ×2 (E6) |
| 4 | route named, executability undecidable and §3.1.2 offers no branch — `M0` ×2, `D1` ×2 (E5, §J item 2) |
| 2 | `A-INDEPENDENT-VARIABLE` unevidenced under §9's three-verb list — `omega` ×2 (E2, §J item 1) |
| 2 | bare convention with no axis-A route — `time_convention` ×2 (E16) |

---

## K. Coverage statement

**Selection rule.** This pass classifies **every occurrence in the transitive input closure of the three
reported physical results of stage023, in both engines** — entered at the value-bearing loci of §A and never
at a verdict token — under the granularity rules G1–G5 of §B. That is **121 occurrences / 67 proposed QIDs**,
and it is complete for that rule: no in-closure binding site was skipped.

**What is NOT covered, with counts:**

1. **Out-of-scope occurrences in the two engines — 211 `reached-by-no-reported-result` (exact, enumerated in
   §F.1) plus ≈ 78 controls/deliberate negatives and ≈ 40 test-scaffolding/display objects (§F.2, §F.3,
   counted at named-object granularity, by group, with loci).** The two ≈ figures are group counts, not
   exhaustive per-temporary enumerations; that is stated rather than presented as exact.
2. **The stage note as an artifact.** §7.2 admits a stage note as an artifact and "one markdown table row is
   one binding site"; the directive scopes this pass to the two engines. `NOTE` is therefore used **only** as
   an evidence locus, and its own binding sites are **not classified** — most visibly the 34-row table at
   `NOTE:405-440`, whose rows are dimension records and would in any case be `reached-by-no-reported-result`.
   **Uncounted: the stage note's binding sites (≥ 34 in that table alone).**
3. **Corpus-wide obligations that a single-stage pass cannot discharge:** the full 152-ID stage043
   reconciliation (only the 7-ID stage023 slice is given, §I), Route C `unvalued-in-universe` (§I direction
   3), and §10.1's rank denominator (`OPEN-METHOD`, §13).
4. **Identity adjudication.** Every cross-engine merge is PROPOSED (§8.3); §G reports both the merged and
   unmerged QID levels so no count depends on a merge this leg is not permitted to make.

⛔ Nothing was truncated silently, and no count above was steered toward the stage's `24/0/10` or `27/0/7`
dimensional tallies — those are dimensional-correctness verdicts on a different axis (§11) and were used only
as a cross-check that the row universe of §C/§D contains the objects the stage itself treats as real.
