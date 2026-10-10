# Lookups — polarization survey review, round 3 (generated 2026-10-09 19:52)

Generator: `_scratch/polarization/gen/survey_r3_lookups.sh` (sha256sum/grep/sed/cat only).

```
$ (cd _scratch/polarization && sha256sum -c survey_review_baseline_r3.sha256)
POLARIZATION_SURVEY.md: OK
```

```
$ grep -n -o 'Verdict:\*\* [A-Z][A-Z ]*\|Verdict: [A-Z][A-Z ]*\|\*\*Verdict: [A-Z][A-Z ]*' _scratch/polarization/survey_review_r3_claude.txt _scratch/polarization/survey_review_r3_grok.txt
_scratch/polarization/survey_review_r3_claude.txt:1:**Verdict: NEEDS REVISION
_scratch/polarization/survey_review_r3_grok.txt:1:**Verdict: CLEAR
```

## Finding 1: E01's class
```
$ grep -n '| E01 ' _scratch/polarization/POLARIZATION_SURVEY.md | sed -n '$p'
214:| E01 | **Reproduced** | **Conditions:** classify the two-dimensional transverse branch content, under the selected homogeneous, isotropic-inertia in-plane \(D=3\) action, positive coefficients, nonzero wavevector and the stated generic strata; separation from out-of-plane fields is supplied. **Rule 1:** **S10_two_transverse_photons.md:73–76** computes “the nonzero root has D − 1 transverse null directions” and hence two at D=3 [Q5](#q5). Its action, field content and dimensional input do **not themselves state the two-direction result**; supplying D=3 is not supplying D−1=2. Thu
```

```
$ grep -n 'R-S1-01' _scratch/polarization/POLARIZATION_SURVEY.md | head -4
157:| Dimension provenance | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:128–145**, R-S1-01: Dbrane as a derived quantity or “explicitly re-affirmed postulate”, status **OPEN**; “Dbrane = 3 went in and D−1 = 2 came out”. [R28](#r28) | **OPEN**, target S6. The record does not prohibit deriving a codimension-one brane from a four-dimensional bulk elsewhere. |
214:| E01 | **Reproduced** | **Conditions:** classify the two-dimensional transverse branch content, under the selected homogeneous, isotropic-inertia in-plane \(D=3\) action, positive coefficients, nonzero wavevector and the stated generic strata; separation from out-of-plane fields is supplied. **Rule 1:** **S10_two_transverse_photons.md:73–76** computes “the nonzero root has D − 1 transverse null directions” and hence two at D=3 [Q5](#q5). Its action, field content and dimensional input do **not themselves state the two-direction result**; supplying D=3 is not supplying D−1=2. Thu
1275:   390	⚠ **`D = 3` is NOT selected here** ⇒ `R-S1-01`, target **S6**, status OPEN.
1798:   128	### R-S1-01 — the brane's spatial dimension
```

```
$ sed -n 134,136p research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
- **on failure** — S10's headline reads *"light having exactly two polarisations is a statement that our
  space is three-dimensional."* Read backwards it says the opposite: `D_brane = 3` went in and `D−1 = 2`
  came out. Without a delivered `D_brane` the sentence is an assumption restated, ⛔ not a result.
```

```
$ sed -n 374,376p research/pde_ledger_v3/V3_STEP_PLAN.md
⭐ **What S10 measured:** the conditional map `D ↦ D − 1` transverse null directions **for the cases it
swept** (`D = 2,3,4,5`), so the `D = 3` member has two — ⛔ **conditional on the supplied action, the
supplied `[u]`, and BOTH structural premises.**
```

```
$ sed -n 182,183p research/pde_ledger_v3/steps/S10_two_transverse_photons.md
- a D-component in-plane displacement u, with its separation from every
  out-of-plane field inherited rather than tested;
```

```
$ sed -n 73,76p research/pde_ledger_v3/steps/S10_two_transverse_photons.md
> For the supplied curl-only in-plane action, at nonzero wavevector and away
> from allowed exceptional strata, the nonzero root has D − 1 transverse null
> directions in every measured MAIN case D = 2, 3, 4, 5. Thus the D = 3 member
> has two transverse directions.
```

```
$ grep -n 'C1' _scratch/polarization/survey_r0_review_disposition.md | head -3
23:| C1 / G1 | E01 is classed "Reproduced", but `D = 3` goes in and `D − 1 = 2` comes out (R-S1-01, OPEN). E01 also omits the longitudinal branch that S10 and S11 record as a departure from Maxwell, conditional on matter coupling. The survey added "Apparent conflict requires an experimentally comparable model prediction", which the task did not ask for. (both) | **ACCEPT.** Lookups: survey L183 (E01) and L179 (the added rule); prompt L53; S10 L172–176 and L249–251; S11 L15 and L173; register L128–136. | The task's class definitions apply as written. "Reproduced" is used only where a st
24:| C2 | E10, E11 and E13–E15 are held to a different standard from E01. The plan states that the two polarizations are degenerate by symmetry in a homogeneous brane. (Claude) | **ACCEPT.** Lookups: plan L1183; survey L192–197. The leg's script also shows the anisotropic control breaks the degeneracy (reviewer evidence). | C1's rule is applied to these rows too, with the conditions (homogeneity, isotropy) stated. |
```

```
$ grep -n 'Reproduced\|supplied input that states it\|E01 and E22' _scratch/polarization/survey_repair3_fresh_author_prompt.md
17:1. **Reproduced.** A v3 step computes the measured fact from inputs that do not themselves state it.
26:   constrains. **A fact that follows only from a supplied input that states it belongs here, and this takes
40:4. **E01 and E22.** Their stated reasons follow the definitions above, including rule 4's precedence, and do not
```

## Finding 2: where the parity condition is located
```
$ grep -n -c 'S11bB_interface_assembly.md:53' _scratch/polarization/POLARIZATION_SURVEY.md
6
```

```
$ grep -n -o 'S11bB_interface_assembly.md:53[^;]*' _scratch/polarization/POLARIZATION_SURVEY.md | head -4
179:S11bB_interface_assembly.md:53–54** says nonreciprocal “odd” couplings occur in driven laboratory media, including active/chiral fluids, and “They are not unphysical”. Lines 55–63 give the standing rule: “A non-passive coupling is admissible only with a NAMED reservoir and a STATED power budget.” [Q4](#q4) | **Correction and conditional admissibility rule**, not a blanket thermodynamic exclusion or a computed birefringence coefficient for this brane. These are the record's statements about admissibility, not adoption of an odd/chiral optical law. |
223:S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget
224:S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget
226:S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget
```

```
$ sed -n '17p;53,55p;76p' research/pde_ledger_v3/steps/S11bB_interface_assembly.md
## ⭐⭐⭐ THE HEADLINE — the velocity-coupled leak costs an ENERGY RESERVOIR
⚠ Non-reciprocal, "odd" constitutive couplings of exactly the excluded form are realised in driven
laboratory media (odd elasticity, odd viscosity, active and chiral fluids). ⛔ They are not unphysical.
⭐ **They have a reservoir.**
structural reason: in-plane parity admits no `e_W ↔ u_T` bilinear.
```

```
$ sed -n '284p;288p' research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md
- **In-plane isotropy _and_ parity** — the full `O(3)` acting on the three in-plane directions.
- ⛔ **Time-reversal is NOT assumed**, and ⛔ no positivity or boundedness is assumed (§0).
```

```
$ sed -n 358p research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md
- **In-plane isotropy _and_ parity** (the full `O(3)` acting on the three in-plane directions).
```

```
$ grep -c 'S11bB\?_SHARED_PHYSICS' _scratch/polarization/POLARIZATION_SURVEY.md
3
```

```
$ sed -n 110,111p research/pde_ledger_v3/steps/O2_steady_brane_balance.md
loads and boundary data remain live with spatial derivatives. Radial profiles impose no constitutive
isotropy, parity, stress symmetry, absent couple/chiral content, derivative cutoff or finite history
```

```
$ cat _scratch/polarization/survey_review_r3_claude_evidence/chiral_passive_check_stdout.txt
K(k) - K(k)^dagger (Hermiticity residual): Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
K(k) - K(-k)^T (Onsager reciprocity residual): Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
K(-k) - K(k) (parity-odd part): Matrix([[0, -2*I*eta*k_z*(k_x**2 + k_y**2 + k_z**2), 2*I*eta*k_y*(k_x**2 + k_y**2 + k_z**2)], [2*I*eta*k_z*(k_x**2 + k_y**2 + k_z**2), 0, -2*I*eta*k_x*(k_x**2 + k_y**2 + k_z**2)], [-2*I*eta*k_y*(k_x**2 + k_y**2 + k_z**2), 2*I*eta*k_x*(k_x**2 + k_y**2 + k_z**2), 0]])
K for k along z: Matrix([[k**2*mu, I*eta*k**3, 0], [-I*eta*k**3, k**2*mu, 0], [0, 0, 0]])
eigenvalue rho*omega^2 = 0  multiplicity 1  eigenvector [[0, 0, 1]]
eigenvalue rho*omega^2 = k**2*(-eta*k + mu)  multiplicity 1  eigenvector [[-I, 1, 0]]
eigenvalue rho*omega^2 = k**2*(eta*k + mu)  multiplicity 1  eigenvector [[I, 1, 0]]
eta=0 eigenvalues: {k**2*mu: 2, 0: 1}
```

## Cause of finding 2: the round-1 disposition and repair-2 brief listed S11bB:53–54 among the parity sources
```
$ grep -n -o 'S11bB:53–54[^;|]*' _scratch/polarization/survey_r1_review_disposition.md
21:S11bB:53–54 (odd couplings are physical)
```

```
$ grep -n 'S11bB_interface_assembly.md:53' _scratch/polarization/survey_repair2_prompt.md
30:   - `research/pde_ledger_v3/steps/S11bB_interface_assembly.md:53–54`;
```

