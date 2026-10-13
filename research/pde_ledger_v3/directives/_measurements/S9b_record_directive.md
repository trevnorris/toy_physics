# Measurements: S9b record directive (generated 2026-10-13 11:27)

Generator: `_scratch/s9b_build/gen/s9b_record_directive_lookups.sh` (git, sed, grep, sha256sum, wc only; this file is
written by it). Each block grounds a statement the directive makes about a source.

## The directive
```
$ sha256sum research/pde_ledger_v3/directives/S9b_record_directive.md _scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md
cdead60d6f39e82969afacb8d6b4140c6729b288279f14d5acb74c2a655b8994  research/pde_ledger_v3/directives/S9b_record_directive.md
9ad9611acf5e71795c0610e9f65f0554cf85570e0872e7e746e7522104285e94  _scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md
```

## Spec at 34a994d5: Parts A–D object lines (Part B's conditions and implied j_n; Part C; Part D)
```
$ git show 34a994d5:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | sed -n '228,256p'
- **Part A.** The branch-existence and path-traversal conditions. Then, for the LAB_HELD profile, as expressions
  or functionals of `δ`, `V` and `ξ_w`, with each retained monomial printed separately:
  - `Δθ`;
  - the round-trip excess time;
  - the two one-way excess times and their nonreciprocal part.
- **Part B.**
  - An effective `γ` from `Δθ`, and one from the round trip's logarithmic slope `𝒮_RT` (observable 2).
  - Their residuals against the references.
  - For each observable, the solution condition on the profiles, relative to `GM`, for the comparison
    residual of the selected grade sum (see "Order counting") at every far-zone `b`. This is an object
    to compute, with its domain and unresolved inputs; it is not a demanded residual value. Other
    retained monomials are printed but not set against the references.
  - For each condition, the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic.
  - The difference between the two `γ`s, printed as an object.
- **Part C.** The Part B conditions rewritten, to first order in `f`, under three supplied responses of `c_γ`
  to the local bulk density, with `V` and `ξ_w` left live in each:
  - `c_γ(x) ≡ c₀` (`δ(x) ≡ 0`);
  - pointwise fixed ratio: `c_γ(x)/c_s(x) = c₀/c_s0`;
  - symbolic power response: `c_γ(x)/c₀ = (ρ(x)/ρ₀)^s`, with `s` symbolic.

  Here `ρ(x)` is the bulk number density at the brane and `ρ₀` is its asymptotic value. For each condition,
  print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
  if it enters at all.

- **Part D.** The Part B conditions re-expressed with the linked neutral-sector steady in-plane
  momentum balance below in force. Print the constructed balance pieces and conditional objects,
  their premises and remaining OPEN inputs, which of `δ`, `V`, `ρ_br` and `j_n` each condition
  determines relative to the independent orbital `GM`, which stay free, and the implied `j_n` for
  each condition. For each condition, also print what the OPEN pieces left in the in-plane balance
```

## Part D held, not accepted
```
$ git show -s --format=%B 97580010 | sed -n '1,16p'
Accept S9b spec v10 repair 3 for Parts A–C; Part D held, not accepted

Repair 3 (fresh Codex author) removes the induced-measure readings: S9b
computes and claims on the coordinate d³x measure only (user, 2026-10-08).
Review round 4 (fresh Claude + Grok, identical prompt): Grok no findings;
Claude one finding, F1, accepted. The coordinate-only rule is unscoped
and, read literally, would drop the native face-area and measure content
that Part D keeps OPEN. The unscoped wording came from my repair brief.

Parts A–C are sound per both legs and are accepted as the governing
physics of the repair build. Part D is held (user, 2026-10-08) until the
brane-material premise behind P1 is settled, and is not accepted; F1's
fix, the stale version line and the narrowing's provenance go with it.
Under the agreed hard stop the user chose this disposition.

Dispositions: directives/_measurements/S9b_spec_v10_r4_review_disposition.md
```

## Decision list at 79055918: D2 premises and D5 obligation
```
$ git show 79055918:research/pde_ledger_v3/directives/S9b_repair_decision_list.md | sed -n '21,34p'
## D2. Adopted premises
Each premise below is user-selected and is an **adopted substrate input to a conditional model**. None is derived.
Every result that depends on one is flagged with it.

**From O2 (user-selected 2026-10-06):** P1–P4, `directives/O2_premise_decision_list.md` (`77d2c39a`).

**New for S9b:**

| Premise | Content | Qualification carried with it |
|---|---|---|
| **P5** (2026-10-07; carriage 2026-10-08) | The brane's in-plane momentum density is `ρ_br V`, and that momentum is carried with the brane material at `V`. Its contribution to the in-plane momentum current is `ρ_br V^i V^j`. | Scoped to Part D. Any other in-plane momentum current stays OPEN. It is counted beside the supplied current and P6's stress, which it does not contain again. Part D states what the OPEN pieces must supply. If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |
| **P6** (2026-10-08, revised the same day) | In steady flow the brane's in-plane stress is an isotropic pressure plus a linear viscous stress: `T^{ij} = −p_br(ρ_br) δ^{ij} + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k`, in the Cauchy convention (traction `T^{ij} n_j`; force density `∂_j T^{ij}`). `p_br` is a general function of `ρ_br` only. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. The shear viscosity `η` and bulk viscosity `ζ` are live general profiles. In the optical regime, light still sees the elastic transverse stiffness `μ_⊥`. | Scoped to Part D. It supplies the steady in-plane part of O2's stress `𝒯_br^live` only. It is an adopted steady form, not a consequence of P1: P1's relaxation response stays OPEN outside it. The linear viscous form excludes power-law creep. The power this stress dissipates is not identified with any O2 energy operand. It stays in the OPEN energy accounting, with its physical supplier and budget OPEN. The first 2026-10-08 row (pressure only, justified by relaxation under steady load) is withdrawn: steady flow keeps straining, so relaxation does not remove the viscous stress (spec v10 review C2). |
| **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥ ∝ ρ_br^α` and `α` one live symbol. This is the user's selected live exponent. | Scoped to Part D. There it supplies `μ_⊥` as a function of `ρ_br` alone. Its other O1 `ℳ_⊥` dependences are set aside by this premise, and results are flagged. Parts A–C keep `μ_⊥(x)` and `ρ_br(x)` independent, as in v8. |
| **w-parity** (2026-10-07) | For an electrically neutral mass, the far field is symmetric under `w → −w`. So `ξ_w` and the material `w`-velocity vanish there. | This is a **labelled neutral-sector restriction**. Parts A–C keep `ξ_w` live as in v8, which covers the charged case. Part D is stated for the neutral sector, and prints the label. |
```

```
$ git show 79055918:research/pde_ledger_v3/directives/S9b_repair_decision_list.md | sed -n '70,87p'
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
```

## Build directive at 4b4db191 (source of the build objects the directive defines; not handed to the record author)
Tag families:
```
$ git show 4b4db191:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sed -n '120,146p'

   | Name | Object |
   |---|---|
   | `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE` | the two conditions and the branch type of the spec's "Branch existence" |
   | `A_DEFLECTION_G<abc>` | Δθ, at that grade |
   | `A_ROUND_TRIP_G<abc>` | the round-trip excess time, at that grade |
   | `A_ONE_WAY_ER_G<abc>`, `A_ONE_WAY_RE_G<abc>` | the one-way excess time from emitter to reflector, and from reflector to emitter, at that grade |
   | `A_NONRECIPROCAL_G<abc>` | their nonreciprocal part, at that grade: half of (emitter-to-reflector time minus reflector-to-emitter time) |
   | `A_NONRECIPROCAL_PATH_DEPENDENCE` | item 8's object |
   | `A_RADAR_LOG_G<abc>` | item 7's slope, for each grade in Part B's comparison sum (`G100`, `G010`, `G020`, `G001`) |
   | `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE` | Part B's effective `γ`s, and their difference: the deflection `γ` minus the radar `γ` |
   | `B_RESIDUAL_<OBS>` | Part B's residuals against the references |
   | `B_CONDITION_<OBS>`, `B_IMPLIED_JN_<OBS>` | Part B's every-`b` condition, and the `j_n` it implies |
   | `C_<RESP>_CONDITION_<OBS>`, `C_<RESP>_IMPLIED_JN_<OBS>` | Part C's rewritten condition, and the `j_n` it implies |
   | `C_<RESP>_N_DEPENDENCE` | where `n` enters that row's conditions |

   Item 10's objects are named by a prefix on the name of the object they restrict:
   - `R_FLOW_ONLY_`, `R_SPEED_ONLY_` and `R_TILT_ONLY_` prefix `A_DEFLECTION_G<abc>`, `A_RADAR_LOG_G<abc>`,
     `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_<OBS>`, `B_CONDITION_<OBS>` and `B_IMPLIED_JN_<OBS>`.
   - `R_BULK_ONLY_` prefixes `C_FIXED_RATIO_CONDITION_<OBS>`, `C_FIXED_RATIO_IMPLIED_JN_<OBS>`,
     `C_POWER_CONDITION_<OBS>` and `C_POWER_IMPLIED_JN_<OBS>`.
   - `F_FLOW_` (the forward case's first stage, `δ ≡ 0`, `ξ_w ≡ 0`) and `F_LIVE_` (its second stage) prefix
     `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE`, every graded `A_` name (`A_…_G<abc>`),
     `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_<OBS>` and `B_CONDITION_<OBS>`. `F_LIVE_` also prefixes
     `C_<RESP>_CONDITION_<OBS>`.
   - `F_MASS_SOLUTION` is the forward case's solved `V`, expressed through `Φ` and `ρ_br`.
5. **No verdicts (D6 B5).** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the
```

Item 9's density-sign premise and labels:
```
$ git show 4b4db191:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sed -n '177,187p'
9. **Part C domains (D6 B6).** Each Part C row's domain predicates are built from that row's own response. They are
   not built from another row, or from the amplitude before the response is applied. Each predicate follows from the
   row's response being real and positive, or from the row's supplied inputs. A supplied input is one the spec
   writes as an equation or an inequality, or the premise this item supplies.
   - **Premise (the user's choice, 2026-10-11):** the asymptotic bulk number density is positive, `ρ₀ > 0`. It is
     supplied here, so it is not a stop event under item 12. Each engine prints it as a supplied premise.
   - **The local bulk density's sign is not a domain predicate (the user's choice, 2026-10-11).** No row's domain
     excludes a point because `ρ = ρ₀(1 + f)` is zero or negative there.
   - The fixed-ratio and power rows read `ρ`. Each of them, with its implied-`j_n` rows and its item 10 copies,
     prints its objects with the part where `ρ > 0` and the part where `ρ ≤ 0` labelled separately. A knife's
     corrupted copy of such a row prints the same two parts. The constant row's response does not read `ρ`.
```

Item 10, the forward case and the single-mechanism restrictions:
```
$ git show 4b4db191:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sed -n '199,240p'
10. **Forward case and single-mechanism restrictions (user requests, 2026-10-08).**
    - **Forward case: no far-zone loss.** This is the user's premise (2026-10-06/07/08): brane material leaves the
      brane only in throats, so `j_n ≡ 0` in the far zone, and the same brane mass crosses every sphere around the
      body.
      - `Φ` is a new live symbol, defined by one equation: `Φ ≡ ∮ ρ_br V·n_out dA`, the brane mass flux through a
        far-zone coordinate sphere with outward unit normal `n_out`. It has no value and no sign, and its relation
        to `GM` is open.
      - Solve the supplied mass balance for `V` under this premise, with `ρ_br` live, and express the integration
        constant through `Φ`.
      - Substitute that `V` into Part A's observables, first with `δ ≡ 0`, `ξ_w ≡ 0`, then with `δ` and `ξ_w`
        live.
      - First print the branch-existence and path-traversal conditions after the same substitution. Gate the
        forward observables on them as Part A does.
      - Print the deflection, the round trip's logarithmic slope (item 7), both effective `γ`s, their
        difference and both Part B residuals against the references, each as a function of `b`. Also print the round-trip excess time, both
        one-way excess times and their nonreciprocal part, with `b`, `Z_E` and `Z_R` kept.
      - At both stages, also print Part B's condition for each observable under the same substitution, reduced as
        in item 6. Its free inputs are `Φ` and the profiles live at that stage: `ρ_br`, plus `δ` and `ξ_w` at the
        second stage. `j_n ≡ 0` is the case's premise, not a free input, so item 6's implied-`j_n` print does not
        apply. The condition is a computed object. It is not a supplied relation between `Φ` and `GM`.
      - Also print the second-stage condition rewritten under each of Part C's three responses, as Part C rewrites
        Part B's conditions. Reduce each as in item 6 and build its domain as in item 9. `Φ`, `ρ_br`, `f`, `ξ_w` and
        the row's `n` or `s` stay symbolic. These are computed objects too, and the implied-`j_n` print does not
        apply to them.
      - Label every object with the premise.
    - For the deflection and for the radar slope (item 7), print each restriction below of:
      - Part A's computed observable;
      - Part B's effective `γ` and residual;
      - the difference between the two effective `γ`s;
      - Part B's condition, reduced as in item 6, with its implied `j_n`.
    - The restrictions:
      - **flow only:** `δ ≡ 0`, `ξ_w ≡ 0`, `V` live;
      - **speed only:** `V ≡ 0`, `ξ_w ≡ 0`, `δ` live;
      - **tilt only:** `δ ≡ 0`, `V ≡ 0`, `ξ_w` live.
    - For Part C's fixed-ratio and power responses, also print each condition restricted to `V ≡ 0`, `ξ_w ≡ 0`
      (bulk density alone), with its domain and the `j_n` it implies through the supplied mass balance, with
      `ρ_br` symbolic.
    - Each restriction is obtained by substituting into the computed general object, then reducing. Its tag and
      label name the restriction and the profiles it sets to zero, never an outcome. Restrictions add no premise.
    - Everything in this item is additional labelled output. This directive supplies the forward case's premise for
      that label only. Building this item is not a stop event under item 12, and it does not replace the live `j_n`
      in Parts A–C.
```

## Engines unchanged since 503736ff (commits touching them after it, then the last commit touching each)
```
$ git log --format=%h 503736ff..HEAD -- research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl | wc -l
0
```

```
$ git log -1 --format='%h %s' -- research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py; git log -1 --format='%h %s' -- research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
bfa034c4 Preserve the S9b repair build after repair round 5 (review r5 baseline; not accepted)
bfa034c4 Preserve the S9b repair build after repair round 5 (review r5 baseline; not accepted)
```

## Production run record exists
```
$ ls -l research/pde_ledger_v3/directives/_measurements/S9b_production_runs.md
-rw-rw-r-- 1 trevnorris trevnorris 4724 Oct 12 08:00 research/pde_ledger_v3/directives/_measurements/S9b_production_runs.md
```

## O2 record at 72866fcf: section headings, and the ξ_w'' obligation lines
```
$ git show 72866fcf:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -n '^## '
26:## 1. Object and evidence convention
72:## 2. Exact domain, adopted premises and supplied inputs
137:## 3. Cross-engine evidence and every printed row class
220:## 4. Complete difference ledger and classification
336:## 5. Printed comparator scope and acceptance limits
368:## 6. Items routed from engine acceptance to the comparator
388:## 7. Remaining OPEN inputs and owners
421:## 8. Register pass and S9b Part D handoff
517:## 9. Source tensions and questions requiring more than retrieval
```

```
$ git show 72866fcf:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | grep -c "xi_w''"
10
```

## Register: entry schema and population calibration
```
$ sed -n '40,60p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
## Entry schema

| field | meaning |
|---|---|
| `id` | stable identifier, cited from both ends |
| `source` | the step that needs it |
| `target` | the step that must deliver it |
| `requirement` | **the object required**, named — not a derivation path for it |
| `on failure` | what breaks in the source step if the target cannot deliver |
| `status` | OPEN · DELIVERED · RETIRED (with the commit that changed it) |

Name the object; do not prescribe how the target should obtain it. A requirement that specifies a
derivation route manufactures arguments about the route — see `CLAUDE.md` rule 3.

Distinguish carefully:
- a **requirement** is an obligation on a *future* step;
- a **postulate** is a value or form assumed *here* and not derived — that belongs in `ANSATZ_LEDGER.md`;
- a **defect** is something already wrong — that belongs in `DEFECT_REGISTER.md`.

A postulate with a named retirement condition generates a requirement. `B_comp` below is exactly that.

```

```
$ sed -n '521,540p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
## Population passes

**Pass 1 — S9 and S10. ⭐ DONE 2026-08-06**, alongside S10's step record and from the same reading.
Seven entries added, against the two that existed. The calibration it produced, for pass 2:

- ⭐ **Consolidate rather than duplicate.** S9 states its shear requirement twice — once for the ordered
  phase and once for the disordered — and a third consumer (the geon's trapped mode) appears elsewhere in
  the record. It is **one object**, `R-S1-02`, with three consumers, ⛔ not three entries.
- ⚠ **Distinguish a requirement from an ansatz by asking whether a retirement condition is LIVE.** `μ_R`'s
  **magnitude** is postulated and defect `B2` closed the polar-`P` route to deriving it ⇒ ansatz, ⛔ no
  entry. Its **form** has no such closure ⇒ `R-S8-01`, a live requirement. The two sit in the same
  sentence of S9's record and go to different registers.
- ⛔ **A requirement is not everything a later step would find useful.** A simulation needs `ρ_br` and
  `μ_R` separately rather than only their ratio; that is a *consumer's* need pointing **forward**, ⛔ not
  an obligation on which a banked result rests. It got no entry.
- ⭐ **The load-bearing ones are the identifications nobody computed.** `R-S8-02` — that the in-plane and
  out-of-plane sectors decouple — is the single entry on which S10's headline number depends, and it was
  flagged inside S10 **at the moment the identification was made**. ⇒ read a record for its *flagged
  identifications* first; they are where the requirements are.

```

## CLAUDE.md step-record row
```
$ grep -n '^| \*\*Step record' CLAUDE.md
58:| **Step record / `.tex` card / physics-bearing prose** | Two; O or Cx; ⛔ never chosen by extension | Source-first fidelity; quote both sides; ⛔ no build directive in the packet; a measured physics claim still needs its command + stdout; for cards check suppressed macro fields | **Review-until-clear** about what may be claimed | Own artifact/version review required; script review alone is not record review. Both reports before commit. |
```

