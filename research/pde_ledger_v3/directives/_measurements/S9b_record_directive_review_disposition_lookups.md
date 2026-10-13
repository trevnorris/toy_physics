# Lookups: S9b record directive G2 review dispositions (generated 2026-10-13 11:27)

Generator: `_scratch/s9b_build/gen/s9b_record_dl_review_lookups.sh` (sed, grep, git, sha256sum only; this file is
written by it). The reviewed text is the frozen copy `_scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md`.

## Reviewed baseline and the legs' reports
```
$ sha256sum _scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md; cat _scratch/s9b_build/s9b_record_directive_review_baseline.sha256
9ad9611acf5e71795c0610e9f65f0554cf85570e0872e7e746e7522104285e94  _scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md
9ad9611acf5e71795c0610e9f65f0554cf85570e0872e7e746e7522104285e94  research/pde_ledger_v3/directives/S9b_record_directive.md
ab37556b43ec18ea0cd27745314b67829b27632eb34f12763469820638089092  _scratch/s9b_build/s9b_record_directive_review_prompt.md
```

```
$ wc -c _scratch/s9b_build/s9b_record_directive_review_codex_final.txt _scratch/s9b_build/s9b_record_directive_review_grok.txt; cat _scratch/s9b_build/s9b_record_directive_review_codex.DONE _scratch/s9b_build/s9b_record_directive_review_grok.DONE
 1647 _scratch/s9b_build/s9b_record_directive_review_codex_final.txt
 9136 _scratch/s9b_build/s9b_record_directive_review_grok.txt
10783 total
EXIT=0
EXIT=0
```

## C1 (Codex): the build directive in the record packet
```
$ sed -n '26,27p' _scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md
- **The build contract:** `directives/S9b_repair_build_directive.md` at `4b4db191`, with amendments 1–10. Item 10
  defines the forward case and the single-mechanism restrictions.
```

```
$ grep -n '^| \*\*Step record' CLAUDE.md
58:| **Step record / `.tex` card / physics-bearing prose** | Two; O or Cx; ⛔ never chosen by extension | Source-first fidelity; quote both sides; ⛔ no build directive in the packet; a measured physics claim still needs its command + stdout; for cards check suppressed macro fields | **Review-until-clear** about what may be claimed | Own artifact/version review required; script review alone is not record review. Both reports before commit. |
```

## G1 (Grok): item 8's entry object, target and status
```
$ sed -n '96,106p' _scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md
8. **The Part D register entry.** The S9b pass adds one register entry for the need Part D was to meet: the brane's
   steady in-plane force balance must deliver the brane flow, density and normal loss that the Part B conditions
   require. The entry:
   - cites the printed conditions and restrictions it rests on;
   - carries P5, P6 and the density link (`c_γ² = μ_⊥/ρ_br`, with `μ_⊥ ∝ ρ_br^α` and `α` live) as conditional inputs
     that a linked brane would use, as D2 states them. Under the density link, the light speed is a function of the
     brane density through `α`; in Parts A–C the speed change and the brane density are independent live profiles.
     The entry states this, and no effect of the link;
   - names the candidate material behaviours raised so far (power-law creep; a two-component supersolid, a superfluid
     through a lattice) as knit candidates, none selected and the list not closed;
   - targets the knit, follows the register's schema and rest-on test, and states its status.
```

```
$ sed -n '44,56p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
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
```

```
$ sed -n '533,535p' research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
- ⛔ **A requirement is not everything a later step would find useful.** A simulation needs `ρ_br` and
  `μ_R` separately rather than only their ratio; that is a *consumer's* need pointing **forward**, ⛔ not
  an obligation on which a banked result rests. It got no entry.
```

```
$ git show 34a994d5:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | sed -n '233,241p'
- **Part B.**
  - An effective `γ` from `Δθ`, and one from the round trip's logarithmic slope `𝒮_RT` (observable 2).
  - Their residuals against the references.
  - For each observable, the solution condition on the profiles, relative to `GM`, for the comparison
    residual of the selected grade sum (see "Order counting") at every far-zone `b`. This is an object
    to compute, with its domain and unresolved inputs; it is not a demanded residual value. Other
    retained monomials are printed but not set against the references.
  - For each condition, the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic.
  - The difference between the two `γ`s, printed as an object.
```

```
$ git show 79055918:research/pde_ledger_v3/directives/S9b_repair_decision_list.md | sed -n '86,87p'
- **Register.** A register entry follows only if a sourced requirement is established. An unresolved difference
  alone does not create one (record's rest-on criterion).
```

```
$ git show 72866fcf:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | sed -n '466p' | grep -o 'No complete action, momentum, normal, work or energy balance is established equal or solved.'
No complete action, momentum, normal, work or energy balance is established equal or solved.
```

## G2 (Grok): Part D premises and named material candidates on the entry
```
$ git show 79055918:research/pde_ledger_v3/directives/S9b_repair_decision_list.md | sed -n '31,33p' | cut -c1-400
| **P5** (2026-10-07; carriage 2026-10-08) | The brane's in-plane momentum density is `ρ_br V`, and that momentum is carried with the brane material at `V`. Its contribution to the in-plane momentum current is `ρ_br V^i V^j`. | Scoped to Part D. Any other in-plane momentum current stays OPEN. It is counted beside the supplied current and P6's stress, which it does not contain again. Part D state
| **P6** (2026-10-08, revised the same day) | In steady flow the brane's in-plane stress is an isotropic pressure plus a linear viscous stress: `T^{ij} = −p_br(ρ_br) δ^{ij} + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k`, in the Cauchy convention (traction `T^{ij} n_j`; force density `∂_j T^{ij}`). `p_br` is a general function of `ρ_br` only. Its compressional s
| **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥ ∝ ρ_br^α` and `α` one live symbol. This is the user's selected live exponent. | Scoped to Part D. There it supplies `μ_⊥` as a function of `ρ_br` alone. Its other O1 `ℳ_⊥` dependences are set aside by this premise, and results are flagged. Parts A–C keep `μ_⊥(x)` and `ρ_br(x)` independent, as in
```

```
$ git show 34a994d5:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | grep -n 'power-law creep'
374:form excludes power-law creep. P6 supplies no normal stress or normal material law, relaxation
539:  form, with general live `η` and `ζ`. P6 is not derived from P1; power-law creep is outside this model
```

```
$ git grep -c -i -e 'two-component' -e 'supersolid' -e 'through a lattice' 4b4db191 -- research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md research/pde_ledger_v3/directives/S9b_repair_decision_list.md research/pde_ledger_v3/directives/S9b_repair_build_directive.md research/pde_ledger_v3/steps/O2_steady_brane_balance.md research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md; echo "git grep exit $?"
git grep exit 1
```

## The folded text
```
$ sha256sum research/pde_ledger_v3/directives/S9b_record_directive.md
cdead60d6f39e82969afacb8d6b4140c6729b288279f14d5acb74c2a655b8994  research/pde_ledger_v3/directives/S9b_record_directive.md
```

```
$ grep -n -i -e 'build directive' -e 'build contract' -e 'creep' -e 'supersolid' -e 'item 10' research/pde_ledger_v3/directives/S9b_record_directive.md
9:claimed. No build directive is in the record author's or the record legs' packet (`CLAUDE.md`, step-record row). The
```

```
$ sed -n '124,133p' research/pde_ledger_v3/directives/S9b_record_directive.md
8. **What Part D was to meet.** If the register's rest-on test admits it, the S9b pass adds an entry whose
   requirement is the printed Part B condition for each observable, with its Part C rewrites: the combination of
   live profiles it constrains relative to `GM`, and the `j_n` it implies, with `ρ_br` symbolic and every quantity
   the print leaves free kept free. The entry:
   - cites the printed conditions and restrictions it rests on;
   - names no means of meeting the requirement: no balance, material, premise or candidate;
   - states no effect of the density link;
   - has a step id or "no owner named" as its target, and a status from the register's schema.

   If the rest-on test does not admit it, the record says why.
```

