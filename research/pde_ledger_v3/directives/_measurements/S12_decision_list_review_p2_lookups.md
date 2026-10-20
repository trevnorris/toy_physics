# Measurements: S12 decision list, review pass 2 (generated 2026-10-19 22:53)

Generator: `_scratch/s12/gen/dl_p2_lookups.sh`. Literal retrieval only. The list is read from the working tree;
M1 shows its hash equals the reviewed baseline. The reviewed revision is preserved in the commit that adds this file.

## M1. The pass: prompt, baseline, exits, verdict lines
```
$ sha256sum /var/projects/toy_physics/_scratch/s12/dl_review_prompt_p2.md | cut -c1-64
b35adc417e170c92869e847e99e55bf2c96053fd0c14499578384d8e6f52ebc2
```

```
$ sha256sum -c /var/projects/toy_physics/_scratch/s12/dl_baseline_p2.sha256 | sed 's#/var/projects/toy_physics/##'
research/pde_ledger_v3/directives/S12_decision_list.md: OK
```

```
$ cat /var/projects/_review_scratch/s12-dl-p2-codex/DONE /var/projects/_review_scratch/s12-dl-p2-grok/DONE
EXIT=0
EXIT=0
```

```
$ grep -n 'FINDINGS\|^[0-9]\. \*\*' /var/projects/_review_scratch/s12-dl-p2-codex/REVIEW.md /var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt | cut -c1-170
/var/projects/_review_scratch/s12-dl-p2-codex/REVIEW.md:93:**FINDINGS**
/var/projects/_review_scratch/s12-dl-p2-codex/REVIEW.md:95:1. **MEDIUM — The K1 repair enforces U1's exterior support condition but omits its settling/reordering requir
/var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt:14:**FINDINGS**
/var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt:16:1. **Medium.** The fast-ordering specialization is not one source. The list leaves the profile as a functi
/var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt:46:2. **Medium.** The S14a open item widens the controlled return law into "the return law." The plan's fourt
/var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt:71:3. **Medium.** The boundary inventory cites the IR/return scheme and then states a smaller scheme. The cit
/var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt:94:4. **Medium.** The bulk-direction carry is filed as open with S12, together with S12's additional momentum
/var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt:122:5. **Low.** The S9b/O2 domain bullet cites lines 74–75 and omits restrictions those lines state. The re
```

## M2. U1's re-ordering clause (Codex 1)
```
$ grep -n 'Outside it,\|Inside it,\|on order fields U1 allows' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
58:- **Outside it,** `χ_B` vanishes.
59:- **Inside it,** `χ_B` is a live field, between 0 and 1. It may depend on history or strain. No profile, rate or sign
69:1. **The source the committed law requires on order fields U1 allows.** The object is `n Γ_B`, with `n`, `u`, `J_χ`,
121:The source the committed order law requires on order fields U1 allows, with every field live, and its labelled
```

```
$ grep -n 'settles into the slab re-orders' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
44:  de-orders it. Bulk material that settles into the slab re-orders. The user confirmed this one rule.
```

```
$ sed -n '5p' /var/projects/toy_physics/docs/oriented_throat_return_source_correction.md | grep -o 'The committed meaning is.\{0,110\}'
The committed meaning is **bulk material re-ordering into the brane state of the same medium**. A local loop back to the original thro
```

## M3. The fast-ordering specialization's arguments (Grok 1)
```
$ grep -n 'equivalently a live centre and thickness\|a given function of' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
56:  live faces `w₋(x,t) < w₊(x,t)`, equivalently a live centre and thickness. The faces are properties of the slab, set
73:2. **A labelled conditional specialization.** The same object where the interior order is a given function of
```

## M4. The S14a return law (Grok 2)
```
$ grep -n 'the return law (S14a' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
115:- the projected order-loss source, `J`, the `J → Q` map and the return law (S14a; `V3_STEP_PLAN.md:636–637`);
```

```
$ sed -n '636,637p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
⇒ **Derive, from the dynamical system:** the projected order-loss source · the definition of `J` · the
profile-dependent `J → Q` map · the controlled `ω → 0` return law.
```

```
$ sed -n '606p' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | cut -c1-260
| **S14a / S14** | Projected order-loss source, drain label `J`, profile-dependent far-field flux map and controlled DC return bridge remain S14a obligations. S14 is **CONDITIONAL ON S14a**, with its recorded frozen DC-sink/zero-mode completion domain. Those c
```

## M5. The IR/return scheme (Grok 3)
```
$ grep -n 'one- or two-sided return' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
86:  (`V3_STEP_PLAN.md:612`). The IR/return scheme stays OPEN: one- or two-sided return, the return or absorber law, and
```

```
$ sed -n '340p' /var/projects/toy_physics/docs/two_throat_simulation_handoff_spec.md
| IR / return scheme | One- vs two-sided domain, return/absorber law, far-field truncation/extrapolation, treatment of zero modes/moduli, and (for moving/radiative rungs) the outgoing/absorbing prescription. |
```

```
$ sed -n '50p' /var/projects/toy_physics/research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md
- **The phase / moduli zero-mode quotients** (the return controller lives in the zero-mode complement).
```

## M6. The bulk-direction carry (Grok 4)
```
$ grep -n 'the bulk-direction carry, which is not' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
95:  - the bulk-direction carry, which is not `j_n` times one velocity (`directives/O2_SHARED_PHYSICS.md:243`);
```

```
$ grep -n 'OPEN .𝒥_map. and the normal material response\|with .𝒥_map. and .𝒩_br^live. named' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | cut -c1-260
246:The reduced component therefore depends on the OPEN `𝒥_map` and the normal material response
584:| Carried exchange-momentum components and any separately attributable source-partner actions | `Π_n`, premise 3's in-plane `j_n V^i`, the bulk-direction component as the OPEN reduction of §5 with `𝒥_map` and `𝒩_br^live` named, O6/native relative f
```

```
$ grep -n '^| .𝒩_br^live.' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | cut -c1-120
134:| `𝒩_br^live`; **OPEN**, C §3 | Normal material response, keeping embedding, centre and thickness content distin
```

```
$ sed -n '396p' /var/projects/toy_physics/research/pde_ledger_v3/steps/O2_steady_brane_balance.md | cut -c1-260
| O3 `Π_n` | Exchange operand: only local in-plane carried momentum is closed by adopted premise 3; bulk carry/reactions remain OPEN. | S12 additional partners/reaction system; native reduction through O6/normal response. |
```

## M7. The S9b/O2 domain sentence (Grok 5)
```
$ grep -n 'The S9b/O2 domain, as recorded' -A1 /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
89:- **The S9b/O2 domain, as recorded.** One isolated mass at rest with its drain flowing, Eulerian steady profiles,
90-  and lab time in the brane's far-field rest frame (`steps/O2_steady_brane_balance.md:74–75`). Its restrictions
```

```
$ sed -n '74,76p' /var/projects/toy_physics/research/pde_ledger_v3/steps/O2_steady_brane_balance.md
The supplied setting is one isolated, spherically symmetric mass at rest with its drain flowing:
far field, Eulerian steady material profiles, lab time in the brane's far-field rest frame, linear
optical waves and leading eikonal. There are three in-plane coordinate directions `x^i` and bulk
```

