# S12 decision list review pass 4: mechanical lookups

Generator: `_scratch/s12/gen/dl_p4_lookups.sh`. Verbatim retrieval and hashes only.

## P1. Baseline, prompt, working tree, exits
```
$ sha256sum -c _scratch/s12/dl_baseline_p4.sha256
research/pde_ledger_v3/directives/S12_decision_list.md: OK
$ sha256sum _scratch/s12/dl_review_prompt_p4.md
c207d863cf2c8238a94da6690355f59b83466634431326960609dab817bccd9d  _scratch/s12/dl_review_prompt_p4.md
$ git status --porcelain | diff _scratch/s12/status_before_dl_p4.txt -  (empty = unchanged)
2a3,4
> ?? research/pde_ledger_v3/directives/_measurements/S12_decision_list_review_p4_disposition.md
> ?? research/pde_ledger_v3/directives/_measurements/S12_decision_list_review_p4_lookups.md
diff exit: 1
$ cat /var/projects/_review_scratch/s12-dl-p4-codex/DONE /var/projects/_review_scratch/s12-dl-p4-grok/DONE
EXIT=0
EXIT=0
```

## P2. Verdicts
```
$ tail -1 /var/projects/_review_scratch/s12-dl-p4-codex/REVIEW.md
**CLEAR**
$ grep -n -o 'Revision 3 does not clear[^.]*\.' /var/projects/_review_scratch/s12-dl-p4-grok/grok_stdout.txt
1:Revision 3 does not clear.
$ grep -n -o '^\*\*1\. High\.\*\*[^.]*\.' /var/projects/_review_scratch/s12-dl-p4-grok/grok_stdout.txt
17:**1. High.** Item 1 identifies `n Γ_B` with every face-concentrated part.
```

## P3. The two legs on S4-3
```
$ grep -o '| \*\*S4-3\*\* |[^|]*|' /var/projects/_review_scratch/s12-dl-p4-codex/REVIEW.md
| **S4-3** | **Holds.** The enumerated input fields have no concentrated component, and the OPEN item says “it is not in the object at this step, and it stays an open possibility.” The source still includes derivative-generated face terms. The class declaration's explicit subject is the four named input fields; it does not remove the face-supported part of the computed `n Γ_B` (E1/E2). |
$ grep -o 'The reviewed object is the combined source[^.]*\.' /var/projects/_review_scratch/s12-dl-p4-codex/REVIEW.md
The reviewed object is the combined source `n Γ_B`; this review does not divide it by discontinuous `n` or infer a separate singular `Γ_B` representative.
$ grep -n -A1 '^\*\*S4-3 does not hold' /var/projects/_review_scratch/s12-dl-p4-grok/grok_stdout.txt
7:**S4-3 does not hold.** The OPEN item says a part of a field concentrated on a face is not in the object. Item 1 says the object covers every part concentrated on a face.
8-
$ grep -o 'A build that follows item 1[^.]*\.' /var/projects/_review_scratch/s12-dl-p4-grok/grok_stdout.txt
A build that follows item 1 either writes the order-flux jump into `n Γ_B`, which is a face-concentrated field part the OPEN item excludes, or enforces only the material requirement and drops the order one.
```

## P4. The face coefficient each leg computed for the order balance
```
$ grep -o 'the computed delta coefficients are[^.]*\.' /var/projects/_review_scratch/s12-dl-p4-codex/REVIEW.md
the computed delta coefficients are `m_R-m_L` for total material and `c_R*m_R-c_L*m_L+j_R-j_L` for the order source.
$ grep -o 'U1 gives the lower-face source[^.]*\.' /var/projects/_review_scratch/s12-dl-p4-codex/REVIEW.md
U1 gives the lower-face source `C*m_R+j_R-j_L` and upper-face source `-C*m_L+j_R-j_L`.
$ grep -n 'C_O\|C_M\|n_Gamma' /var/projects/_review_scratch/s12-dl-p4-grok/evidence/09_independent_face_requirements.stdout.txt
9:lower_C_M 0
10:lower_C_O 4
11:lower_n_Gamma_below 0
12:lower_n_Gamma_above 0
13:upper_C_M 0
14:upper_C_O -4
15:upper_n_Gamma_above 0
16:upper_n_Gamma_below 0
17:density_jump_C_M 29
18:density_jump_C_O 0
19:matched_C_M 0
20:matched_C_O 0
```

## P5. The list: revision 3 text, and revision 2's item 1 and U1 at HEAD
```
$ grep -n 'every part concentrated on a face\|Fields at a face\|no field has a part\|any part of a field concentrated' research/pde_ledger_v3/directives/S12_decision_list.md
19:order carried along a face out for now. Revision 3 states that in item 1, under "Fields at a face", and adds one
83:   distinguishes: outside the region, inside it, and at its faces, including every part concentrated on a face. No
86:   **Fields at a face.** This holds for `n`, `u` and `J_χ`, and for the order inside the region (U1 sets it to zero
91:   - no field has a part concentrated on a face itself.
128:- any part of a field concentrated on a face, such as order carried along a face. The user chose "Leave it out for
$ git show HEAD:research/pde_ledger_v3/directives/S12_decision_list.md | grep -n -A1 'distinguishes: outside the region'
77:   distinguishes: outside the region, inside it, and at its faces, including any part concentrated on a face if the
78-   order jumps there. No flow profile is imposed, and the river profile does not enter.
$ grep -n 'past that thickness\|de-orders it\|may be continuous or may jump' research/pde_ledger_v3/directives/S12_decision_list.md
52:  inside the slab's thickness. A throat, or a pinch hard enough, that carries material past that thickness
53:  de-orders it. Bulk material that settles into the slab re-orders. The user confirmed this one rule.
70:- **At a face,** `χ_B` may be continuous or may jump. Neither is chosen.
```
