# S12 build review r1: mechanical lookups

Generator: `_scratch/s12/gen/build_r1_lookups.sh`. Verbatim retrieval, hashes and existence checks only.

## C1. Baseline, prompt and working tree
```
$ sha256sum -c _scratch/s12/build_review/s12_build_baseline_r1.sha256
research/pde_ledger_v3/scripts/S12_order_source.py: OK
research/pde_ledger_v3/mathematica/S12_order_source.wl: OK
research/pde_ledger_v3/scripts/out/S12_order_source.out: OK
research/pde_ledger_v3/mathematica/out/S12_order_source.out: OK
research/pde_ledger_v3/directives/_measurements/S12_build_run_record.md: OK
$ sha256sum _scratch/s12/build_review/s12_build_review_prompt_r1.md
531de6604b755f5b7d575c9f76bc74c188a0a9170da3dad8c6943ea5a57783c1  _scratch/s12/build_review/s12_build_review_prompt_r1.md
$ git status --porcelain | diff _scratch/s12/build_review/status_before_r1.txt -  (empty = unchanged)
6a7,8
> ?? research/pde_ledger_v3/directives/_measurements/S12_build_review_r1_disposition.md
> ?? research/pde_ledger_v3/directives/_measurements/S12_build_review_r1_lookups.md
diff exit: 1
```

## C2. The two legs' verdict lines
```
$ grep -n '^## Verdict\|^\*\*F[0-9]' /var/projects/_review_scratch/s12-build-r1-claude/REVIEW.md
3:## Verdict: NEEDS REVISION
9:**F1 (Low). Both engines hold n and u differentiable across both faces.**
30:**F2 (Medium). The repair adds order flux concentrated on each face, which takes J beyond the list's "live functions", and the headline face objects include 
$ grep -n -o '\*\*Verdict: [A-Z ]*\*\*\|\*\*Findings:\*\* [a-z]*' /var/projects/_review_scratch/s12-build-r1-grok/grok_stdout.txt
3:**Verdict: CLEAR**
5:**Findings:** none
$ cat /var/projects/_review_scratch/s12-build-r1-grok/DONE
EXIT=0
```

## C3. The decision list on the live fields
```
$ sed -n '75,78p;104p;113p;117,118p' research/pde_ledger_v3/directives/S12_decision_list.md
1. **The source the committed law requires on order fields satisfying U1's support condition.** The object is `n Γ_B`, with `n`, `u`, `J_χ`,
   the two faces and the interior order all live functions. It covers each part of the domain the condition
   distinguishes: outside the region, inside it, and at its faces, including any part concentrated on a face if the
   order jumps there. No flow profile is imposed, and the river profile does not enter.
- **Freezes:** none intended. `J_χ` stays live (v2's `J_χ = 0` default is not adopted). Item 2 assumes fast ordering
- the interior order's law, including any dependence on history or strain; the faces' law; the form of `J_χ`;
- the flow profile (each force's requirement; S21);
- the slab's material law (S8/S22).
```

## C4. The engines' stated calculus, and the face-supported flux
```
$ sed -n '9,13p' research/pde_ledger_v3/scripts/S12_order_source.py
Distribution calculus retains the original differentiable n and u. J has
independent restrictions on all three open regions, independent one-sided
traces, and arbitrary distributions supported at the faces. Derivatives of
the supported distributions remain formal. No form, continuity, exterior
vanishing, boundary condition, or pointwise face value is selected for J.
$ sed -n '209,214p;227,228p' research/pde_ledger_v3/scripts/S12_order_source.py
    supported_flux = [[fun(*coordinates) for fun in row] for row in face_j_heads]
    distributed_flux = [
        s.Add(
            *(s.Mul(mask, row[c]) for mask, row in zip(regional_masks, regional_flux)),
            *(row[c] for row in supported_flux),
        )
    p.obj("flux/face_supported_fields", supported_flux)
    p.obj("flux/face_support_domains", list(zip(supported_flux, [s.Eq(w, face) for face in faces])))
$ sed -n '7,12p' research/pde_ledger_v3/mathematica/S12_order_source.wl
   Distribution convention: n and u retain their original differentiable
   bulk-field class. J has independent differentiable restrictions on the
   three open regions, independent one-sided traces, and an arbitrary
   distribution supported at each face. Derivatives of those distributions
   stay formal; no face-flux form, continuity, or boundary condition is imposed.
   The order can jump. No point value at a face is selected.
$ sed -n '117,120p' research/pde_ledger_v3/mathematica/S12_order_source.wl
  faceFlux = Map[(# @@ coordinates) &, faceFluxHeads, {2}];
  markers = Table[Unique["faceFluxMarker$"], {Length[faces]}];
  markerRules = Thread[markers -> ConstantArray[1, Length[faces]]];
  markedFlux = regionIndicators . regionFlux + markers . faceFlux;
$ grep -n '^flux/face_support_domains:' research/pde_ledger_v3/scripts/out/S12_order_source.out | cut -c1-200
25:flux/face_support_domains: [([j1_face1(t, x1, x2, x3, w), j2_face1(t, x1, x2, x3, w), j3_face1(t, x1, x2, x3, w), j4_face1(t, x1, x2, x3, w)], Eq(w, wminus(t, x1, x2, x3))), ([j1_face2(t, x1, x2, x
$ grep -n '^flux/face_support_domains:' research/pde_ledger_v3/mathematica/out/S12_order_source.out | cut -c1-200
25:flux/face_support_domains: {{{j1Face1[t, x1, x2, x3, w], j2Face1[t, x1, x2, x3, w], j3Face1[t, x1, x2, x3, w], j4Face1[t, x1, x2, x3, w]}, w == wm[t, x1, x2, x3]}, {{j1Face2[t, x1, x2, x3, w], j2Fa
```

## C5. The Grok leg's description of the repair
```
$ grep -o 'The diff opens[^.]*\.[^.]*\.' /var/projects/_review_scratch/s12-build-r1-grok/grok_stdout.txt
The diff opens `J_χ` into an independent differentiable restriction on each open region, with free one-sided traces, and adds an arbitrary face-supported flux whose divergence is left formal. The supplied balances, the two live faces, `n`, `u`, the interior order, the two-argument `G`, and the dimension assignment stay as they were.
```

## C6. The Claude leg's evidence files
```
present, 9487 bytes: /var/projects/_review_scratch/s12-build-r1-claude/derive_first.py
present, 37385 bytes: /var/projects/_review_scratch/s12-build-r1-claude/derive_first.stdout.txt
present, 2088 bytes: /var/projects/_review_scratch/s12-build-r1-claude/face_supported_flux_check.py
present, 2379 bytes: /var/projects/_review_scratch/s12-build-r1-claude/face_supported_flux_check.stdout.txt
$ sed -n '23p;47p' /var/projects/_review_scratch/s12-build-r1-claude/derive_first.stdout.txt | cut -c1-160
B_nu_cont_J_jump/balance1_face_m: 0
A_minus_B_face_m_with_B_fields_read_as_inside_traces: 0
$ grep -n '^coefficient_of_delta' /var/projects/_review_scratch/s12-build-r1-claude/face_supported_flux_check.stdout.txt | cut -c1-200
2:coefficient_of_delta: Derivative(K1(t, x1, x2, x3), x1) + Derivative(K2(t, x1, x2, x3), x2) + Derivative(K3(t, x1, x2, x3), x3)
3:coefficient_of_delta_prime: -K1(t, x1, x2, x3)*Derivative(a(t, x1, x2, x3), x1) - K2(t, x1, x2, x3)*Derivative(a(t, x1, x2, x3), x2) - K3(t, x1, x2, x3)*Derivative(a(t, x1, x2, x3), x3) + Kw(t, x1, 
4:coefficient_of_delta_double_prime: -L1(t, x1, x2, x3)*Derivative(a(t, x1, x2, x3), x1) - L2(t, x1, x2, x3)*Derivative(a(t, x1, x2, x3), x2) - L3(t, x1, x2, x3)*Derivative(a(t, x1, x2, x3), x3) + Lw(
```

## C7. Round 0: the Claude leg's note on n and u, and the repair wording
```
$ sed -n '36,37p' /var/projects/_review_scratch/s12-build-r0-claude/REVIEW.md
- I also checked whether jumps in `n` and `u` change anything. They do not: the face coefficient stays `f` times the
  interior-trace material flux. The gap is specific to `J_χ`.
$ grep -n -A3 'What must be true.\*\* The face objects' research/pde_ledger_v3/directives/_measurements/S12_build_review_r0_disposition.md
49:- **What must be true.** The face objects cover every face-concentrated part the committed law produces with `J_χ`'s
50-  face behaviour left open. That includes its jump across each face. Nothing about that behaviour is chosen: no
51-  continuity, no vanishing outside, no boundary condition. Every other object, field and label stays as it is.
52-- **Grok's CLEAR** does not bear on this. Its own derivation builds a single flux over the whole domain, the same
$ grep -n 'Choose nothing' _scratch/s12/build/resume_s12_repair_1.md
6:- The face objects cover every face-concentrated part the committed law produces when `J_χ`'s behaviour at each face is left open, including its jump across that face. Choose nothing about that behaviour: no continuity, no vanishing outside the region, no boundary condition.
```
