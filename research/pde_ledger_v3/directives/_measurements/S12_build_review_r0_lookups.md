# S12 build review r0: mechanical lookups

Generator: `_scratch/s12/gen/build_r0_lookups.sh`. Verbatim retrieval, hashes and existence checks only.

## B1. Baseline, prompt and working tree
```
$ sha256sum -c _scratch/s12/build_review/s12_build_baseline_r0.sha256
research/pde_ledger_v3/scripts/S12_order_source.py: OK
research/pde_ledger_v3/mathematica/S12_order_source.wl: OK
research/pde_ledger_v3/scripts/out/S12_order_source.out: OK
research/pde_ledger_v3/mathematica/out/S12_order_source.out: OK
research/pde_ledger_v3/directives/_measurements/S12_build_run_record.md: OK
$ sha256sum _scratch/s12/build_review/s12_build_review_prompt_r0.md
b6377445b066dd17705fa3346063313cfe7d516bb325493da02721921163d98b  _scratch/s12/build_review/s12_build_review_prompt_r0.md
$ git status --porcelain | diff _scratch/s12/build_review/status_before_r0.txt -  (empty = unchanged)
1a2,3
> ?? research/pde_ledger_v3/directives/_measurements/S12_build_review_r0_disposition.md
> ?? research/pde_ledger_v3/directives/_measurements/S12_build_review_r0_lookups.md
diff exit: 1
```

## B2. The two legs' verdict lines
```
$ grep -n '^## Verdict' /var/projects/_review_scratch/s12-build-r0-claude/REVIEW.md
7:## Verdict: NEEDS REVISION
$ grep -n -o '\*\*Verdict: [A-Z ]*\*\*' /var/projects/_review_scratch/s12-build-r0-grok/grok_review.txt
1:**Verdict: CLEAR**
$ grep -n -o '\*\*Findings:\*\* [a-z]*' /var/projects/_review_scratch/s12-build-r0-grok/grok_review.txt
3:**Findings:** none
$ cat /var/projects/_review_scratch/s12-build-r0-grok/DONE
EXIT=0
```

## B3. What the decision list says about J_χ and the faces
```
$ grep -n 'J_χ' research/pde_ledger_v3/directives/S12_decision_list.md
32:∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B ,   Γ_B = Γ_return − Γ_drain ,   [Γ_B] = T⁻¹
36:  `χ_B` the ordered fraction, `J_χ` the order flux not carried by `u`. `Γ_B` is signed: positive is bulk material
75:1. **The source the committed law requires on order fields satisfying U1's support condition.** The object is `n Γ_B`, with `n`, `u`, `J_χ`,
104:- **Freezes:** none intended. `J_χ` stays live (v2's `J_χ = 0` default is not adopted). Item 2 assumes fast ordering
113:- the interior order's law, including any dependence on history or strain; the faces' law; the form of `J_χ`;
$ sed -n '75,78p' research/pde_ledger_v3/directives/S12_decision_list.md
1. **The source the committed law requires on order fields satisfying U1's support condition.** The object is `n Γ_B`, with `n`, `u`, `J_χ`,
   the two faces and the interior order all live functions. It covers each part of the domain the condition
   distinguishes: outside the region, inside it, and at its faces, including any part concentrated on a face if the
   order jumps there. No flow profile is imposed, and the river profile does not enter.
```

## B4. How each engine constructs J_χ (lines the Claude leg cites)
```
$ sed -n '9,13p;127,134p;146,149p' research/pde_ledger_v3/scripts/S12_order_source.py
Distribution calculus uses differentiable primitive bulk functions admitting
face traces. It permits an order jump; the face's point value is unspecified.
Additional singularities of the primitive bulk fields are outside this formal
calculus. No live field is frozen.
"""
    )
    p.emit(
        "calculus_scope",
        (
            "distributional order; differentiable bulk fields with face traces",
            "pointwise face values unspecified",
            "interior range and strict face order are supplied domain conditions",
        ),
    density = n_head(*coordinates)
    chi = chi_head(*coordinates)
    velocity = [fun(*coordinates) for fun in u_heads]
    flux = [fun(*coordinates) for fun in j_heads]
$ sed -n '6,10p;75,78p' research/pde_ledger_v3/mathematica/S12_order_source.wl
   expressions are differentiation, elimination, restriction, or extraction.
   Distribution convention: the live bulk functions admit differentiation and
   traces; the order can jump. No point value at a face is selected. Additional
   singularities of the primitive bulk fields are not modelled by this calculus.
   This script computes and prints; it has no expected-value acceptance guard. *)
  density = n @@ coordinates;
  velocity = (#[Sequence @@ coordinates] & /@ {u1, u2, u3, um});
  flux = (#[Sequence @@ coordinates] & /@ {j1, j2, j3, jm});
  fraction = chi @@ coordinates;
```

## B5. The Claude leg's evidence files for Finding 1 exist
```
present, 2189 bytes: /var/projects/_review_scratch/s12-build-r0-claude/face_jump_check.py
present, 1176 bytes: /var/projects/_review_scratch/s12-build-r0-claude/face_jump_check.stdout.txt
present, 13153 bytes: /var/projects/_review_scratch/s12-build-r0-claude/copies/py_P1_J_confined.py
present, 152164 bytes: /var/projects/_review_scratch/s12-build-r0-claude/out/py_P1_J_confined.stdout
present, 1128 bytes: /var/projects/_review_scratch/s12-build-r0-claude/out/py_P1_J_confined.extract.txt
present, 5209 bytes: /var/projects/_review_scratch/s12-build-r0-claude/derive_first.py
present, 4164 bytes: /var/projects/_review_scratch/s12-build-r0-claude/derive_first.stdout.txt
$ grep -n 'J_continuous\|J_jumps' /var/projects/_review_scratch/s12-build-r0-claude/face_jump_check.stdout.txt | cut -c1-200
2:J_continuous/balance2_face_coefficient: (n_in(t, x, a(t, x))*uw_in(t, x, a(t, x)) - n_in(t, x, a(t, x))*ux_in(t, x, a(t, x))*Derivative(a(t, x), x) - n_in(t, x, a(t, x))*Derivative(a(t, x), t))*f(t,
3:J_continuous/balance2_face_coefficient_minus_f_times_interior_mass_flux: 0
4:J_jumps/balance2_face_coefficient: (n_in(t, x, a(t, x))*uw_in(t, x, a(t, x)) - n_in(t, x, a(t, x))*ux_in(t, x, a(t, x))*Derivative(a(t, x), x) - n_in(t, x, a(t, x))*Derivative(a(t, x), t))*f(t, x, a
5:J_jumps/balance2_face_coefficient_minus_f_times_interior_mass_flux: Jw_in(t, x, a(t, x)) - Jw_out(t, x, a(t, x)) - Jx_in(t, x, a(t, x))*Derivative(a(t, x), x) + Jx_out(t, x, a(t, x))*Derivative(a(t,
```

## B6. How the Grok leg's own derivation builds the order flux
```
$ sed -n '14p;42p' /var/projects/_review_scratch/s12-build-r0-grok/derive_s12_object.py
No flow profile. J_chi, n, u live on the whole 4D domain.
Jw = F("J_w", ALL_ARGS)
```

## B7. The engines' flux lines, by line number
```
$ grep -n 'flux = ' research/pde_ledger_v3/scripts/S12_order_source.py research/pde_ledger_v3/mathematica/S12_order_source.wl
research/pde_ledger_v3/scripts/S12_order_source.py:149:    flux = [fun(*coordinates) for fun in j_heads]
research/pde_ledger_v3/mathematica/S12_order_source.wl:77:  flux = (#[Sequence @@ coordinates] & /@ {j1, j2, j3, jm});
```
