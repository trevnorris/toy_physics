# S12 comparator decision list, review pass 1: mechanical lookups

Generator: `_scratch/s12/gen/cmp_dl_p1_lookups.sh`. Verbatim retrieval, hashes and existence checks only.

## E1. Baseline, prompt, working tree, exits
```
$ sha256sum -c _scratch/s12/cmp/cmp_dl_baseline_p1.sha256
research/pde_ledger_v3/directives/S12_comparator_decision_list.md: OK
$ sha256sum _scratch/s12/cmp/cmp_dl_review_prompt_p1.md
3588b049789086af303efb8b66a61028ad639eaba85e1baaf7a20ffb73942778  _scratch/s12/cmp/cmp_dl_review_prompt_p1.md
$ git status --porcelain | diff _scratch/s12/cmp/status_before_cmp_dl_p1.txt -  (empty = unchanged)
2a3,4
> ?? research/pde_ledger_v3/directives/_measurements/S12_comparator_decision_list_review_p1_disposition.md
> ?? research/pde_ledger_v3/directives/_measurements/S12_comparator_decision_list_review_p1_lookups.md
diff exit: 1
$ cat /var/projects/_review_scratch/s12-cmpdl-p1-codex/DONE /var/projects/_review_scratch/s12-cmpdl-p1-grok/DONE
EXIT=0
EXIT=0
```

## E2. Verdicts and finding headings
```
$ grep -n '^[0-9]\. \*\*' /var/projects/_review_scratch/s12-cmpdl-p1-codex/REVIEW.md | cut -c1-140
48:1. **High — Specify the face restriction when comparing a delta part. Raw coefficients are not unique distributional data.**
94:2. **High — Several functional forms do not ensure that G's two arguments are tested independently. The numerical fallback can lose an 
159:3. **Medium — Shared tag spelling does not establish that the payloads denote the same object on these streams.**
$ grep -n '^\*\*FINDINGS\*\*' /var/projects/_review_scratch/s12-cmpdl-p1-codex/REVIEW.md
46:**FINDINGS**
$ grep -n '^### [0-9]\.\|^\*\*FINDINGS\*\*' /var/projects/_review_scratch/s12-cmpdl-p1-grok/grok_stdout.txt | cut -c1-140
9:### 1. High. Item 3’s delta coefficient is not the value on the face, so equal distributions fail on the printed domain
53:### 2. High. Item 6 does not catch a comparator that only evaluates inside the supplied domain, so the step is hidden and every required 
84:**FINDINGS**
```

## E3. The engines distinguish a raw delta coefficient from its trace on the face
```
$ sed -n '347,349p;363,364p' research/pde_ledger_v3/scripts/S12_order_source.py
            coefficient = source.coeff(delta)
            jacobian = s.Abs(s.diff(delta.args[0], w))
            trace = s.expand(coefficient.subs(w, face) / jacobian)
            p.obj(name + "/coefficient_operand", coefficient)
            p.obj(name + "/trace_coefficient", trace)
$ sed -n '266,267p;275,276p' research/pde_ledger_v3/mathematica/S12_order_source.wl
      coefficient = Coefficient[source, delta];
      trace = Expand[coefficient /. w -> faces[[k]]];
      object[label <> "/face/" <> ToString[k] <> "/coefficient_operand", coefficient];
      object[label <> "/face/" <> ToString[k] <> "/trace_coefficient", trace];
```

## E4. Shared tag names, differently constituted payloads
```
$ sed -n '151,154p;260p;323p' research/pde_ledger_v3/scripts/S12_order_source.py
    material_current = (density,) + tuple(density * v for v in velocity)
    order_current = (chi * density,) + tuple(
        chi * density * v + j for v, j in zip(velocity, flux)
    )
    p.obj("material/conserved_current_operand", mass_current_field)
        p.obj(label + "/material_elimination_operand", bulk_multiple)
$ sed -n '100,101p;181,182p;242p' research/pde_ledger_v3/mathematica/S12_order_source.wl
  materialDensity = density; materialCurrent = density velocity;
  orderDensity = fraction density; orderCurrent = fraction density velocity + flux;
  object["material/conserved_density_operand", materialDensityField];
  object["material/conserved_current_operand", materialCurrentField];
    object[label <> "/material_elimination_operand", timeRules];
```

## E5. The legs' evidence (existence, and the lines each finding quotes)
```
present, 2597 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-codex/06_review_quotes.py
present, 30877 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-codex/06_review_quotes.stdout.txt
present, 5404 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-codex/07_distribution_and_G_diagnostics.py
present, 4033 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-codex/07_distribution_and_G_diagnostics.stdout.txt
present, 1449 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-codex/08_admissibility.py
present, 686 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-codex/08_admissibility.stdout.txt
present, 13280 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-grok/scripts/18_distribution_residuals.py
present, 8337 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-grok/evidence/18_distribution_residuals.txt
present, 2384 bytes: /var/projects/_review_scratch/s12-cmpdl-p1-grok/evidence/19_control_gap_and_dimensions.txt
$ grep -n 'partial_1-minus-partial_2: \|^[a-z_]*_sum: G=\|admissible separating' /var/projects/_review_scratch/s12-cmpdl-p1-codex/07_distribution_and_G_diagnostics.stdout.txt | cut -c1-120
21:rational_sum: G=(a + b)/(a + b + 1)
23:exponential_sum: G=1 - exp(-a - b)
25:hyperbolic_sum: G=tanh(a + b)**2
27:zero_face_rational_sum: G=(a + b - 3)**2/((a + b - 3)**2 + 1)
30:zero_face_exponential_sum: G=1 - exp(-(a + b - 3)**2)
33:zero_face_hyperbolic_sum: G=tanh(a + b - 3)**2
36:admissible separating two-offset profile: (a + 2*b)/(a + 2*b + 1)
37:partial_1-minus-partial_2: -1/(a + 2*b + 1)**2
$ grep -n 'components=' /var/projects/_review_scratch/s12-cmpdl-p1-codex/06_review_quotes.stdout.txt
272:scripts line=362 components=5
274:mathematica line=78 components=4
277:scripts line=2059 components=5
279:mathematica line=180 components=4
282:scripts line=2057 components=3
284:mathematica line=176 components=3
```
