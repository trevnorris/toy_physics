# Independent review — S11c-a WL density-grounding fix DIRECTIVE

Artifact: `research/pde_ledger_v3/directives/S11c_a_wl_density_grounding_fix_directive.md`

## Commands

```bash
python3 research/pde_ledger_v3/_measurements/s11c_a_wl_density_grounding_fix_directive_review_grok/derive_density_normal_jet.py
```

Literal stdout: `derive_density_normal_jet.stdout`.
Engine path dump: `engine_density_paths.txt`.

## Verdict

Diagnosis of the ungrounded `rhoBulkBackground` free premise is correct; §3c forbids it; both §2b
members force `∂_w ρ_4D,bg⁰ = 0` by differentiation (script stdout). Object/property framing, leakage
controls, form control, keying requirement, and blindness instructions are sound.

One directive defect survives the filter: acceptance/out-of-scope under-specify the FACE_SHIFT-keyed
co-emits in the control-form and uniform-limit packages (see findings in the review return).
