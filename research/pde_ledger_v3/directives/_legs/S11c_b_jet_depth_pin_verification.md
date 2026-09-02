# Verification — S11c-b jet-depth reconciliation + spec-pin (B) records vs committed reality + governing text

You are verifying that an orchestrator's committed reconciliation + spec-pin is CORRECT and honestly recorded — no
overclaim, no inaccuracy, no missing caveat — AND that the decisive governing citations the pin rests on actually
say what is claimed. This is a correctness/faithfulness check against git + the source spec, not a re-run of the
physics adjudication (that was leg-gated). Cite file:line / git evidence for every claim. Commit under review:
`be743891` (`git show --stat be743891`; `git log --oneline -6`).

## Background (what was decided)

The #89b PY-check flagged a strong-row jet-depth mismatch (WL order-3 vs PY order-2). It was reconciled as a
REPRESENTATION mismatch (not a jet-depth freeze, not a physics error) and spec-PINNED (B): the slab U-momentum row
is the constraint-reduced in-plane equation carrying `−∇μ_θ`; PY (currently raw held-fixed `δU/δu`) is the engine
that must change; WL is already correct. Records:
- `research/pde_ledger_v3/directives/_measurements/S11c_b_strong_row_jet_depth_reconciliation.md`
- `research/pde_ledger_v3/directives/S11c_b_jet_depth_spec_pin_decision.md` (FOLDED VERDICT section)
- `STATUS.md` top section; advisor/leg computational evidence under
  `…/_measurements/s11c_b_jet_depth_consult_{codex,grok}/` (raw `.log` transcripts are gitignored, in
  `~/.s11_build/S11c_b_jet_depth/`).

## What to verify

1. **The pin rests on real governing text — re-read it yourself and confirm or refute.** Quote from the actual
   source, do not trust the records' paraphrase:
   - `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` around `:325-352` and `:420-430` — does it say the
     virtual displacements `{δ_vθ, δ_v(δW), δ_vu}` are NOT independent, "Do NOT vary U with θ held fixed", and that
     the in-plane equation "must carry `−∇(δU/δθ)` … varying at fixed θ removes this contribution … selects the
     convention uniquely"? If it does, (B) is supported and (A) is refuted. If it does NOT, say so — that would
     overturn the pin.
   - `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md:124-148,272-288` — is §1c's "θ may not be
     eliminated through a constraint before this derivative" genuinely scoped to the CONSTITUTIVE `μ_θ` (not the
     U-row), and does §3b invoke S11b's balance-law method / binding virtual-displacement rule? Is the separate
     `S11CB_MU_THETA_OPERATOR` a constitutive operand (face affinity) rather than an independent `θ=0` EOM?
   - Spot-check the corroboration `S11b_interface_coupling_law_mathematica_audit.wl:~280`
     (`constrainedUL = Expand[explicitUL + I k muTheta]` or equivalent) — does the S11b engine's in-plane EOM carry
     the `μ_θ` gradient?
2. **The records match the committed reality.** Do the reconciliation record, decision list, and STATUS accurately
   describe the pin, the mechanism (order-3 = `∇μ_θ` constraint reaction, raw EL order-2 on both engines), and the
   engine consequence (PY changes: fold constraint + raise `STRONG_ROW_JET_DEPTH` 2→3 + keep `μ_θ` separate +
   θ-row = mass evolution; WL unchanged; §3b amendment; #88 re-adjudication)? Any overclaim or misstatement? Do the
   committed advisor/probe scripts + stdout support the "three consistent computations" claim (spot-check a stdout)?
3. **Housekeeping.** Did `be743891` accidentally include unrelated work (e.g. the other session's `memory/*`)?
   (`git show --stat be743891`.) Are the moved-`.log` citations consistent (records point to `~/.s11_build/…` for
   the gitignored raw transcripts, not to a committed path that does not exist)? Is `STATUS.md`'s new section
   consistent with the reconciliation record and the memory?

## Output

A numbered list of any INACCURACY, OVERCLAIM, MISSING CAVEAT, or — most important — any way the pin (B) is NOT
actually supported by the governing text, each as `severity — file:line — what the record says — what the source /
git actually shows`. Severity ∈ {MUST-FIX, SHOULD-FIX, NIT}. Then `VERDICT: N issues (M must-fix)` or
`VERDICT: RECORDS FAITHFUL AND PIN (B) SUPPORTED`. If you find the S11b text does NOT support (B), that is a
MUST-FIX and the single most valuable thing you can report.
