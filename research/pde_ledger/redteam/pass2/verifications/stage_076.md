---
unit_id: 076
batch: III.4
verdict: verified
overall: verified
material_change: false
verified_at: 2026-06-05
---

# Verification — unit 076 (pass 2)

**Verdict: verified — refresh-only (stale_output).** The single finding (F1) was a
`stale_output`: both committed transcripts predated their scripts and carried the
pre-renumber self-banner (`STAGE 59` / `STAGE 059`) while the current source already emits
`STAGE 076` (py:30, wl:26). No source edit; the math (n=5 wall-depth lock closed forms, the
`n=3` guard residual `3*K*rho**2/4`, all `= 0` residuals) is unchanged.

Orchestrator independent re-run regenerated both committed transcripts: SymPy exit 0,
Mathematica exit 0. The only committed-output change is the banner refresh
`STAGE 59/059 → STAGE 076`; all result/PASS/closed-form lines unchanged. Arbiter grep
confirms no residual self-epoch (59) banner/closing. `material_change: false`.
