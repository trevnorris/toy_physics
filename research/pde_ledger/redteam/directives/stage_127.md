---
unit_id: 127
batch: IV.4
created_at: 2026-05-27T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 127

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_127.tex:1` quote: `\section[Stage 127]{Stage 127: Geometric Mouth-Penetration Families}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage127_penetration_families.md:2` quote: `# Moving-Throat PDE — Stage 229: Geometric Mouth-Penetration Families` (note: notes file title says "Stage 229" rather than "Stage 127" — yet another label inconsistency in this unit's documentation chain)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage127_penetration_families_sympy_audit.py:12` quote: `banner("STAGE 110 — GEOMETRIC MOUTH-PENETRATION FAMILIES")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:26` quote: `banner["STAGE 110 — GEOMETRIC MOUTH-PENETRATION FAMILIES"];`

## Resolve before fix_loop

The script banner strings say "STAGE 110", the paper card and Mathematica end-message say "Stage 127", and the source notes file title says "Stage 229". Three different stage numbers reference the same unit. Which number is canonical for this stage?

Possible directions (the user picks one):
- (a) The paper card and filename are authoritative (Stage 127) → update both script banners at SymPy L12 and Mathematica L26 to "STAGE 127 — GEOMETRIC MOUTH-PENETRATION FAMILIES"; also update notes title to "Stage 127". No numerical change.
- (b) The notes file's "Stage 229" is canonical (this stage was re-numbered) → update paper card title, filename references, and both script banners to "Stage 229" — large refactor; would normally require renaming files. Almost certainly NOT the intended direction.
- (c) The "STAGE 110" banner came from a sibling file by copy-paste and is purely clerical → same outcome as (a).

Recommended direction: (a) or (c) — both lead to setting all banners to "STAGE 127". The notes title's "Stage 229" appears to be an old internal numbering left over from before the renumbering; the user should confirm and (optionally) also clean up the notes title in a separate paper-side patch.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:26–48`

**Issue:**
The Mathematica script is a structural line-by-line port of the SymPy script. Both engines take the closed-form bias factors `g_slab(x) = (2/(pi x)) sin(pi x / 2)` and `g_exp(x) = 2(2 + pi x e^{-1/x}) / ((4 + pi^2 x^2)(1 - e^{-1/x}))` from the notes verbatim, without independently deriving them from the source integrals `int_0^{xL} (1/(xL)) cos(pi z / (2 L)) dz` and `int_0^L (e^{-z/(xL)} / (xL (1 - e^{-1/x}))) cos(pi z / (2 L)) dz`. A typo or sign error in either closed form would slip through both audits identically.

**Required change:**

Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl`. Between the existing line 29 (`$Assumptions = Element[x, Reals] && x > 0;`) and line 31 (`rDisc = Sqrt[4107 - 100*Pi^2];`), insert an independent symbolic-integration check.

Before (current L28–35):
```
Clear[x];
$Assumptions = Element[x, Reals] && x > 0;

rDisc = Sqrt[4107 - 100*Pi^2];
gMinus = N[(2*rDisc - 37*Sqrt[3])/(20*Pi), 80];

gSlab = FullSimplify[2*Sin[Pi*x/2]/(Pi*x), Assumptions -> $Assumptions];
gExp = FullSimplify[2*(2 + Pi*x*Exp[-1/x])/((4 + Pi^2*x^2)*(1 - Exp[-1/x])), Assumptions -> $Assumptions];
```

After:
```
Clear[x, z];
$Assumptions = Element[x, Reals] && x > 0 && Element[z, Reals];

rDisc = Sqrt[4107 - 100*Pi^2];
gMinus = N[(2*rDisc - 37*Sqrt[3])/(20*Pi), 80];

(* Independent derivation of g_slab from the slab source integral.
   With L = 1 (the bias factor is L-independent after the change of variables u = z/L). *)
gSlabFromIntegral = FullSimplify[Integrate[(1/x)*Cos[Pi*z/2], {z, 0, x}], Assumptions -> x > 0];
gSlab = FullSimplify[2*Sin[Pi*x/2]/(Pi*x), Assumptions -> $Assumptions];
slabClosedFormResidual = FullSimplify[gSlabFromIntegral - gSlab, Assumptions -> x > 0];
If[TrueQ[slabClosedFormResidual === 0],
   pass["slab closed-form matches source integral"],
   fail["slab closed-form matches source integral", slabClosedFormResidual]];

(* Independent derivation of g_exp from the exponential source integral. *)
gExpFromIntegral = FullSimplify[Integrate[(Exp[-z/x]/(x*(1 - Exp[-1/x])))*Cos[Pi*z/2], {z, 0, 1}], Assumptions -> x > 0];
gExp = FullSimplify[2*(2 + Pi*x*Exp[-1/x])/((4 + Pi^2*x^2)*(1 - Exp[-1/x])), Assumptions -> $Assumptions];
expClosedFormResidual = FullSimplify[gExpFromIntegral - gExp, Assumptions -> x > 0];
If[TrueQ[expClosedFormResidual === 0],
   pass["exp closed-form matches source integral"],
   fail["exp closed-form matches source integral", expClosedFormResidual]];
```

(The `Clear[x, z]` is needed because the new code introduces `z` as a fresh symbol; `$Assumptions` is updated to allow `z` to be a real integration variable.)

The rest of the script (lines 40–53) is unchanged: same `FindRoot` calls, same `expectApprox` residual checks, same `Stage 127 Mathematica audit passed.` message.

**Claim manifest** (not a missing-script finding, but documenting what the new checks verify):

M1: `Integrate[(1/x)*Cos[Pi*z/2], {z, 0, x}] === 2*Sin[Pi*x/2]/(Pi*x)`, i.e., the closed-form `g_slab(x)` is the true value of the slab source integral.

M2: `Integrate[(Exp[-z/x]/(x*(1 - Exp[-1/x])))*Cos[Pi*z/2], {z, 0, 1}] === 2*(2 + Pi*x*Exp[-1/x])/((4 + Pi^2*x^2)*(1 - Exp[-1/x]))`, i.e., the closed-form `g_exp(x)` is the true value of the exponential source integral.

Both integrals are evaluable in closed form by Mathematica (slab: standard `Cos` antiderivative; exp: standard `Exp * Cos` antiderivative with finite-interval evaluation). The `FullSimplify` step is needed to compare against the algebraically-rearranged form in the notes.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 127` and confirm:
1. Two new `PASS:` lines appear in the output: `PASS: slab closed-form matches source integral` and `PASS: exp closed-form matches source integral`.
2. The existing `PASS: slab compensation root` and `PASS: exponential compensation root` lines remain.
3. The script exits 0.
4. SymPy script is unchanged (the transliteration fix applies only to the second engine).
