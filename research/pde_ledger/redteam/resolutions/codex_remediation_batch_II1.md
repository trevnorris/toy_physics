You are remediating two blocked findings from batch II.1 v2 fix_loop: **stage 024 F1 and F2**, both blocked on a sign-convention question.

# Background

In the original directive `redteam/directives/stage_024.md`:
- **F1** (mathematica_transliteration) prescribed replacing Section III/V of `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl` with a matrix-inverse-based independent re-derivation using `mPair = {{omegaU^2 - omega^2, +rPair}, {+rPair, omegaW^2 - omega^2}}`.
- **F2** (insufficient_verification) prescribed adding an analogous SymPy anchor at `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:240-286` using the same matrix.

You correctly identified that the directive's positive off-diagonal `+rPair` is inconsistent with the paper's `Q_r = G_U²Ω_W² + 2 G_U G_W R + G_W²Ω_U²` (positive mixed term). Inverting the `+rPair` matrix produces a `-2 g_U g_W R` mixed term in `g^T M^{-1} g`, not the paper's `+2`. You blocked both findings pending sign resolution.

# Resolution (math-authority decision: yours to verify)

**The correct conservative matrix has `-R` off-diagonal.** Independent derivation:

For the paper's spring-matrix definitions (`paper/stages/stage_024.tex:103,108,113`):
- `Δ_r = Ω_U² Ω_W² − R²` (det of [[Ω_U², ±R], [±R, Ω_W²]] is the same)
- `Q_r = G_U² Ω_W² + 2 G_U G_W R + G_W² Ω_U²` (positive 2 G_U G_W R)
- `P_r = Ω_U² G_W + R G_U` (positive R G_U)

With `M = [[Ω_U² − ω², −R], [−R, Ω_W² − ω²]]`:
- `det M = (Ω_U² − ω²)(Ω_W² − ω²) − R²` ✓
- `M^{-1} = (1/det) × [[Ω_W² − ω², +R], [+R, Ω_U² − ω²]]`
- `g^T M^{-1} g = (g_U²(Ω_W² − ω²) + 2 g_U g_W R + g_W²(Ω_U² − ω²)) / det` — at ω=0 gives `Q_r/Δ_r` ✓
- `(M^{-1} g)_W = (+R g_U + (Ω_U² − ω²) g_W) / det` — at ω=0 gives `(R G_U + Ω_U² G_W)/Δ_r = P_r/Δ_r` ✓

The paper's Lagrangian (search `paper/stages/stage_024.tex` and the part-level appendix for the moving-throat conservative spring matrix or `+R U W` Lagrangian term) is the upstream authority. Confirm by reading. Then proceed.

# Authorized edits

**Mathematica** (`mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`):
- Replace Section III (lines 113-150 or whatever the III block currently spans) per the F1 directive's prescription, BUT with `mPair = {{omegaU^2 - omega^2, -rPair}, {-rPair, omegaW^2 - omega^2}}` (negative off-diagonal). Keep `coupling = {gU, gW}`, `qRef = gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2` (paper's positive form), `pRef = omegaU^2*gW + rPair*gU` (paper's positive form), `deltaRef = omegaU^2*omegaW^2 - rPair^2`. The `expectZero` anchor checks for `zFromMatrix - zRefRational` and `nFromMatrix - nRefRational` should now succeed.
- Section V replacement (if any) follows the same sign convention.

**SymPy** (`scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:240-286`):
- Insert the analogous F2 anchor block using `Mpair = sp.Matrix([[OmegaU**2 - omega**2, -Rmix], [-Rmix, OmegaW**2 - omega**2]])` with `g = sp.Matrix([gU, gW])`. Verify `g.T * Mpair.inv() * g` at ω=0 matches paper's `Q_r/Δ_r` form.

# Process

1. Read `paper/stages/stage_024.tex` lines ~80-130 to find the upstream Lagrangian / matrix definition. Confirm the sign convention. Cite the paper line you used to make the call.
2. Apply F1 (Mathematica Section III) and F2 (SymPy anchor) with the `-R` off-diagonal.
3. Run both scripts to confirm exit 0 — RESPECT the Mathematica single-seat rule (one `math -script` at a time across the whole system; the orchestrator is NOT running any other `math -script` right now, but DO check before invoking).
4. Update `redteam/directives/stage_024.md`: replace the `## Blocked: F1` and `## Blocked: F2` blocks with `## Applied: F1` and `## Applied: F2` blocks (preserving the `## Blocked` text as-is by appending to the directive — keep the historical record, just add Applied blocks below).
5. Report each block's resolution as `files_changed`, `summary`, `deviation: <how the sign convention was resolved>`.

# Critical rules

- **Read paper first.** Don't apply blindly; confirm the `-R` sign by reading the paper section.
- **Single-seat Mathematica.** One `math -script` at a time.
- **No JSON.** Markdown/YAML only.
- **No fake commentary scripts.** Don't write python3 -c commentary blocks; read and reason.
- **You're the math authority for this remediation.** If your independent derivation says the directive's sign convention is wrong, override it with the right one and document the override.

# Performance pitfall — Section IV in current 024 mathematica hangs

Important: the orchestrator's most recent `math -script` run of `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl` hit a hang inside `Table[tripleOverlap[basis[[i]], qMat, basis[[j]]], {i, 1, 5}, {j, 1, 5}]` at Section IV (line 202). 25 × 729 calls to a non-memoized `i6` sphere integral, each wrapped in `FullSimplify[..., Assumptions -> $Assumptions]`, blew up at >18 min CPU with no output. The v1 script (without F3/F4) ran in seconds; symbol context contamination from earlier sections is the likely cause.

When you do the Section III rewrite for F1, also add at the very top of Section IV (line 199 region, just after `banner["SECTION IV — AXISYMMETRIC SPLITTING MATRIX"];` and before line 200's `$Assumptions = True;`):

```
ClearAll[gU, gW, rPair, omegaU, omegaW, mPair, zFromMatrix, nFromMatrix,
  qRef, hRef, pRef, deltaRef, sRef, zRefRational, nRefRational,
  zSeries, z0, z2, z4, nSeries, n0, n2, n4, dSeries, d0, d2, d4,
  gUiso, gWiso, rrIso, omU, omW];
```

(or whichever symbols actually leaked — the exact set may depend on your Section III rewrite). Then `$Assumptions = True;` and proceed. Add a one-line comment: `(* Reset symbol context before sphere-integral table — F3/F4 additions leaked symbols. *)`.

Optionally also memoize `i6[i_, j_, k_, l_, m_, nn_] := i6[i,j,k,l,m,nn] = Integrate[...]` so re-derivation runs don't pay 18000-call cost. Same for `i4`. Add a `(* memoized *)` comment.

Confirm by running the modified script and observing Section IV completes in seconds, not minutes.

# Working directory

`/var/projects/toy_physics/research/pde_ledger`
