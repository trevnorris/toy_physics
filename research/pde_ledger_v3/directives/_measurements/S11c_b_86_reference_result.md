# S11c-b #86 — the corrected §1d variable-coefficient energy basis is 40 (RESULT)

## Headline (settled)
The §3a variable-coefficient energy basis, computed with the background scalars carrying their full jet tower
(the §1d-correct quotient), has dimension **40 = 10 uniform + 15 `∂W_bg`-spurion + 15 `∂μ_R,bg`-spurion**, and
it reduces to the engine's committed frozen basis **26** when the background Hessian is frozen. Both engines are
undercomplete at **8 spurion invariants/source**, but by **DIFFERENT mechanisms** (code-verified, not just
count agreement): **PY** enumerates 15 candidates and `basis_euler_signatures` freezes the Hessian (treats
`∂W_bg`/`∂μ_R,bg` as constant) ⇒ a spurious total-divergence quotient collapses 15→8; **WL**
(`mathematica/…_audit.wl` `newInvariantExpressions` L417 / `independentRepresentativeIndices` L599) HAND-CODES
only 8 invariant forms and takes literal linear independence — **no divergence quotient at all**, hence no
frozen-Hessian defect; it is undercomplete by *incomplete enumeration*. The correct per-source count is **15**
(8 → 15, **+14** total; carrier-independent). ⇒ the #89 repair differs per engine (PY: correct the quotient to
retain the Hessian; WL: complete the enumeration — since the correct quotient has nullity 0, linear independence
of all 15 gives 15). This is task #85's deferred §1d energy-basis quotient, root-caused and quantified.

## Established four independent ways (not one probe)
1. **Engine anchor (committed ground truth).** `grep S11CB_ENERGY_BASIS_COUNT scripts/out/S11c_b_brane_operator_sympy_audit.out` → `Integer(26)` (both anchorings; cross-engine agreed — WL emits 26 too). The corrected quotient MUST reduce to this in the frozen limit; that forces the 1-carrier uniform reading (below).
2. **Codex decision-leg** (`~/.s11_build/S11c_b_86_probes/codex_leg_probes/pre_directive_quotient_probe.stdout`): frozen combined `WM_NO_HESS rank 26`, full combined `rank 40`, `COMBINED_MINUS_SEPARATED 0` (no cross-family drops), per-source frozen `[1,4,5,6,7,9,10,13]`.
3. **Claude build-leg, from-scratch reconstruction** (imports nothing; `~/.s11_build/S11c_b_86_probes/claude_buildleg_probes/decomp.out`): the carrier table below; the **1-carrier** column is `FULL=40 FROZEN=26 uniform=10 spurionW=15 spurionMu=15`.
4. **Orchestrator crux check** (own CAS; `~/.s11_build/S11c_b_86_probes/orchestrator_crux_check.py`): a pair the frozen quotient merges (`c01+c02`) is genuinely independent once the Hessian is retained — `c01+c02 ≡ −H_ij u_i u_j (mod div)`, nonzero; at `H=0` it vanishes (the spurious frozen merge). Confirms the mechanism.

## The decisive decomposition — the §1d fix is carrier-independent; the reference's defect was not
`claude_buildleg_probes/decomp.out`:
```
[1-carrier {const}       ] combined FULL=40 FROZEN=26 | uniform=10 | spurionW=15 spurionMu=15   <- correct (engine)
[2-carrier {W_bg,μ_R_bg}  ] combined FULL=50 FROZEN=36 | uniform=22 | spurionW=15 spurionMu=15   <- the defective reference
[3-carrier {const,W,μ}    ] combined FULL=60 FROZEN=46 | uniform=32 | spurionW=15 spurionMu=15
```
The genuine §1d Hessian correction is **+14 in every column** (carrier-independent: spurion 8→15/source). The
**uniform coefficient-carrier count** is a separate axis adding **+10 per extra carrier** — it is NOT the #86
fix. The engine uses 1 carrier (one free coefficient per independent uniform invariant), giving frozen 26; the
frozen anchor therefore forces the 1-carrier reading, and the §1d-corrected basis is **40**.

## The reference SCRIPT was defective (removed) — and the two-leg build gate caught it
The Codex-built reference `scripts/S11c_b_energy_basis_reference.py` (rev-2 directive) enumerated the uniform
family with **two carriers** (each contraction carried separately with `W_bg` and `μ_R,bg`), giving uniform 22,
combined 50, frozen 36 ≠ engine 26. Cause: my rev-2 directive's "coefficient `W_bg` **or** `μ_R,bg`" was read as
per-source duplication (rule-15 note: my directive mis-specified the uniform family twice — rev-1 "functions of
them" ill-defined ring, rev-2 the carrier doubling). The **Grok build-leg missed it** (its own derivation shared
the 2-carrier misreading — the matching-number trap); the **Claude build-leg caught it** via a carrier-sensitivity
probe + the engine anchor. Both build legs confirmed the script's machinery is otherwise sound (Hessian retention
load-bearing, certificates genuine, no hard-coding, no asserts). Per the user's call the defective 2-carrier
script was removed (preserved as `~/.s11_build/S11c_b_86_probes/DEFECTIVE_2carrier_reference.{py,out}`); we accept
the triangulated 40 without rebuilding a 1-carrier reference (its independent-verification purpose is served by
route 3, the Claude leg's from-scratch 1-carrier reconstruction, which reduces to the engine anchor).

## Corrections to `S11c_b_coupling_84_basis_verification.md` (this record supersedes those two claims)
- ⛔ "The CORRECT variable-coefficient count is between 8 and 15 … rank 15 over-counts" — **WRONG**. Nullity is 0:
  the would-be null-Lagrangians `g·h` (`∂·h=0`) are the cross-product of gradients, which is parity-**odd** and
  excluded from the parity-even basis, so per-source **15 is exact**, not an over-count. Corrected total = 40.
- ⛔ The implied "second facet" (uniform and spurion quotiented in separate passes miss cross-sector merges) —
  **WRONG**. Combined vs family-separated is delta **0** (`COMBINED_MINUS_SEPARATED 0`, 1-carrier): the engine's
  separate-pass structure is fine. The sole defect is the frozen Hessian.
- ✓ STILL TRUE: the frozen §3a basis is undercomplete and the ~118-term bulk-core coupling residual is genuine
  §1d physics (PY reps vs WL reps differ at variable coefficients), not relabeling — but the mechanism is
  per-engine undercompleteness (PY: frozen Hessian; WL: incomplete hand-enumeration), per-source 8→15, and the
  corrected basis is 40.

## Scope note
#86 fixes the SPURION undercompleteness (PY frozen Hessian; WL incomplete hand-enumeration). The uniform family
inherits the engine's 1-carrier convention (the cross-engine-agreed frozen 26). Whether the uniform coefficient
should be a larger module (3-carrier → 60) is a SEPARATE question, out of #86's scope, with no current evidence —
both engines already agree on the 1-carrier 26.

## Next (NEXT = #88)
- **#88 (NEXT)** blast radius: does correcting 26→40 disturb the already-adjudicated families
  (admissibility/kinetic/advective/coupling — which may have agreed across engines only because BOTH were
  undercomplete at 8/source, by whichever mechanism)?
- #87 (largely resolved by the code-read above, not a pending blocker): WL is undercomplete by *incomplete
  hand-enumeration* — a DIFFERENT mechanism than PY's frozen quotient (both code-verified). Remaining CAS check:
  confirm WL's 8 hand-coded invariants span a strict subspace of the correct 15 (genuinely undercomplete, not 8
  different correct forms).
- #89 both-engine §3a repair — the two fixes DIFFER: **PY** retain the Hessian in the quotient; **WL** complete
  the enumeration to the full 15 (nullity 0 ⇒ linear independence of all 15 gives 15). Repaired engines must
  emit **40** (checked against this reference; 40 withheld from the builder).
- #90 PY §3c content fix (face+response) on the corrected basis.
