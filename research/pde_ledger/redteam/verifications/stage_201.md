---
unit_id: 201
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T03:30:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 201

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the previously-absent second engine at the directive-registered path
`mathematica/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_mathematica_audit.wl`
(new untracked file, 306 lines, mtime 2026-06-02 00:36:08). No other script was touched — the SymPy
reference `scripts/...stage201..._sympy_audit.py` is unchanged (the directive forbade touching it; the
empty `stage_201_diff.patch` and absent `stage_201_sympy.log` are consistent with a single new `.wl` and
no SymPy edit). The runner regenerated
`mathematica/output/...stage201..._mathematica_audit.txt` (mtime 2026-06-02 03:08:11, newer than the
`.wl`), which exits 0 with 22 PASS lines and 0 FAIL lines.

**Assessment:**
Each claim-manifest item M1–M8 appears as an explicit, non-tautological residual/PASS check, with
residual vectors printed as `MatrixForm[{{0,...,0}}]` before each PASS:

- M1 (lines 126–151) — right-section identity `M_* . S - I_3 = 0` and cancellation `M_* . Δx_rep + q = 0`. (2 PASS)
- M2 (lines 153–167) — three mismatch-chart identities via `Exp[...]`/`PowerExpand`. (3 PASS)
- M3 (lines 169–187) — `Δx_rep - closedform = 0`, free rows {1,2,3,4,6} vanish, dependent entries (Kη,μ,T) match. (3 PASS)
- M4 (lines 189–201) — `Δx_rep = 0` at R_tr=R_nt=R_eta=1, and eta-only perturbation yields nonzero Kη repair. (2 PASS)
- M5 (lines 203–227) — componentwise `Log[x_proj/x] - Δx_rep = 0` and fixed-point reduction `x_proj - x = 0` at target ratios. (2 PASS)
- M6 (lines 229–259) — `Det[block] - (1+chi0) = 0`, det ≠ 0, unique solution, dT/dKη/dμ match closed form, residual `M_*·Δx_dep + q = 0`. (6 PASS)
- M7 (lines 261–280) — pairwise witness packet equals intrinsic packet under the witness substitution. (1 PASS)
- M8 (lines 282–296) — `M_* . (Δx + Δx_rep^lin) = 0` and `q^lin - M_*·Δx = 0`. (2 PASS)

Total 21 manifest checks + the bare `Δx_rep at target = 0` print already counted; the `grep -c '^PASS:'`
yields 22 (M1:2, M2:3, M3:3, M4:2, M5:2, M6:6, M7:1, M8:2 = 21, plus the M4 `expectTrue` boolean PASS
also emitting a PASS line). Every manifest item M1–M8 is represented; none passes by collapsing the
object under test.

Non-tautology spot-checks of the traps the directive flagged:
- **M8 / Δx not collapsed:** `generalDrift = {Dl,Dc,Dg,DU,DKeta,DW,Dmu,DT}` (line 284) keeps all eight
  independent drift symbols. Output line 104 shows `q^lin` depending on all eight, and line 106 shows
  `Δx_rep^lin` is a full 8-vector function of all eight symbols. The residual `M_*·(Δx+Δx_rep^lin)=0`
  therefore vanishes only because the section is a genuine right-inverse — not by construction. Trap avoided.
- **M3 support:** `dxRep[[freeSlots]]` (line 183) prints the actual {1,2,3,4,6} entries and shows them
  zero, rather than asserting support by fiat.
- **A "definitional" exp∘log identity** (the SymPy A2 the auditor flagged as harmless bookkeeping) is NOT
  carried over as a load-bearing M2 check — M2 here checks `m_* - R_*^exponent` directly.

**Independent-derivation check (load-bearing, not a rubber-stamp):**
The `.wl` is a GENUINELY INDEPENDENT derivation, not a transliteration of the `.py`. A reviewer could not
line-pair the two files. Corresponding sections by a DIFFERENT route:

1. **Section / repair construction (M1) — different decomposition.**
   - `.py` 66–74: `Pdep = Matrix.hstack(...); Pdep_inv = Pdep.inv(); Sdep = Edep * Pdep_inv` — inverts the
     3×3 pivot block and left-multiplies by an embedding matrix.
   - `.wl` 121–144: builds an 8×8 `constrainedSystem = Join[Mstar, freeSelector]` (M_* stacked with five
     unit-vector rows pinning the free coordinates to zero) and obtains both `Ssection` and `dxRep` via
     `LinearSolve[constrainedSystem, ...]`. No pivot-block inverse is formed. This is the directive's
     prescribed independent route, satisfied.

2. **Dependent-triple solve (M6) — directly on the linear system, NOT by reusing S.**
   - `.py` 162–170: `sp.solve([Eq(dep_equations[i],0)...], [dT,dKeta,dmu])`.
   - `.wl` 234–238: `Solve[Thread[Mstar . dxDep == -q], {dT,dKeta,dmu}, Reals]`. The `.wl` M6 does not
     touch `Ssection`; it solves the raw system and independently confirms `Det[block]-(1+chi0)=0` and
     uniqueness. Matches the directive's "must be solved directly, not by reusing the section S."

3. **Log-ratio identities (M2/M5) — exponent arithmetic vs expand_log.**
   - `.py` `expect_zero`/`simplify_log` (21–38) route everything through `sp.expand_log(..., force=True)`.
   - `.wl` `cleanScalar` (29–32) routes through `FullSimplify[PowerExpand[Together[Expand[...]]]]` —
     `PowerExpand` over explicit exponents, the directive's prescribed alternative, not a log-expansion mirror.

Variable choreography differs throughout: `.py` uses `Sdep/Pdep/Edep`, 0-based indices `T_idx=7,
Keta_idx=4, mu_idx=6`; `.wl` uses `Ssection/constrainedSystem/freeSelector`, 1-based `depSlots={8,5,7}`,
`AssociationThread` row construction. The shared `M_*` matrix literal is explicitly permitted by the
directive ("you may state the same matrix, but at least one load-bearing result must come from an
independent route") and the independent-route condition is met by M1 and M6. Classification: **independent
derivation confirmed**, not transliteration.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy log was captured and `stage_201_diff.patch` is empty — both expected, since
the directive forbade touching the SymPy script and only added a `.wl`. The SymPy `.py` reads unchanged.

**Mathematica:** exit=0. Notable lines from the captured log:
- `PASS: M1 residual: M_* . S - I_3` and `PASS: M1 residual: M_* . Delta x_rep + q` (residuals printed as
  the 3×3 / 1×3 zero matrices).
- `PASS: M6 residual: Det[dependent block] - (1+chi0_star)` with `dependent solution = {dT -> ...,
  dKeta -> Log[Reta], dmu -> ...}` matching the closed form.
- `PASS: M8 residual: M_* . (Delta x + Delta x_rep^lin)` with `q^lin` and `Δx_rep^lin` printed as full
  functions of all eight drift symbols (confirms the non-tautology).
- Final `# exit_code: 0`.

**Output freshness:** confirmed. Output `.txt` mtime 2026-06-02 03:08:11 is newer than the `.wl` mtime
2026-06-02 00:36:08 — the saved output was regenerated by the orchestrator's post-fix re-run.

## Material-change assessment

`material_change`: false. This finding only ADDS a second engine; it changes no derived result. The SymPy
reference is byte-for-byte unchanged, and the new `.wl` independently reproduces the same repair vector,
mismatch chart, orbit projection, and linearized compiler that the SymPy script already established and
that downstream units consume. No downstream unit's inputs changed; no `upstream_stale` propagation is
warranted on math grounds.

## Side observations (non-blocking)

- The `.wl` banner and ledger text read "STAGE 201"; the SymPy `.py` banner/ledger still read "STAGE 184"
  (lines 41, 263). This is a pre-existing cosmetic mislabel in the SymPy script (the file is the stage-201
  audit). It is outside this finding's scope, does not affect any check, and the directive forbade editing
  the `.py` — flagging only, not blocking.
- M4's eta-only-perturbation check is phrased as `Exp[dxRep[[5]] /. {...Reta->r}] != 1` with `r != 1`
  asserted in `$Assumptions`; this is a legitimate "nonzero under perturbation" witness, not a tautology.

## Verdict justification

The single finding F1 (missing Mathematica engine under the dual-engine rule) is fully resolved: the `.wl`
exists at the registered path, exits 0 with 22 PASS and 0 FAIL, and every manifest item M1–M8 is an
explicit non-tautological residual check. Critically, the `.wl` is an independent derivation — it obtains
the section and repair vector by a constrained `LinearSolve` (not pivot-block inversion + multiplication),
solves the dependent triple directly with `Solve` (not by reusing the section), and verifies the log-ratio
identities via `PowerExpand` (not `expand_log` mirroring) — so a reviewer could not line-pair it against
the SymPy script. No regressions, no material change. Verdict: verified.
