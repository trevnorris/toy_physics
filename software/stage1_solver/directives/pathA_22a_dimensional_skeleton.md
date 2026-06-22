# Directive pathA_22a — Dimensional homogeneity + dimensionless-knob ledger (scopes the m̂0²·S_port verdict)

**Status:** DRAFT v2 (post Codex design-review, SOUND-WITH-FIXES → 10 fixes applied) — for Codex CONFIRM-PASS (read-only) →
execute. **Owner split:** Codex derives + codes + iterates; Claude reviews. Orchestrator owns this directive + decisions.
**Gate context:** FIRST, CHEAP step of the pathA_22 critical path (user-gated 2026-06-21: "dimensional skeleton first"). It is a
SCOPING gate that runs BEFORE the hard derivations (emergent `G`, mass-bridge, `W_eff` kernel, source factor `m̂0`, port
convention `S_port`/`χ_Q`). Resume: `decisions/13_emergent_constants_derivation.md` §0 (5g) + §1–§5.

**⚠ HONEST SCOPE (Codex design-review #2 — do NOT overclaim):** dimensional analysis ALONE CANNOT prove the comparison is a pure
number (it cannot rule out arbitrary dimensionless functions `f(c_γ/c_s, W_eff/a, R0/a, χ_Q, …)`). So this step does NOT, and
must NOT promise to, deliver a clean `FAIL_ABLE` verdict. It delivers: (i) a dimensional HOMOGENEITY result, (ii) a resolved
question of whether `P0=N0/D0` is truly dimensionless or hides a radiation-normalization factor, (iii) a LEDGER of the
dimensionless knobs/residuals the eventual verdict hinges on — each classified — and (iv) the explicit catch of any free knob
ALREADY visible in the symbolic skeleton. The decisive next step it scopes is deriving the MINIMAL COMBINATION
`ξ := m̂0²·S_port / [G·c_s⁵/(a⁵·c⁵)]` (NOT a broad symbolic exercise, and NOT necessarily `G`/`m̂0`/`S_port` separately).

---

## 1. The question
With physical units RESTORED (do NOT trust the `a = c_s = ħ = m = 1` pins), and on the model's own normalization law
```
code:   R_norm = m̂0²·S_port·P0  −  54·G·c_s⁵/(5·a⁵·c⁵)         (patha_extraction.py:544)
paper:  m̂²·χ_Q·P0 = 54·G·c_s⁵/(5·a⁵·c⁵)                        (pde.tex eq:outgoing-factorized-normalization, ~line 2080)
```
1. Is the comparison dimensionally HOMOGENEOUS once units are restored (catching any dropped `a⁵`/volume-Jacobian/`ρ₀` power)?
2. Which dimensionless KNOBS/residuals does `R_norm = 0` depend on, and which of them are (a) fixed by a prior derivation,
   (b) branch-determined by a target-blind solve, (c) underived residuals, or (d) TRUE FREE calibration knobs (= tunability
   channels that would make the GR-quadrupole "match" not a real external prediction)?

This is [[feedback-dimensional-consistency-check]] + [[feedback-decisive-test-not-tautological]] applied to the test itself.

## 2. The structural skeleton (VERIFY each — REFUTE if wrong)
Claude's hand-bookkeeping; Codex must INDEPENDENTLY re-derive each dimension from its definition:
- Target `54·G·c_s⁵/(5·a⁵·c⁵)`: 3D `G = L³T⁻²M⁻¹`, `c_s = c = L·T⁻¹`, `a = L` ⇒ `M⁻¹·L⁻²·T⁻²` (pre-reg ~451; decision-13 §2).
- `m̂` (= `m̂0`) `= L⁻¹·T⁻¹·M⁻¹ᐟ²` (dimensional table; `pathA_21b` report ~129) ⇒ `m̂0² = M⁻¹·L⁻²·T⁻²` = the FULL target dimension.
- ⇒ `S_port·P0` (paper: `χ_Q·P0`) must be DIMENSIONLESS for homogeneity; `R_norm = 0 ⟺ ξ·(S_port·P0) = 54/5`,
  `ξ := m̂0²/[G·c_s⁵/(a⁵·c⁵)]` dimensionless.
- **`S_port` (code) ↔ `χ_Q` (paper):** the paper's law uses `χ_Q` (dimensionless OUTGOING-normalization scalar), `= 1` ONLY on
  the canonical compact outgoing branch (`eq:canonical-chiQ`); the code's `S_port` occupies that slot. RECONCILE: is
  `S_port ≡ χ_Q`, does it absorb `χ_Q`, or is it a separate convention? `χ_Q` is currently `OUTGOING_DTN_BRANCH_UNDERIVED` — if
  NOT canonically fixed/derived, it is a FREE dimensionless knob (a tunability channel). This reconciliation is a cheap,
  likely-decisive deliverable.

If Codex's independent derivation contradicts any of the above (e.g. `m̂0² ≠` target dim, or `S_port·P0` not dimensionless),
that is a `HOMOGENEITY_FAILURE` finding — REPORT it (name the missing power), do NOT patch the algebra to force homogeneity.

## 3. Required analysis (REQUIREMENTS + ACCEPTANCE CRITERIA; Codex designs the route — do not read this as a script spec)
SymPy with explicit `(M, L, T)` symbols (units RESTORED), as an ADDITIVE group in / sibling to
`src/stage1_solver/dimensional_check.py` (do NOT modify the pathA_18/19/20/20b/21 groups or any faithful operator). A narrow
Mathematica `.wl` cross-check of the dimensional algebra is required ([[feedback-dual-engine-required]]) — noting (fix #10) it
only verifies algebra under the SAME assumptions; it cannot independently prove pure-number closure.

**(A) Homogeneity check — must be able to FAIL, two distinct outputs.** Emit a `HOMOGENEITY_PASS/FAIL` separate from the
free-scale audit (a homogeneous comparison does NOT prove fail-able). Re-derive `[P0]`, `[N0]`, `[D0]`, `[K]`, `[B0]`, `[Z0]`,
`[S_port]`/`[χ_Q]`, `[m̂0]`, `[G]`, `[c_s]`, `[c]`, `[a]` FROM THEIR DEFINITIONS:
- `D0 = K − B0 − Z0`, `P0 = N0/D0`, with `N0 = P²/Δ²` (mixed-port) etc. (`patha_extraction.py`).
- **Do NOT assume `P0` is dimensionless from code names.** Explicitly resolve whether `N0`/`P0` is ALREADY the
  normalized stiffness-dimension coefficient or a RAW coefficient needing a hidden radiation-normalization factor
  (`Γ_port`, canonical `σ_Q`, or `a⁵/c_s⁵`). Trace the integrands + measures; watch for volume `a⁵` factors the pins hide.
- Then verify `[m̂0²·S_port·P0] == [54·G·c_s⁵/(5·a⁵·c⁵)]`.

**(B) Free-scale / dimensionless-knob LEDGER — the core; genuinely two-sided.** Represent every UNDERIVED form by an explicit
arbitrary dimensionless multiplier on its dimensional monomial (fix #5): `G = G_monomial·g_G(·)`, `m̂ = m̂_monomial·g_m(·)`,
`S_port = g_S(·)`, etc. — do NOT set any `g_*` to 1 unless DERIVED. Then substitute every in-hand model relation
(`c_s²=5Kρ⁴/m`; the dimensional monomials; `P0=N0/D0`) and:
- Express `ξ·(S_port·P0)` in fundamental params `{K, m (m_GNLS), ħ, a, ρ₀}` + the `g_*` factors.
- Enumerate and CLASSIFY each dimensionless knob/residual (fix #6 — broaden beyond `c/c_s`): at least `λγ = c_γ/c_s`,
  `χ_Q` (≈ `S_port`), `W_eff/a`, `R0/a`, `L/a`, `Θ_Q`, the `J`-selector / flux law, `α_J`, branch-kernel choices, and the
  `g_*` residuals. Class each as: **(a) fixed-by-prior-derivation**, **(b) branch-determined (target-blind)**, **(c) underived
  residual (form pending)**, **(d) TRUE FREE calibration knob**.
- (fix #7) Dependence on branch data (`R0`, `J`, profile) is NOT automatically "tunable" — it is `CONDITIONAL /
  PENDING_BRANCH_DERIVATION` UNLESS the directive/calibration permits CHOOSING it to hit `54/5`. Only class (d) = tunable.
- Catch any class-(d) free knob ALREADY explicit in the current skeleton (most likely `χ_Q`/`S_port` if not canonically fixed).

**(C) Order-of-magnitude (only if a class-(a/b) reduction already makes `ξ` numeric).** Report the rough order (≈1 vs ≈1e9).
Do NOT reverse-engineer `ξ` from `10.8/P0` (forbidden post-hoc value — the determination must be independent).

## 4. Outputs (structured — NOT a single hardcoded label; the report carries all three)
- **Homogeneity:** `HOMOGENEITY_PASS` (dims match) or `HOMOGENEITY_FAILURE` (name the missing power) — plus the resolved
  `P0`-dimensionless-or-not finding.
- **Knob ledger:** the table of every dimensionless knob/residual with its class (a/b/c/d) + provenance line ref.
- **Headline + next step (whichever the ledger supports — all must be REACHABLE):**
  - `HOMOGENEITY_FAILURE` → fix the bookkeeping first (cheap finding).
  - `TUNABILITY_CHANNEL_PRESENT` → a class-(d) free knob (e.g. `χ_Q`/`S_port` not canonically fixed) is already visible ⇒ the
    test is conditional/tunable UNLESS that knob is independently fixed; name it + bring to the user BEFORE deriving `G`.
  - `INDETERMINATE_NEEDS_FORMS` → homogeneous, no explicit class-(d) knob, but pure-number closure needs the forms ⇒ the
    decisive next step is deriving the minimal combination `ξ = m̂0²·S_port/[G·c_s⁵/(a⁵·c⁵)]`; list exactly which forms
    (`g_G`, `g_m`, `χ_Q` value) it hinges on. **(This is the honest EXPECTED outcome — do not force a cleaner one.)**

## 5. Anti-trap requirements (DERIVED-FORM GATE binds)
- (A) and (B) must be able to FAIL; a `HOMOGENEITY_FAILURE` / `TUNABILITY_CHANNEL_PRESENT` is a VALID, informative outcome — do
  NOT patch to force homogeneity or pure-number closure.
- **Negative control (fix #8):** include a control where an unresolved dimensionless factor multiplies `ξ·S_port·P0`; the
  machinery MUST emit `TUNABILITY_CHANNEL_PRESENT`/`INDETERMINATE`, NOT cancel it by a dimensional assumption. Also a planted
  dimensional mismatch MUST be caught by (A). These prove the checks can fail.
- NO `x == x` identity posing as a check; NO hardcoded verdict; NO symbol-assumption that declares an underived quantity
  dimensionless/unity to make it cancel (every underived form keeps its `g_*` multiplier).
- Units RESTORED throughout; flag every place a natural-unit pin hid a power of `a`/volume Jacobian/`ρ₀`.
- `timeout 600` APPLIES (SymPy/`.wl` DERIVATION script, NOT a solver run; the lifted cap is solver-only). Timeout (exit 124) ⇒
  reformulate, never raise ([[feedback-script-timeout-policy]]).
- Source paths: PDE source is `research/pde/paper/pde.tex` (fix #9 — NOT `research/pde_ledger/paper/pde.tex`).

## 6. Deliverables
- `reports/pathA_22a_dimensional_skeleton.md` — the dimension table (each ingredient + provenance), the homogeneity result, the
  `P0`-dimensionless finding, the classified knob ledger, and the headline + scoped next step (§4).
- Additive group in `dimensional_check.py` (or a new sibling module) implementing (A) + (B) + the negative controls; a narrow
  `tools/pathA_22a_dimensional_skeleton_crosscheck.wl`.
- Tests: a planted dimensional mismatch is CAUGHT by (A); a planted unresolved dimensionless factor yields
  `TUNABILITY_CHANNEL_PRESENT`/`INDETERMINATE` from (B).
- Run logs; do NOT commit (orchestrator commits after review, only when the user asks).

## 7. Review plan (after Codex exits 0)
Dimensional-fidelity audit ([[feedback-transliteration-fidelity-audit]]: code-vs-source dimensions term-by-term, incl. the
`S_port`↔`χ_Q` reconciliation) + adversarial review (is the knob ledger complete? is the free-scale audit genuinely two-sided,
or does it bias away from class (d)? do the negative controls actually fire?) — separate clean agents. Claude reads the
residual; brings the homogeneity result + knob ledger + scoped next step to the user. If `TUNABILITY_CHANNEL_PRESENT` or
`HOMOGENEITY_FAILURE`, escalate before any further derivation.
