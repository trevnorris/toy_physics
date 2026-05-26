---
unit_id: 021
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md"]
  paper_appendix: present
---

# Audit unit 021 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_021.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (rows 9 and 64; stage `\input` on row 121)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt`

## What the paper claims

The card's Output paragraph states verbatim: "Stage~021 exports the reduced one-port self-energy \eqref{eq:app-stage021-self-energy}, the transfer factor \eqref{eq:app-stage021-transfer-factor}, and the wall-level outgoing quadrupole coefficient \eqref{eq:app-stage021-wall-odd}." Three distinct deliverables: (1) the exact conservative Σ_l^(EM+mix,cons)(ω) = [g_A² W_l(ω) + 2 g_A g_W R_l + g_W² A_l(ω)] / [A_l W_l − R_l²]; (2) the first-order outgoing transfer factor N_l(ω) = [A_l(ω) g_W,l + R_l g_A,l]² / [A_l(ω) W_l(ω) − R_l²]²; (3) the composed wall-level odd quadrupole coefficient δD_2^(odd)(ω) = −i N_2(0) (a^5/(27 c_s^5)) ω^5 + O(ω^7). Body equations additionally pin down gauge invariance of E_w, C_a; the z0/z2/z4 low-frequency coefficients; the manifestly nonnegative N_l(0); and the compact outgoing fingerprint Ŷ_2^(out)(ω) = 1 + a²ω²/(9c_s²) + 4a⁴ω⁴/(81c_s⁴) + i a⁵ω⁵/(27c_s⁵) + O(ω⁶). The notes additionally elaborate the scalar derivative-coupling check (notes §8, mirrored in part-appendix row 64's "scalar-compatibility criterion").

## What the script claims to verify

Both scripts (.py and .wl) verify, in five sections: (I) gauge invariance of the mixed Maxwell fields E_w and C_a under the scalar gauge transformation; (II) Euler–Lagrange equations for the reduced Lagrangian, the exact conservative self-energy Σ_cons obtained by eliminating A, W (sympy uses hand-written A_sol/W_sol with residual checks; mathematica uses LinearSolve and additionally asserts Σ_cons = g_A A_sol + g_W W_sol), and the closed forms for z0, z2, z4 via series-coefficient extraction matched against the notes formulas; (III) the first-order-in-Π_out correction giving the compact transfer factor N(ω) = (A_ker g_W + R g_A)² / Δ², and its zero-frequency limit; (IV) the compact outgoing l=2 fingerprint Ŷ_2^(out)(ω) via the spherical Hankel solution (sympy uses explicit j_2/y_2 series; mathematica uses the built-in SphericalHankelH1[2,·]) and the extracted Γ_5^port = a^5/(27 c_s^5); (V) the scalar derivative-coupling case (g_A = 0, g_W = ηω) showing the wall-level odd term starts at iω^3 rather than iω.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Gauge invariance of E_w, C_a (eq:app-stage021-mixed-fields) | sympy `mixed_field_gauge_invariance_audit`; .wl §I | match |
| Σ_l^(EM+mix,cons) (eq:app-stage021-self-energy) | sympy §II A/W residuals + closed form; .wl §II LinearSolve + `sigmaConsDerived - sigmaCons` | match |
| z0, z2, z4 coefficients (eq:app-stage021-z-coeffs) | sympy/`.wl` §II.3 series-coefficient comparisons | match |
| Transfer factor N_l(ω) (eq:app-stage021-transfer-factor) | sympy `N(omega) compact formula`; .wl `N(omega) compact formula` (derivative method) | match |
| N_l(0) closed form (eq:app-stage021-n0-positive) | sympy/`.wl` `N(0) positive-square form` | match |
| Ŷ_2^(out)(ω) fingerprint (eq:app-stage021-outgoing-fingerprint) | sympy/`.wl` §IV `Y2_hat minimal branch` and `Gamma5_port - a^5/(27 c_s^5)` | match |
| Composed δD_2^(odd)(ω) (eq:app-stage021-wall-odd) — Output item #3 | sympy/`.wl` §III prints `Dcorr = -i Gamma_port ω^5 N0`, but `Gamma_port` is a generic symbol; never substituted to `a^5/(27 c_s^5)` and never asserted as identity against the paper's specific eq | partial / missing |
| Scalar-compatibility criterion (notes §8; part-appendix row 64) | sympy/`.wl` §V `N_scalar leading term` and `scalar odd order` | match |

Set `paper_alignment: partial` — seven of eight deliverables have direct assertions; the third Output item (composed δD_2^(odd)) is computed in pieces (N_2(0) and Γ_5^port separately, in different sections) but never multiplied together and asserted as the literal paper formula.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 88 | `expect_zero("E_w gauge variation", ...)` | gauge invariance of E_w | yes |
| A2 | sympy | 89 | `expect_zero("C_a=F_aw gauge variation", ...)` | gauge invariance of C_a | yes |
| A3 | sympy | 128–130 | `expect_zero` on Q/A/W EL equations vs hand-written form | EL reduction (notes §3) | yes (notes-anchored) |
| A4 | sympy | 148–149 | `A/W exact solution residuals` (Schur-complement check) | self-energy derivation foundation | yes |
| A5 | sympy | 161–163 | `z0/z2/z4 formula` from toy rational | z-coeffs formulas (eq:app-stage021-z-coeffs) | yes |
| A6 | sympy | 175–177 | `Sigma z0/z2/z4` series match closed forms | z-coeffs specialized | yes |
| A7 | sympy | 222–225 | `N(omega) compact formula` | transfer factor N_l(ω) (eq:…transfer-factor) | yes |
| A8 | sympy | 228–231 | `N(0) positive-square form` | N_l(0) (eq:…n0-positive) | yes |
| A9 | sympy | 278–287 | `Y2_hat minimal branch` | Ŷ_2 fingerprint (eq:…outgoing-fingerprint) | yes |
| A10 | sympy | 290 | `Gamma5_port - a^5/(27 c_s^5)` | Γ_5^port quoted in notes §6 | yes |
| A11 | sympy | 323–326 | `N_scalar leading term` | scalar-compatibility (notes §8) | yes |
| A12 | sympy | 332–335 | `scalar odd order` (iω^3) | scalar-compatibility (notes §8) | yes |
| — | sympy | 239–242 | `Dcorr` printed only, no assert | δD_2^(odd) (eq:…wall-odd) Output item #3 | NO — print-only, generic Γ_port not specialized |
| B1 | .wl | 55–56 | `expectZero` on E_w/C_a gauge variation | same as A1/A2 | yes |
| B2 | .wl | 77–79 | `expectZero` on Q/A/W EL equations | same as A3 | yes |
| B3 | .wl | 91 | `sigmaCons from LinearSolve matches closed form` | Σ_cons consistency (independent derivation route vs. sympy) | yes |
| B4 | .wl | 95–96 | A/W exact solution residuals | foundation | yes |
| B5 | .wl | 104–106 | `z0/z2/z4 formula` | z-coeffs | yes |
| B6 | .wl | 117–119 | `Sigma z0/z2/z4` | z-coeffs specialized | yes |
| B7 | .wl | 141 | `N(omega) compact formula` (derivative method) | transfer factor | yes |
| B8 | .wl | 142 | `N(0) positive-square form` | N_l(0) | yes |
| B9 | .wl | 162–165 | `Y2_hat minimal branch` (built-in SphericalHankelH1) | Ŷ_2 fingerprint | yes |
| B10 | .wl | 166 | `Gamma5_port - a^5/(27 c_s^5)` | Γ_5^port | yes |
| B11 | .wl | 186 | `N_scalar leading term` | scalar-compatibility | yes |
| B12 | .wl | 187 | `scalar odd order` | scalar-compatibility | yes |
| — | .wl | 135, 140 | `dCorr` printed only, no assert | δD_2^(odd) Output item #3 | NO — print-only, generic gammaPort |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** script_missing_paper_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_021.tex:77-81` — Output paragraph
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_021.tex:71-75` — eq:app-stage021-wall-odd
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py:238-244` — Dcorr is computed and printed, no assertion
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:135, 140` — dCorr is computed and printed, no assertion

**What's wrong:**

The paper's Output paragraph (stage_021.tex:77–81) explicitly enumerates three exports: (1) the reduced one-port self-energy, (2) the transfer factor, and (3) "the wall-level outgoing quadrupole coefficient \eqref{eq:app-stage021-wall-odd}". Equation `eq:app-stage021-wall-odd` (lines 71–75) reads:

> δ D_2^{odd}(ω) = -i N_2(0) a^5/(27 c_s^5) ω^5 + O(ω^7)

with the explicit numeric coefficient `a^5/(27 c_s^5)` substituted in for Γ_5^port.

In the sympy script, Section III computes:

```
Gamma_port = sp.symbols("Gamma_port", positive=True, real=True)
Dcorr = sp.simplify(-I * Gamma_port * omega**5 * N0)
print("delta D_wall^(odd) =")
sp.pprint(Dcorr)
```

`Gamma_port` is a free symbolic constant; it is never specialized to `a^5/(27 c_s^5)` and `Dcorr` is never asserted against `-i N_2(0) * a^5/(27 c_s^5) * omega**5`. The matching Mathematica block (lines 122–142) has the same shape: `dCorr = -I gammaPort omega^5 n0` printed but not asserted. The .wl `gammaPort` is similarly never tied to `radius^5/(27 cS^5)`.

Sections III and IV separately verify the two pieces (N_l(0) closed form and Γ_5^port = a^5/(27 c_s^5)), and the FINAL ledger comments narrate "the first wall-level odd quadrupole coefficient is the outgoing-port coefficient a^5/(27 c_s^5) multiplied by the conservative mixed-sector transfer factor N(0)" — but the composition is asserted only in prose, not in code. A reader running the scripts would never see an assertion fail (or pass) for the paper's third Output deliverable.

Quoting paper card:

> Stage~021 exports the reduced one-port self-energy \eqref{eq:app-stage021-self-energy}, the transfer factor \eqref{eq:app-stage021-transfer-factor}, and the wall-level outgoing quadrupole coefficient \eqref{eq:app-stage021-wall-odd}.

Quoting sympy (line 239):

> `Dcorr = sp.simplify(-I * Gamma_port * omega**5 * N0)`

At the reduced one-port level the index l is just a label, so the natural script-side check is: substitute `Gamma_port → a^5/(27 c_s^5)` (no other index substitution required since N(0) is already the per-lane formula), then assert that the result equals `-i * (Ω_A² g_W + R g_A)² / (Ω_A² Ω_W² − R²)² * a^5/(27 c_s^5) * ω^5`.

**Why this matters:**

The paper Output is a three-item list, and one item has no script-side check. The composition is algebraically trivial (a substitution + multiplication), but absence of the assertion means a future edit that mis-substitutes one of the two pieces would silently pass the audit. More importantly, the paper Output line tells the reader the script proves three things; a faithful audit must prove all three.

**Required change:**

(Routed to user via the directive's `## Resolve before fix_loop` block — Codex must not auto-edit either paper or scripts to resolve this until the user picks a direction.)

The natural fix is to add one assertion in each script. For sympy, add to Section III (after `Dcorr` is defined) imports of the constants from Section IV, or restructure as a Section VI; the assertion shape is:

```python
expect_zero(
    "delta D_2^(odd) composed from Section III N(0) and Section IV Gamma5_port",
    Dcorr.subs(Gamma_port, a**5 / (27 * c_s**5))
    - (-I * (OA**2 * gW + R * gA)**2 / (OA**2 * OW**2 - R**2)**2 * a**5 / (27 * c_s**5) * omega**5),
)
```

For mathematica, mirror the same composed-identity assertion at the end of §IV or in a new §VI.

**Verification:**

After the fix, both `.py` and `.wl` should contain at least one `expect_zero` / `expectZero` assertion naming the composed wall-level odd coefficient. The new check should appear in the saved output as a `= 0` line. The paper's eq `eq:app-stage021-wall-odd` is then traceable to a specific script-side line.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. Three concrete points of structural divergence:

1. **Schur reduction of Σ_cons**: sympy hand-writes `A_sol`, `W_sol` and verifies them via residuals; mathematica uses `LinearSolve[matEAW, {gA, gW}]` and additionally asserts `sigmaConsDerived = gA aSol + gW wSol = sigmaCons` (line 91 — a check that sympy does not perform). Different derivation routes, with mathematica strictly stronger at this step.

2. **Transfer factor N(ω) derivation**: sympy expands `Sigma_full` in a Taylor series in `Pi` and divides out: `N_omega = (Sigma_first - Sigma_cons) / Pi` (line 214); mathematica takes the analytic derivative directly: `nOmega = D[sigmaFull, piOut] /. piOut -> 0` (line 133). Different mathematical operations producing the same formula — exactly the kind of cross-check the second-engine policy requires.

3. **Outgoing l=2 fingerprint**: sympy builds `j_2`, `y_2` from explicit closed forms (`j2a = (3/za^3 - 1/za) sin(za) - 3 cos(za)/za^2`, similarly y_2, then `h2a = j2a + I*y2a`); mathematica uses the built-in `SphericalHankelH1[2, za]`. Different sources for the same special function.

Conclusion: not a `mathematica_transliteration` finding.

## Engine cross-check

Both scripts pass all assertions (all `... = 0` lines in saved outputs). Final printed identities match across engines:

| Quantity | sympy output | mathematica output |
|---|---|---|
| Σ_cons | `(-2 R gA gW + gA² (-Ω_W²+ω²) + gW² (-Ω_A²+ω²)) / (R² - (Ω_A²-ω²)(Ω_W²-ω²))` | `(gW²(-oA²+omega²) + gA²(omega-oW)(omega+oW) - 2 gA gW r)/((oA-omega)(oA+omega)(omega-oW)(omega+oW) + r²)` |
| z0^(EM+mix) | `(Ω_A² g_W² + Ω_W² g_A² + 2 R g_A g_W)/(Ω_A² Ω_W² - R²)` | `(gW² oA² + gA² oW² + 2 gA gW r)/(oA² oW² - r²)` |
| N(0) | `(Ω_A⁴ g_W² + 2 Ω_A² R g_A g_W + R² g_A²)/(Ω_A⁴ Ω_W⁴ - 2 Ω_A² Ω_W² R² + R⁴)` | `(gW oA² + gA r)²/(-(oA² oW²) + r²)²` |
| Γ_5^port | `a^5/(27 c_s^5)` | `radius^5/(27 cS^5)` |
| Ŷ_2^(out)(ω) | `1 + a²ω²/(9 c_s²) + 4 a⁴ω⁴/(81 c_s⁴) + i a⁵ω⁵/(27 c_s⁵)` | `1 + omega² radius²/(9 cS²) + 4 omega⁴ radius⁴/(81 cS⁴) + I omega⁵ radius⁵/(27 cS⁵)` |

Engines agree (modulo symbol naming). `engines_agree: true`.

## Verdict justification

The two scripts pass all assertions they make, and those assertions cover seven of eight paper-side deliverables faithfully and non-tautologically. The cross-engine independence is real (different reduction routes, different ways to extract N(ω), different sources for the spherical Hankel function). The outputs are fresh.

The one gap: the paper Output paragraph names three exports, and the third (the composed wall-level odd quadrupole coefficient, eq:app-stage021-wall-odd) is computed only as a symbolic placeholder (`Dcorr = -i Gamma_port ω^5 N0`) and never substituted+asserted against the paper's specific value. This is `paper_misalignment / script_missing_paper_claim`: the paper says three things are exported; the scripts assert two and a half. The composition is algebraically trivial, but the auditor's job is to flag missing assertions even when the algebra is straightforward.

I also noted but did not file as findings: the cosmetic stale labels "STAGE 004 — …" and "FINAL STAGE-4 LEDGER" in the script banners and module docstring (sympy line 3, line 352; mathematica line 35, line 190) — these reflect the script's pre-renumbering history. They are documentation drift, not math findings; none of the ten finding categories cleanly apply. Worth fixing on the next pass but not blocking.

Attacks I tried that failed:
- Verified by hand that `dΣ/dΠ |_{Π=0} = (Aker gW + R gA)² / Δ²` from the Schur quotient. Matches.
- Checked sign convention: notes line 200 says `δD_l^(odd) = −i N_l(0) Γ_l^port ω^{2l+1}`. Script's Dcorr = `−i Γ_port ω^5 N0`. Sign consistent.
- Checked symbol assumptions: `OA, OW, M, K > 0`; `R, gA, gW` real (no positivity assumed). This matches the physical setup (frequencies positive, couplings can be either sign). No `simplify` under aggressive assumptions hides anything.
- Checked variable independence for `sp.diff` / `D[..., piOut]`: Sigma_full does depend on piOut/Pi; derivative is non-trivial.
- Checked parity / domains: no unbounded-domain integrals here; no parity traps.
- Checked the toy rational `(N0 - G2 ω²)/(D0 - S2 ω² + ω⁴)`: it is a non-tautological model whose ω-series coefficients are computed and then compared to the claimed closed forms. Not tautological.
- Checked Y2_hat construction: sympy uses explicit j_2/y_2 formulas, normalizes by Y2_static, then expands in k → ω/c_s. Mathematica uses built-in SphericalHankelH1. Both give the same minimal branch. Cross-engine check passes.

## Self-test notes

- **Variable independence**: For Section III's `sp.series(Sigma_full, Pi, 0, 2)` and `.wl`'s `D[sigmaFull, piOut]`, confirmed Sigma_full depends explicitly on Pi/piOut (sympy line 212 / .wl line 132). Derivatives are non-degenerate.
- **Parity / symmetry**: no unbounded-domain integrals; the only ω-series steps operate over a polynomial truncation, no parity-zero traps.
- **Trivial-case pre-check for proposed F1 fix**: substituting `Gamma_port = a^5/(27 c_s^5)` into `Dcorr = -I Gamma_port omega^5 N0` and subtracting `-I * N0 * a^5/(27 c_s^5) * omega^5` gives identically zero by inspection — assert_zero will pass for the correct composition. There is no risk of trivially-zero collapse since both N0 (a non-trivial rational function of OA, OW, R, gA, gW) and the constant a^5/(27 c_s^5) are present in the asserted side.
- **Paper round-trip**: the proposed F1 fix uses the same constant `a^5/(27 c_s^5)` that the paper card states verbatim and the script's Section IV already proves. No new `paper_misalignment` introduced.
- **Stale-output check**: sympy output mtime (May 21) > script mtime (Apr 1); .wl output mtime (May 21 15:16) > .wl mtime (May 21 15:14). Both fresh; banner "STAGE 004" is in the script source itself, so the output legitimately reflects it. No `stale_output` finding.
