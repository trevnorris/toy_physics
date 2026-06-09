---
unit_id: 197
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T18:51:52Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage197_conditional_packetA_closure_theorem.md]
  paper_appendix: present
---

# Audit unit 197 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_197.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage197_conditional_packetA_closure_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/paragraphs at lines 119–125, 1216–1326, 1467 referencing this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage197_conditional_packetA_closure_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage197_conditional_packetA_closure_theorem_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` states verbatim: "Proves \(\Delta_{\rm branch}=0\iff \chi_Q=1\) and hence \(\Delta_Q=\chi_Q-1\) is the sole Packet-A scalar." This is the original-stage audit card for *Conditional Packet-A closure* (anchor MTDC-T9.5). The notes (the authoritative, fuller derivation) enumerate five deliverables: (1) the exact collapse of the full Stage-191 Packet-A residual `(a₂,b₂,a₄,b₄,a_{P0},b_{P0},Δ_pole,Δ_norm)` to a single scalar slot once the carried isotropic one-pole front end is imposed, i.e. `Δ_branch = (0,…,0,Δ_norm)`; (2) the exact equivalence `Δ_branch=0 ⟺ χ_Q=1` on the natural point-particle source-map branch; (3) the equivalent forms `χ_Q=1 ⟺ N_Q=1 ⟺ Δ_Q:=χ_Q−1=0`; (4) the deformation-algebra form `χ_Q=1 ⟺ 3S(β⁵−1)+Σ₀+27Σ₅=0`, with `χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` carried from Stage 194; (5) the linearized finish-line map `χ_Q−1 = 5ε_β + δΣ₀/(3S) + 9δΣ₅/S + O(2)` plus the higher-odd (`O(ω⁷)`) irrelevance of the theorem. The carried constant `P₀^target = 54Gc_s⁵/(5a⁵c⁵) > 0` (notes line 152/184) is what makes the iff bracket-sign argument go through.

## What the script claims to verify

Both scripts march through the same seven labeled sections (I–VII) and the docstring/banner "EXACT CONDITIONAL PACKET-A CLOSURE THEOREM". Section I proves the isotropic outgoing lane kills both grouped-anisotropy coordinates `a_{P0}=b_{P0}=0`; II asserts the carried residual is `(0,…,0,Δ_norm)`; III extracts `χ_Q` from a `z⁵` coefficient of a normalized DtN response, confirms it equals the Stage-194 deformation-algebra formula `3(Sβ⁵+9Σ₅)/(3S−Σ₀)`, applies the source map `N_Q=1/χ_Q`, and forms `Δ_norm = P₀^target(1/χ_Q−1)`; IV checks the iff via algebraic identities (SymPy) and `Reduce`/`Equivalent` (Mathematica), plus a negative control at `χ_Q=6/5` giving a nonzero defect; V verifies the deformation-numerator form `(3S−Σ₀)(χ_Q−1)=3S(β⁵−1)+Σ₀+27Σ₅`; VI checks the linearized slope `5ε_β+δΣ₀/(3S)+9δΣ₅/S`; VII checks `∂χ_Q/∂L7 = ∂Δ_norm/∂L7 = 0` (higher-odd irrelevance). Every section maps onto one of the five paper deliverables.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) Residual collapse to `(0,…,0,Δ_norm)` | py II (L64–67); wl II (L91–94) | match |
| isotropic lane `a_{P0}=b_{P0}=0` (notes §1.2) | py I (L50–56); wl I (L76–86) | match |
| (2) `Δ_branch=0 ⟺ χ_Q=1` | py IV (L124–142); wl IV Reduce/Equivalent (L177–218) | match |
| (3) `χ_Q=1 ⟺ N_Q=1 ⟺ Δ_Q=0` | py IV (L133–138); wl IV (L199–214) | match |
| (4) deformation gate `3S(β⁵−1)+Σ₀+27Σ₅=0`, `χ_Q=3(Sβ⁵+9Σ₅)/(3S−Σ₀)` | py V (L149–166) + III extractor (L94–95); wl V (L223–240) + III extractor (L138–150) | match |
| (5) linearized map + higher-odd irrelevance | py VI–VII (L176–210); wl VI–VII (L252–271) | match |
| `P₀^target=54Gc_s⁵/(5a⁵c⁵)` | py literal (L76); wl literal + `Γ₅`-native re-derivation (L99–103) | match |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological script-side check on both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55–56 | `expect_zero(a_{P0}), expect_zero(b_{P0})` from grouped projector of `(P0,P0,P0)` | isotropic lane (notes §1.2) | yes |
| A2 | sympy | 67 | first-seven-slots `== 0` | deliverable 1 | partial (asserts the constructed zeros; collapse itself is upstream) |
| A3 | sympy | 108 | `chi_from_series − chi_from_def == 0` | deliverable 4 (χ_Q formula) | yes (series extractor vs closed form) |
| A4 | sympy | 109–112 | `Δ_norm − P₀^target(1/χ_Q−1) == 0` | deliverable 2/3 | yes |
| A5 | sympy | 124–131 | identities encoding `Δ_norm=0 ⟺ χ_Q=1` | deliverable 2 | yes |
| A6 | sympy | 139–142 | negative control: `Δ_norm(χ_Q=6/5) ≠ 0` | deliverable 2 (sharpness) | yes |
| A7 | sympy | 155–166 | deformation-numerator identity | deliverable 4 | yes |
| A8 | sympy | 185–202 | linearized χ_Q and Δ_norm slopes | deliverable 5 | yes |
| A9 | sympy | 209–210 | `∂χ_Q/∂L7 = ∂Δ_norm/∂L7 = 0` | deliverable 5 (higher-odd) | yes |
| B1 | mathematica | 85–86 | `expectZero` on native trace-projector `a_{P0},b_{P0}` | isotropic lane | yes |
| B2 | mathematica | 93–94 | first-seven slots + survivor slot `==0` | deliverable 1 | partial (same as A2) |
| B3 | mathematica | 103 | `p0FromGamma − p0Target == 0` (native `Γ₅` derivation) | `P₀^target` value | yes (independent re-derivation) |
| B4 | mathematica | 134–135 | native-Solve `Σ₂,Σ₄` agree with SymPy closed forms | deliverable 4 (matching) | yes (cross-engine) |
| B5 | mathematica | 149–150 | Series-coeff vs 5th-derivative vs SymPy `χ_Q` | deliverable 4 | yes (triple cross-check) |
| B6 | mathematica | 164–169 | source-map solve, `N_Q=1/χ_Q`, `Δ_norm` | deliverable 2/3 | yes |
| B7 | mathematica | 199–206 | `Equivalent[Reduce[Δ_norm==0],χ==1]`, `Equivalent[Reduce[N_Q==1],χ==1]` | deliverable 2/3 | yes (genuine zero-set logic) |
| B8 | mathematica | 207–218 | `Δ_Q` forms + negative control `χ=6/5` | deliverable 2/3 | yes |
| B9 | mathematica | 229–240 | deformation-numerator identities | deliverable 4 | yes |
| B10 | mathematica | 258–271 | linearized slopes + `∂/∂tail7=0` | deliverable 5 | yes |

A2/B2 are flagged "partial" because the residual-collapse vector is *constructed* with seven literal zeros and the assertion only re-confirms those zeros; the physical collapse argument (`a₂=…=Δ_pole=0`) is carried from Stages 191/193 and not re-derived here. This is legitimate carry-forward (the stage's own scope is the *normalization slot* finish line), not a finding — sections III–VII carry the real, non-tautological content. The isotropic-lane kill (A1/B1) is the one piece of the collapse this stage *does* derive in-script, and both engines do it non-trivially (the grouped trace/anisotropy projector applied to `(P0,P0,P0)`).

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.txt:3` and `:151`

**What's wrong:**
The committed SymPy output carries a stale pre-renumber banner `STAGE 180 — EXACT CONDITIONAL PACKET-A CLOSURE THEOREM` (line 3) and `STAGE 180 LEDGER` (line 151), whereas the current `.py` (line 42) prints `STAGE 197 …`. `180 = 197 − 17`, the known +17 renumber offset. Independently, the SymPy `.txt` mtime is `2026-06-01 12:03:34`, older than the `.py` mtime `2026-06-03 15:59:11`, so the captured transcript predates the current script. (The Mathematica `.txt` mtime `2026-06-01 12:03:34` is newer than its `.wl` `2026-06-01 12:02:15` and its banner already reads `STAGE 197`, so the Mathematica side is fresh.) The numerical/symbolic content of the stale SymPy transcript still matches what the current script would produce (all sections agree with the fresh Mathematica run), so the staleness is cosmetic/banner-level only.

**Why this matters:**
A reader trusting the committed transcript sees the wrong stage number. No math is affected.

**Required change:**
None for Codex (no source edit). The orchestrator's independent re-run refreshes `scripts/output/…_sympy_audit.txt`; after re-run, line 3 should read `STAGE 197 — EXACT CONDITIONAL PACKET-A CLOSURE THEOREM` and line ~151 `STAGE 197 LEDGER`.

**Verification:**
After re-run, `grep "STAGE 180" scripts/output/moving_throat_pde_stage197_conditional_packetA_closure_theorem_sympy_audit.txt` returns nothing.

## Independent-derivation check (Mathematica)

**GENUINELY INDEPENDENT — not a transliteration.** The two engines extract the same `χ_Q` but from structurally distinct routes:

1. **DtN fingerprint source.** The `.wl` derives the outgoing operator from the *actual physics object*:
   `lambdaOut = FunctionExpand[z*D[SphericalHankelH1[2, z], z]/SphericalHankelH1[2, z]]` (wl L105), i.e. `z d/dz ln h₂^(1)(z)` from the native spherical Hankel function, then `Series` to `z⁵` (wl L106). The `.py` never touches a Hankel function — it *hardcodes* the already-expanded `Λ`-coefficients `L0=-3S+Σ0, L2=S β²/3+Σ2, L4=S β⁴/9+Σ4, L5=S β⁵/9+Σ5` (py L81–84). So the `.wl` regenerates from first principles exactly the coefficients the `.py` carries as literals.

2. **Canonical-even matching.** The `.wl` *Solves* `{coeff(z²)==1/9, coeff(z⁴)==4/81}` for `{sigma2,sigma4}` (wl L117–124) and asserts uniqueness (`Length===1`, wl L125). The `.py` *hardcodes* the closed-form solutions `Sigma2_match = -(3Sβ²-3S+Σ0)/9`, `Sigma4_match = -(3Sβ⁴-3S+Σ0)/27` (py L90–91). The `.wl` then cross-checks its native solve against the SymPy literals (`expectZero["Sigma_2 agreement with SymPy result", …]`, wl L134–135) — a cross-engine confirmation, not an echo.

3. **`χ_Q` extraction + iff logic.** The `.wl` extracts `χ_Q` two independent ways — `Coefficient[matchedSeries,z,5]/(I/27)` and the 5th-derivative route `(D[…,{z,5}]/.z->0)/5!/(I/27)` — and asserts they agree (wl L138–149). It proves the iff with genuine logic: `Equivalent[Reduce[Δ_norm==0, χ, Reals], χ==1]` and the same for `N_Q==1` (wl L177–206). The `.py` instead encodes the iff as hand-built *algebraic identities* (`χ_Q*(Δ_norm/P₀+1)−1==0`, py L124–127) — a different proof technique for the same fact. The `.wl` also adds a value the `.py` lacks entirely: re-deriving `P₀^target` from the `Γ₅=2G/(5c⁵)` normalization (`27 c_s⁵ Γ₅/a⁵`, wl L99–103) versus the `.py` literal (py L76).

Evidence quoted (2–3 corresponding sections): py L81–94 (hardcoded `Λ` coeffs + hardcoded `Σ`-match) vs wl L105–135 (native Hankel + `Solve`); py L124–131 (algebraic iff identities) vs wl L177–206 (`Reduce`/`Equivalent` iff). Different choreography, different primitives, mutual cross-checks. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical bottom-line forms and all checks PASS:
- `χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` — sympy `.txt` L59–63; mathematica `.txt` L35.
- `Δ_norm(natural) = 18Gc_s⁵(−3Sβ⁵+3S−Σ₀−27Σ₅)/(5a⁵c⁵(Sβ⁵+9Σ₅))` — sympy `.txt` L71–76; mathematica `.txt` L42/L53 (algebraically identical, written in the `(−1+(3S−σ0)/(3β⁵S+27σ5))` factored form).
- Negative control `Δ_norm(χ_Q=6/5) = −9Gc_s⁵/(5a⁵c⁵)` — sympy `.txt` L106; mathematica `.txt` L70 (identical).
- Linearized χ_Q and Δ_norm — sympy `.txt` L131–142; mathematica `.txt` L89–94 (identical).
- Higher-odd `∂/∂L7 = 0` — both `.txt` (sympy L147–148, mathematica L99–102).

`engines_agree: true`.

## Verdict justification

`findings`, single informational `stale_output` only. Attacks tried and failed: (a) **transliteration probe** — the `.wl`'s native `SphericalHankelH1` derivation and `Solve`-based even-matching defeat any port allegation; the engines reach the same `χ_Q` from genuinely different primitives and cross-check each other. (b) **tautology probe on the iff** — the `.py` uses hand-built algebraic identities and the `.wl` uses `Reduce`+`Equivalent` over `Reals`; both are backed by a negative control at `χ_Q=6/5` that returns a nonzero defect (so the checks can fail). (c) **derivative-zero trap (VII)** — `∂χ_Q/∂L7` and `∂Δ_norm/∂L7`: `χ_Q` genuinely does not depend on `L7`/`tail7` (the `z⁷` tail does not reach the `z⁵` coefficient), so these `=0` results are physically meaningful (higher-odd irrelevance), not vacuous derivatives of a constant w.r.t. an absent variable — `tail7`/`L7` IS present in the constructed `deformedDtn`/`Y_stage194_hi` and simply drops out at theorem order. (d) **value reconciliation** — every emitted deliverable reconciles to the notes/appendix (see below). Read the card, notes, and appendix in full; the script's verified claim matches the paper's `\stagefield{Output}` exactly. One paper-side card-text lag is noted (not a script finding): card line 11 still says "Mathematica audit: none yet" despite the present, passing `.wl` — this is a paper/prose file the red-team does not edit; flag for the orchestrator's paper-card sync, not Codex.

## Self-test notes

Checked: (1) Variable-independence — VII's `D[chi, tail7]`/`D[Δ_norm, tail7]`: `tail7`/`L7` is genuinely a coefficient in the constructed series (py L86, wl L112), so the zero result reflects real higher-odd irrelevance, not a derivative w.r.t. a symbol absent from the expression. (2) Trivial-case — negative control `χ_Q=6/5` gives nonzero `−9Gc_s⁵/(5a⁵c⁵)` on both engines, confirming the iff checks are falsifiable. (3) Symbol domains — `χ_free != 0`, `3S−σ0 != 0`, `Sβ⁵+9σ5 != 0`, `mHat0>0`, `1+δQ != 0` (wl L68–71) are exactly the non-degeneracy conditions the algebra needs (denominators of `χ_Q`, `N_Q`); positivity of `P₀^target` constituents (`radius,soundSpeed,lightSpeed,bigG>0`) matches the notes' physical-branch sign argument. No numbering self-label in the card requiring a fix (`+17 of 197 = 214` absent from the `\stagefield{Purpose}`). Only fix is the single informational `stale_output`.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` | py L95; wl L146; sympy.txt L59–63 / math.txt L35 | notes L232 (`\frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0}`); appendix L1297–1299 | MATCH |
| `χ_Q−1 = (3S(β⁵−1)+Σ₀+27Σ₅)/(3S−Σ₀)` | py L149,157; wl L223,231; sympy.txt L116–119 / math.txt L77,79 | notes L237; appendix L1301–1302 | MATCH |
| closure gate `3S(β⁵−1)+Σ₀+27Σ₅=0` | py L149; wl L223; math.txt L77 | notes L62,244,362; card `\stagefield{Output}` (anchor); appendix L1467 | MATCH |
| `Δ_norm^pt = P₀^target(1/χ_Q−1) = 18Gc_s⁵(−3Sβ⁵+3S−Σ₀−27Σ₅)/(5a⁵c⁵(Sβ⁵+9Σ₅))` | py L98,160; wl L158,224; sympy.txt L71–76 / math.txt L42 | notes L160–164 (`Δ_norm = P₀^target(χ_Q⁻¹−1)`), L249–255 (full deformation form) | MATCH |
| `P₀^target = 54Gc_s⁵/(5a⁵c⁵)` | py L76; wl L101; sympy.txt L53–57 | notes L152, L184 | MATCH |
| linearized `χ_Q−1 = 5ε_β + δΣ₀/(3S) + 9δΣ₅/S` | py L182; wl L250; sympy.txt L131–135 / math.txt L89 | notes L274–278 | MATCH |
| equivalence `χ_Q=1 ⟺ N_Q=1 ⟺ Δ_Q:=χ_Q−1=0` | py L133–138; wl L199–214; math.txt L62–69 | notes L51–57, L210–219; card `\stagefield{Output}`; appendix L1323–1325 | MATCH |

INTERNAL scaffolding (accounted for, no finding): `gamma5Target = 2G/(5c⁵)` (wl L99, an intermediate used solely to re-derive `P₀^target`); `Σ₂_match = -(3Sβ²-3S+Σ₀)/9`, `Σ₄_match = -(3Sβ⁴-3S+Σ₀)/27` (canonical-even matching intermediates, py L90–91 / wl L128–135 — derived/cross-checked, not stage deliverables); `Δ_norm(χ_Q=6/5) = -9Gc_s⁵/(5a⁵c⁵)` (negative-control probe value); the constructed residual vectors and pass/fail flags.

All six stated deliverable values plus the equivalence reconcile to the notes and/or card and appendix. No `value_mismatch` and no `script_missing_paper_claim` arise from the reconciliation; the only finding remains the informational `stale_output` (F1).
