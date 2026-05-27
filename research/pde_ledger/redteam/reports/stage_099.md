---
unit_id: 099
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T16:58:52Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage099_reduced_finish_line.md]
  paper_appendix: present
---

# Audit unit 099 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_099.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage099_reduced_finish_line.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read only the sections referencing stage 099: lines 20-90, 230-340, 1196-1235; section `\input{stages/stage_099}` at line 1232 and the part-IV audit-path map at line 25)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.txt`

## What the paper claims

Stage 099 ("Stage~116: The Reduced Finish Line After the Geometry-Lane Check", label `stage:099`) is the geometry-lane firewall ledger step that caps Part IV's first audit path (Stages 091–099, see `stage_appendix_part04.tex:25`). The card's bottom-line claim (quote block, `stage_099.tex:16`) is verbatim:

> The only remaining reduced theorem gate is `N_Q=1`, equivalently actual passive/outgoing quadrupole normalization.

The card identifies the conservative quadrupole module `Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` (line 13) and lists three explicit `\stagefield{Checks}` (`stage_099.tex:21-25`):
(a) static-limit check `epsilon_2 = epsilon_4 = 0` returns `c_pole = 1/4`;
(b) `l=0` and `l=2` orthogonality before applying the geometry firewall;
(c) any support/source success statement still carries the minimal-module hypothesis.

The card also lists Inputs (`stage_099.tex:9`): the Part III minimal isotropic module, the grouped real P_2 carrier, the static/dynamic geometry split, and the branch identity `K_0 K_4 = 4 K_2^2`. The notes (`moving_throat_pde_stage099_reduced_finish_line.md`) reinforce the gate: section 2 defines `N_Q := Kbar_0/Kbar_0^target`, `Kbar_0^target = 64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` after `Omega_Q = 3 c_s/(2a)`, and section 3 calls the open task equivalently "compute Kbar_0", "compute Gammabar_5", or "compute the scalar defect N_Q - 1." The appendix row at `stage_appendix_part04.tex:25` summarizes the audit path as "geometry-lane firewall and the conservative 3/4+1/4 module."

## What the script claims to verify

Both engines (sympy `*_audit.py` and Mathematica `*_audit.wl`) claim, per docstring/banner ("STAGE 82 — REDUCED FINISH LINE"), to verify that "the reduced finish line is a single normalization defect." The bottom-line assertions are:
- the geometric form of `K0_target` (i.e., `64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` after `Omega_Q = 3 c_s/(2a)`) — `*.py:36`, `*.wl:38`;
- `R_n - (N_Q - 1) = 0` for `n in {0, 2, 4, 5}` where `R_n := K_n/K_n_target - 1` (sympy `*.py:52-55`) and likewise in Mathematica `*.wl:63-66`.

The Mathematica script additionally derives `K_0, K_2, K_4` by `Series` expansion of `nQ * k0Target * yhatCons` in `omega`, then computes `Gamma5 = 9 Sqrt[K2^5/K0^3]`. The sympy script bypasses the series and writes `K_n = (N_Q * K_0_target) / (4 Omega_Q^{2n})` directly. There is no script-side test of `c_pole = 1/4`, `l=0 / l=2` orthogonality, the minimal-module hypothesis carry, the branch identity `K_0 K_4 = 4 K_2^2`, or the gate `N_Q = 1` itself.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Reduced gate `N_Q = 1` as final closure | (none — script just confirms `K0=N_Q*K0_target`-style defect propagation) | missing (paper does not require N_Q=1 to be derived here, only stated; see verdict) |
| Equivalence "compute Kbar_0 / Gammabar_5 / N_Q-1" (notes §3) | `R0=R2=R4=R5=N_Q-1` chain at `*.py:52-55`, `*.wl:63-66` | mismatch — assertions are tautological by construction (see F1) |
| `Kbar_0^target = 64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` after `Omega_Q = 3 c_s/(2a)` | `expect_zero("K0_target geometric form", ...)` at `*.py:36`, `*.wl:38` | match |
| Conservative module form `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` | `Yhat_Q_cons = sp.Rational(3,4) + sp.Rational(1,4)/...` at `*.py:31`, `*.wl:33` | match (printed, never asserted, but trivially the same expression) |
| Check (a) static limit `eps_2=eps_4=0 → c_pole = 1/4` | none | missing |
| Check (b) `l=0` and `l=2` orthogonality before geometry firewall | none | missing |
| Check (c) support/source success carries minimal-module hypothesis | none | missing |
| Input identity `K_0 K_4 = 4 K_2^2` | none (satisfied by construction inside the script; not asserted) | missing |

Set `paper_alignment: partial` — the script honors the algebraic core (K0_target geometric form) but does not exercise any of the three explicit Checks listed in the card, and the equivalence chain it does test is tautological by construction.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 36 | `expect_zero("K0_target geometric form", K0_target_geom - 54 G c_s^5/(5 a^5 c^5))` | `Kbar_0^target` two-form identity (notes §2) | yes — substitution identity, can fail if either form is wrong |
| A2 | sympy | 52 | `expect_zero("R0 - (N_Q - 1)", R0 - (N_Q - 1))` where `R0 = K0/K0_target - 1` and `K0 = N_Q * K0_target` | "equivalence" of K0 defect with N_Q-1 (notes §3) | no — `K0 = N_Q*K0_target` makes `R0 = N_Q - 1` algebraically inevitable |
| A3 | sympy | 53 | `expect_zero("R2 - (N_Q - 1)", R2 - (N_Q - 1))` where `K2 = K0/(4 Omega_Q^2) = N_Q*K0_target/(4 Omega_Q^2) = N_Q * K2_target` | "equivalence" of K2 defect with N_Q-1 | no — tautological by definition of K2 |
| A4 | sympy | 54 | `expect_zero("R4 - (N_Q - 1)", R4 - (N_Q - 1))` | "equivalence" of K4 defect with N_Q-1 | no — tautological by definition of K4 |
| A5 | sympy | 55 | `expect_zero("R5 - (N_Q - 1)", R5 - (N_Q - 1))` where `Gamma5 = 9 K2^(5/2)/K0^(3/2) = 9 (N_Q K2_target)^(5/2)/(N_Q K0_target)^(3/2) = N_Q * 9 K2_target^(5/2)/K0_target^(3/2)` | "equivalence" of Gamma_5 defect with N_Q-1 (notes §3) | no — `Gamma5/Gamma5_target = N_Q^(5/2)/N_Q^(3/2) = N_Q` is forced by algebra |
| A6 | mathematica | 38 | `expectZero["K0_target geometric form", k0TargetGeom - 54 G cS^5/(5 a^5 c^5)]` | same as A1 | yes |
| A7 | mathematica | 63 | `expectZero["R0 - (N_Q - 1)", r0 - (nQ - 1)]` — here `k0 = (kbar /. omega→0)` with `kbar = nQ*k0Target*yhatCons` and `yhatCons|_{omega=0} = 1` | same as A2 | no — `k0 = nQ*k0Target` by construction; tautological |
| A8 | mathematica | 64 | `expectZero["R2 - (N_Q - 1)", r2 - (nQ - 1)]` where `k2 = Coefficient[Series[kbar,{omega,0,4}], omega, 2]` | same as A3 | no — `Coefficient[nQ*k0Target*yhatCons, omega, 2] = nQ*k0Target/(4 omegaQ^2) = nQ*k2Target`; tautological |
| A9 | mathematica | 65 | `expectZero["R4 - (N_Q - 1)", r4 - (nQ - 1)]` | same as A4 | no — tautological |
| A10 | mathematica | 66 | `expectZero["R5 - (N_Q - 1)", r5 - (nQ - 1)]` | same as A5 | no — tautological |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py:38-55`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl:40-66`

**What's wrong:**
The four "equivalence" assertions `R_n - (N_Q - 1) == 0` (sympy lines 52-55; Mathematica lines 63-66) are algebraically guaranteed by the script's own definitions and cannot fail no matter what `Yhat_Q^cons` actually is.

In sympy: `K0 = sp.simplify(N_Q * K0_target)` (line 38), `K2 = sp.simplify(K0 / (4 * Omega_Q**2))` (line 39), `K4 = sp.simplify(K0 / (4 * Omega_Q**4))` (line 40), `Gamma5 = sp.simplify(9 * K2 ** sp.Rational(5, 2) / K0 ** sp.Rational(3, 2))` (line 41). Targets are constructed by the same recipes (lines 43-45). Then `R0 = K0/K0_target - 1 = N_Q - 1` identically; `R2, R4` likewise; and `Gamma5/Gamma5_target = N_Q^(5/2)/N_Q^(3/2) = N_Q` exactly, giving `R5 = N_Q - 1` exactly. So the assertion `R_n - (N_Q-1) == 0` is forced before any physics enters.

In Mathematica: `kbar = nQ*k0Target*yhatCons` (line 40), then `series = Normal[Series[kbar, {omega,0,4}]]` (line 41). Because `yhatCons|_{omega=0} = 1`, `Coefficient[yhatCons, omega, 0] = 1`, and the Series of `yhatCons` is `1 + omega^2/Omega_Q^2 + omega^4/Omega_Q^4 + O(omega^6)` while the targets use `K0_target/(4 Omega_Q^{2n})` — the n=0 coefficient of `nQ*k0Target*yhatCons` series is `nQ*k0Target * 1 = nQ*k0Target`, n=2 coefficient is `nQ*k0Target/Omega_Q^2 = 4*nQ*k2Target`. Wait — re-check: `yhatCons = 3/4 + 1/(4(1 - omega^2/Omega_Q^2))`. Expanding the pole: `1/(1-u) = 1 + u + u^2 + ...` with `u = omega^2/Omega_Q^2`, so `yhatCons = 3/4 + 1/4 + (1/4) omega^2/Omega_Q^2 + (1/4) omega^4/Omega_Q^4 + ... = 1 + omega^2/(4 Omega_Q^2) + omega^4/(4 Omega_Q^4) + ...`. Thus `Coefficient[kbar, omega, 0] = nQ*k0Target`, `Coefficient[kbar, omega, 2] = nQ*k0Target/(4 Omega_Q^2) = nQ*k2Target`, `Coefficient[kbar, omega, 4] = nQ*k0Target/(4 Omega_Q^4) = nQ*k4Target`. So in Mathematica too, the n=0,2,4 ratios reduce to `nQ` by construction; the assertions cannot fail.

The output transcripts (`*_sympy_audit.txt:14-22`, `*_mathematica_audit.txt:21-32`) print "0" / "PASS" for all four R-checks, exactly as forced.

**Why this matters:**
The script docstring (`*.py:3`) claims to verify "the reduced finish line is a single normalization defect." The intended non-trivial content is that the *four physically distinct radiative coefficients* (`Kbar_0, Kbar_2, Kbar_4, Gammabar_5`) all share the same normalization defect on the canonical branch — i.e., that the structural relations `K_2 = K_0/(4 Omega_Q^2)`, `K_4 = K_0/(4 Omega_Q^4)`, `Gammabar_5 = 9 Kbar_0/(32 Omega_Q^5)` (or equivalently the branch identity `K_0 K_4 = 4 K_2^2` listed in `\stagefield{Inputs}`) genuinely come out of `Yhat_Q^cons` and the canonical odd-coefficient formula. By writing `K_n` and `K_n_target` with the *same recipe* `R_target * (single relation in Omega_Q)`, both engines reduce the test to `(nQ*X)/X - 1 = nQ - 1` — the same scalar `nQ` cancels in both numerator and denominator. The script asserts nothing the algebra hasn't already forced.

**Required change:**
Rebuild the assertions so that the `K_n` and `Gammabar_5` are extracted from independently stated physical formulas rather than from `K_n = nQ * K_n_target` by construction. Two concrete restructurings (either is acceptable; Codex picks one and applies consistently across both engines):

- (a) Define the *actual-branch* moments by series-expanding `Kbar_0_actual * yhatCons` (with `Kbar_0_actual = nQ * Kbar_0_target` a free symbol relation only at the top), and define the *target* moments by series-expanding `Kbar_0_target * yhatCons`. Then the ratios `K_n/K_n_target = nQ` follow only because both series share the *same* `yhatCons` — the test still uses `yhatCons` non-trivially. Additionally, assert the branch identity `K_0 K_4 - 4 K_2^2 == 0` (which is genuinely a property of `yhatCons` and is listed as an Input in the card at `stage_099.tex:9`) and the odd-coefficient relation `Gammabar_5/(9 Kbar_0/(32 Omega_Q^5)) - 1 == 0` (a property of the assumed odd-pole structure, not of nQ).

- (b) Replace the four `R_n - (N_Q-1) == 0` assertions with a *substantive identity* the notes §3 actually claims: that, given the canonical relations `K_2 = K_0/(4 Omega_Q^2)`, `K_4 = K_0/(4 Omega_Q^4)`, `Gammabar_5 = 9 K_0/(32 Omega_Q^5)`, the ratio `Gammabar_5/(2G/(5c^5))` equals `(Kbar_0)/(Kbar_0^target)` (call it `N_Q`). Concretely, with `Kbar_0_target = 64 G Omega_Q^5/(45 c^5)`, assert `expect_zero("Gamma5_canonical - (2G/(5 c^5)) * (K0/K0_target)", 9*K0/(32*Omega_Q**5) - (2*G/(5*c**5)) * (K0/K0_target))` — this is a non-trivial check of the `Gammabar_5 = chi_Q N_Q (2G/(5 c^5))` factorization shown at `stage_appendix_part04.tex:286-288`. Here `K0` is a free symbol, not `nQ*K0_target`.

Codex should pick (b) — it is closer to the appendix's stated factorization and produces a clearly non-tautological check.

**Verification:**
After the fix, in the new outputs, at least one assertion should fail if the canonical odd-coefficient relation `Gammabar_5 = 9 K_0/(32 Omega_Q^5)` is changed (e.g., to `8 K_0/(32 Omega_Q^5)`), and the branch identity assertion `K_0 K_4 - 4 K_2^2 == 0` should fail if `K_4` is set to `K_0/(8 Omega_Q^4)`. The verifier confirms by checking the assertion expressions are *not* of the form `nQ*X - nQ*X` after symbolic expansion.

### F2 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_099.tex:21-25` (Checks list)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py` (no static-limit, no orthogonality, no minimal-module assertions anywhere)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl` (same)

**What's wrong:**
The stage card lists three explicit `\stagefield{Checks}` (`stage_099.tex:21-25`):

> \begin{verificationchecklist}
> \item Check the static limit `eps_2 = eps_4 = 0` returns `c_pole = 1/4`.
> \item Check `l=0` and `l=2` orthogonality before applying the geometry firewall.
> \item Check that any support/source success statement still carries the minimal-module hypothesis.
> \end{verificationchecklist}

None of these three are exercised in either engine. The script does print `Yhat_Q^cons(omega) = (Omega_Q^2 - 3*omega^2/4)/(Omega_Q^2 - omega^2)` (output line 13), from which a reader could read off `c_pole = 1/4` as the pole residue, but no assertion enforces it. The `l=0/l=2` orthogonality and the minimal-module carry are not addressable in the present scripts because no `l=0` Legendre weight, no `l=2` Legendre weight, and no support-side hypothesis appear in either file.

**Why this matters:**
The card's Checks block is meant to be the executable contract for the stage. A stage card that lists Checks the audit scripts do not exercise is exactly the failure mode the v2 paper-grounded audit is designed to catch: the paper claims more than the scripts verify. Either the scripts need to grow these checks, or the card needs to remove/relocate them to upstream stages (091-098) where the geometry-lane and orthogonality content actually lives.

This is `script_missing_paper_claim` and requires user resolution. The fix may legitimately be either direction:
- (a) Add the three checks to the scripts (e.g., assert `c_pole = 1/4` by extracting `Residue[yhatCons, omega → Omega_Q]` or by direct partial-fraction inspection; add explicit `l=0/l=2` orthogonality integrals; add an explicit upstream hypothesis dependency note); or
- (b) Trim the Checks block to one item — the geometric-form `Kbar_0^target` identity actually tested — and migrate the other two to the upstream stages 091-098 cards (which is where the geometry-lane firewall content originates per `stage_appendix_part04.tex:25`).

This is the user's call. Codex must not silently pick.

**Required change:**
See `## Resolve before fix_loop` block in the directive.

**Verification:**
After user resolution, either three new script assertions appear (with non-trivial residuals — e.g., `Residue` returns `Omega_Q^2/4` so `c_pole := -Residue/(Omega_Q^2) = -1/4` or whatever the agreed sign convention is) and exit-0 transcripts show all three pass; or the card's Checks block is trimmed to match the audit reality. Either outcome the verifier can confirm by re-grepping the script files and reading the new card.

### F3 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_099.tex:1` ("Stage~116") and `\label{stage:099}`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py:3` (`"""Stage 82 SymPy audit..."""`) and line 25 (`banner("STAGE 82 — REDUCED FINISH LINE")`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl:26` (`banner["STAGE 082 — REDUCED FINISH LINE"]`) and `*.wl:73` (`"STAGE 082 AUDIT PASSED"`)
- Output transcripts echo the stale numbers: `output/*_sympy_audit.txt:11,24`, `output/*_mathematica_audit.txt:11,34`.

**What's wrong:**
The stage's three identifiers disagree:
- Paper card title: "Stage~116" (`stage_099.tex:1`), label `stage:099`, filename `stage_099.tex`.
- Audit-unit filename: `..._stage099_...` (in scripts/, mathematica/, notes/).
- Script docstring/banner: "Stage 82" / "STAGE 082".

The script banners and the "STAGE 82 AUDIT PASSED" final line are stale carryovers from an earlier numbering scheme (consistent with the recent reorder commit `0d09ef6 fully reorder the pde ledger`). They contradict both the new paper card number ("Stage~116") and the audit-unit number (099).

**Why this matters:**
Confusion when reading transcripts in isolation. Not a math error; not blocking. But the audit-unit number `099` is what the verifier and orchestrator key on, and the paper card title `Stage~116` is what readers see; the script banner `STAGE 82` matches neither. A grep for the right ID across all artifacts must not surface a stale third number.

**Required change:**
Rename the script banners and final-line printouts to refer to "STAGE 099" (the audit-unit number), or to "STAGE 116" if the card title's number is the canonical user-facing label. The card title and the script banner should agree. Direction of resolution (use 099 or 116) is a user-facing convention question, so this is routed via `## Resolve before fix_loop` together with F2.

**Verification:**
The verifier re-runs `redteam exec-sympy 099` and confirms the banner and final-line in the regenerated `output/*_sympy_audit.txt` print the chosen ID, and likewise Mathematica.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py:31` (`Yhat_Q_cons = ...`; printed but never used or asserted against)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl:33` (`yhatCons = ...`; appears in `kbar` construction but the partial-fraction structure `3/4 + 1/4 * pole` is never asserted as such)

**What's wrong:**
The conservative module `Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` is the *load-bearing* expression for the entire reduced finish line per the notes §1.2 and the appendix `eq:app-part04-intro-3quarter-module` (referenced at `stage_appendix_part04.tex:70`). In sympy, the expression is bound to `Yhat_Q_cons` (line 31) and immediately overwritten by the unrelated direct construction `K0 = N_Q * K0_target` (line 38); the conservative-module form is never used after line 32. In Mathematica it does appear in `kbar = nQ*k0Target*yhatCons` (line 40), but no assertion checks that the static-slot value is 1 (i.e., `yhatCons|_{omega=0} = 1`), nor that the pole residue at `omega = Omega_Q` gives the `1/4 * Omega_Q^2 / (2 * Omega_Q) = Omega_Q/8` factor, nor that the static limit cleanly separates the `3/4` static part from the `1/4` pole part (the static-limit check from `stage_099.tex:22`).

**Why this matters:**
This is the place where the script's own docstring claim — "the reduced finish line is a single normalization defect" — is supposed to crystallize: the *form* `3/4 + (1/4)/(1-omega^2/Omega_Q^2)` is what forces the four moments to share a single defect. The script bypasses this and asserts only the trivial `R_n - (N_Q-1) = 0` post-construction. (This is also covered structurally by F1, but is listed separately because the cleanest fix is to add at least one assertion that *uses* `Yhat_Q^cons` non-trivially.)

**Required change:**
Add to both engines a non-tautological assertion that uses `Yhat_Q^cons` itself:
- Assert `Yhat_Q^cons|_{omega=0} - 1 == 0` (static slot equals 1).
- Assert that the residue of `Yhat_Q^cons` at `omega -> Omega_Q` equals `-Omega_Q/8` (or equivalent partial-fraction extraction confirming the pole coefficient is `1/4`). Concretely in sympy: `pole_residue = sp.residue(Yhat_Q_cons, omega, Omega_Q)` then `expect_zero("c_pole equals 1/4 (residue check)", pole_residue + Omega_Q/8)` (sign per the residue convention `1/4 * 1/(1 - omega^2/Omega_Q^2) = 1/4 * (-Omega_Q^2)/((omega-Omega_Q)(omega+Omega_Q))`, residue at `omega=Omega_Q` is `-Omega_Q/8`).
- Assert the branch identity `K_0 K_4 - 4 K_2^2 == 0` after re-deriving `K_n` from `Yhat_Q^cons` series, not from the multiplicative recipe.

These three assertions exercise the conservative module form non-trivially and address the static-limit check in `stage_099.tex:22`.

**Verification:**
After fix, the new output should print the residue and the static-slot value with explicit `PASS` for the partial-fraction structure and the branch identity. The verifier can perturb `Yhat_Q^cons` (e.g., to `1/2 + 1/2/(1 - omega^2/Omega_Q^2)`) and confirm the new assertions fail; under the canonical form they pass.

## Independent-derivation check (Mathematica)

The Mathematica script is *not* a literal transliteration. It computes `kbar = nQ*k0Target*yhatCons` and then uses `Series` and `Coefficient` to extract `K0, K2, K4`. The sympy script instead writes `K2 = K0/(4 Omega_Q**2)`, `K4 = K0/(4 Omega_Q**4)` by direct algebra. So the two engines take genuinely different routes to the same `K_n` values:

- sympy (`*.py:38-41`):
  ```
  K0 = N_Q * K0_target
  K2 = K0 / (4 * Omega_Q**2)
  K4 = K0 / (4 * Omega_Q**4)
  Gamma5 = 9 * K2**(5/2) / K0**(3/2)
  ```
- Mathematica (`*.wl:40-46`):
  ```
  kbar = nQ*k0Target*yhatCons
  series = Normal[Series[kbar, {omega, 0, 4}]]
  k0 = series /. omega -> 0
  k2 = Coefficient[series, omega, 2]
  k4 = Coefficient[series, omega, 4]
  gamma5 = 9*Sqrt[k2^5/k0^3]
  ```

The Mathematica path is actually closer to a first-principles extraction (read `K_n` off the Taylor expansion of the conservative module). The sympy path hardcodes the algebraic relations `K_2 = K_0/(4 Omega_Q^2)`. Both arrive at the same per-script result because `yhatCons` Taylor expands to `1 + omega^2/(4 Omega_Q^2) + omega^4/(4 Omega_Q^4) + ...`. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines arrive at identical R-values and identical printed-form for `K0_target geometric form = 0`:

- sympy output line 14: `K0_target geometric form = 0`; lines 15-18: all four `R_n - (N_Q - 1) = 0`; lines 19-22: `R_n = N_Q - 1` factored.
- Mathematica output line 14: `K0_target geometric form = 0`; lines 21-28: all four `R_n - (N_Q - 1) = 0` with `PASS` lines; lines 29-32: `R_n = -1 + nQ`.

Numerical/symbolic results agree. `engines_agree: true`.

The Mathematica output additionally prints intermediate `K_n` values from the series:
- `K0 = (64*G*nQ*omegaQ^5)/(45*c^5) = nQ * K0_target` (line 17)
- `K2 = (16*G*nQ*omegaQ^3)/(45*c^5) = nQ * K0_target/(4 Omega_Q^2)` ✓ (line 18)
- `K4 = (16*G*nQ*omegaQ)/(45*c^5) = nQ * K0_target/(4 Omega_Q^4)` ✓ (line 19)
- `Gamma5 = (2*G*nQ)/(5*c^5) = nQ * 9 * K0_target/(32 Omega_Q^5)` — check: `9 K0_target/(32 Omega_Q^5) = 9 * 64 G Omega_Q^5/(45 c^5 * 32 Omega_Q^5) = 576 G/(1440 c^5) = 2 G/(5 c^5)` ✓ (line 20)

All cross-checks consistent. Both engines also use compatible positivity assumptions on `G, c, c_s, a, Omega_Q (omegaQ), N_Q (nQ)` (sympy `*.py:27-29`; Mathematica `*.wl:29-31`). No `symbol_assumption_error`.

`outputs_fresh: true` — sympy script mtime Apr 3 < sympy output mtime May 11; Mathematica script mtime Apr 21 < Mathematica output mtime May 11.

## Verdict justification

`findings`. The scripts pass their printed assertions, but the four `R_n - (N_Q - 1) = 0` checks are tautological by construction in both engines (F1) — the test cannot distinguish the canonical `Yhat_Q^cons` from any alternative that shares its static-slot value, because the construction `K_n = N_Q * K_n_target` (or the equivalent Mathematica series extraction of `nQ * k0Target * yhatCons`) forces the ratios to `N_Q` independently of the underlying physics. The paper-side Checks (static-limit `c_pole = 1/4`, `l=0/l=2` orthogonality, support/source carry of the minimal-module hypothesis) are not exercised at all (F2), so the script does not honor the card's executable contract. The stale "Stage 82/082" banners disagree with the card's "Stage 116" and the audit-unit number 099 (F3). The conservative module `Yhat_Q^cons` is bound but barely used (F4). None of this is `UNFIXABLE` — the right fix is to rebuild the assertions around `Yhat_Q^cons` itself, add the three Checks, and resolve the banner-number convention. None of it is `CRITICAL_DOWNSTREAM` either: the paper's downstream use of `rho_alpha = 4/3`, `zeta_req = 1/3`, and `N_Q` depends on the equivalence-of-defects claim being correct, which it *is* (just not non-trivially verified by the current scripts). Attacks I tried that failed: I checked whether `Gamma_5 = 9*sqrt(K_2^5/K_0^3)` agrees with the appendix's `Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega_Q^5)` and confirmed they match identically given `K_2 = K_0/(4 Omega_Q^2)` (the script implicitly takes `chi_Q = 1`); I checked the Mathematica series expansion gives the correct `1/(4 Omega_Q^{2n})` factors; I checked the `K_0_target` two-form identity `64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` under `Omega_Q = 3 c_s/(2a)` — substitution gives `64 G (3 c_s/(2a))^5/(45 c^5) = 64 G * 243 c_s^5/(32 a^5)/(45 c^5) = 64*243/(32*45) * G c_s^5/(a^5 c^5) = 486/45 * G c_s^5/(a^5 c^5) = 54/5 * G c_s^5/(a^5 c^5)` ✓.

## Self-test notes

1. Variable independence — no `sp.diff` / `D[]` calls in either engine, so the missed-dependency trap does not apply.
2. Symmetry/parity — no integrals over unbounded domains, so the parity trap does not apply.
3. Trivial-case pre-check — the F1 finding *is* the trivial-case detection: substituting `N_Q = 1` gives `R_n = 0` in both engines, and substituting any other `N_Q` gives `R_n = N_Q - 1` exactly. The assertion `R_n - (N_Q-1) == 0` therefore fires identically for all `N_Q`, confirming it is `0 - 0 = 0` after the algebra collapses.
4. Path specifications — F2, F3, F4 prescribe edits to the existing script files. Both target paths are correct (`scripts/*.py` and `mathematica/*.wl`).
5. Paper round-trip — the proposed F1 / F4 assertions (`Yhat_Q^cons|_{omega=0} = 1`, residue at `omega=Omega_Q` equals `-Omega_Q/8`, branch identity `K_0 K_4 = 4 K_2^2`, `Gammabar_5 = 9 K_0/(32 Omega_Q^5)`) are all directly anchored in `stage_appendix_part04.tex:225-289` and `stage_099.tex:9,13,22`. No new paper-side claim is introduced. F2 routes to user resolution per the paper_misalignment policy.
