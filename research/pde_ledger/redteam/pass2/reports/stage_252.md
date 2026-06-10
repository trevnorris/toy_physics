---
unit_id: 252
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 252 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_252.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows: line 102 status row; lines 320-332 main theorem heat-fraction + safe-energy; line 369-371 stage-pair note + `\input`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.txt`

## What the paper claims

Stage 252 is a channel-resolved energy compiler built on the Stage-251 microscopic export kernel `K_exp(s)=Gamma_3 s^3 + Gamma_5 s^5`. It splits each coefficient into vacuum/lattice pieces and proves: (1) the exact exported-energy ledger `E_vac = Gamma3vac I2 + Gamma5vac I3`, `E_lat = Gamma3lat I2 + Gamma5lat I3` and total; (2) the partition fractions reduce to one event-shape quotient `r_V = I3/I2`, `f_lat = (Gamma3lat + Gamma5lat r_V)/(Gamma3 + Gamma5 r_V)`, `f_vac+f_lat=1`; (3) the speed-drift law `df_lat/dr_V = (Gamma5lat Gamma3vac - Gamma3lat Gamma5vac)/(Gamma3+Gamma5 r_V)^2` with endpoint limits `Gamma3lat/Gamma3` and `Gamma5lat/Gamma5`; (4) the single-growth specialization `V=Vin e^{st}` giving `I2 = Vin^2 s^3 (e^{2sT}-1)/2`, `I3 = s^2 I2`, `r_V = s^2`, and event-equivalent rates `gamma^eq = Gamma3 s^2 + Gamma5 s^4`; (5) the safe-edge theorem `E_exp,min^safe = (Vin^2/2)(e^2-1) mu_eta (s0^2-sc^2)` plus `gamma_eff,safe^eq = mu_eta (s0^2-sc^2)/sc`; (6) the Session-IV 3:1 split as the surface `Gamma3lat+Gamma5lat sc^2 = 3(Gamma3vac+Gamma5vac sc^2)`, speed-independent on the phi-family (`phi=3/4`); (7) a Session-IV benchmark calibration (`gamma_eff,safe^eq ≈ 87.26925235`, `gamma_vac/lat ≈ 21.817/65.452`, `Vin_match ≈ 8.21771260e-3`, `E_vac/E_lat ≈ 0.00258365/0.00775095`). The card's `\stagefield{Verification}` line states "SymPy audit: ... Mathematica audit: none yet."

## What the script claims to verify

The SymPy script (docstring items 1-7) and the Mathematica script (M1-M9) each independently assert every one of the seven deliverables above: the ledger and fraction forms (M1), drift law (M2), endpoint limits (M3), exponential shape integrals and `r_V=s^2` (M4), event-equivalent rates (M5), the safe-edge energy theorem with the Stage-251 safe-rule substitution (M6), the safe-edge rate identity (M7), the 3:1 surface and phi-family (M8), and the Session-IV numeric benchmark (M9). All `assert`/`expectZero`/`expectApprox` pass in both saved outputs.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Ledger E_vac/E_lat/E_tot | sympy L45-53, M1 L76-94 | match |
| Fractions f_vac/f_lat, sum=1 | sympy L58-60, M1 L92-94 | match |
| Drift law + endpoints | sympy L66-77, M2/M3 L98-115 | match |
| `I2 = Vin^2 s^3(e^{2sT}-1)/2`, `I3=s^2 I2`, `r_V=s^2` | sympy L104-105, M4 L137-139 | match |
| Event-equiv rates | sympy L106-107, M5 L151-152 | match |
| Safe-edge energy theorem | sympy L137-139, M6 L168-172 | match |
| Safe-edge rate identity | sympy L140,148, M7 L182-185 | match |
| 3:1 surface + phi-family | sympy L178-190, M8 L226-231 | match |
| Session-IV benchmark numbers | sympy L223-232, M9 L263-271 | match |
| Card `\stagefield{Verification}`: "Mathematica audit: none yet" | a working `.wl` (M1-M9, all PASS) now exists | mismatch (stale card text) |

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 58-60 | `simplify(fv/fl - closed form)==0`, sum=1 | claim 2 | yes |
| A2 | sympy | 75-77 | drift law + two limits ==0 | claim 3 | yes |
| A3 | sympy | 104-107 | I2, I3=s^2 I2, E_vac/lat = gamma^eq I1 | claim 4/5 | yes |
| A4 | sympy | 137-139 | I1_safe, Et_safe_expected==Et_safe, reduced form | claim 5 | yes |
| A5 | sympy | 140 | `safe_combo/sc - gamma_eq(sc) ==0` (genuine: sc^3/sc=sc^2) | claim 5 (rate) | yes |
| A6 | sympy | 148 | bridge: safe_combo→mu_eta(...) then /sc == mu_eta(...)/sc | claim 5 (rate) | partial (connective) |
| A7 | sympy | 178,190 | 3:1 numerator==surface; phi-family fractions | claim 6 | yes |
| A8 | sympy | 223-232 | benchmark numerics vs targets | claim 7 | yes (calibration pin) |
| M1-M9 | mathematica | 92-271 | engine-independent mirror of A1-A8 | claims 1-7 | yes (M7 raw quotient = A5 = load-bearing) |

## Findings

### F1 — paper_misalignment (subtype: paper_missing_script_claim)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_252.tex:4`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage252_..._mathematica_audit.wl` (entire file, 9 sections M1-M9, all PASS)

**What's wrong:**
The stage card's verification line reads: `\stagefield{Verification}{SymPy audit: ... Mathematica audit: none yet.}` A complete, passing Mathematica audit (`..._mathematica_audit.wl`, sections M1-M9) now exists and was created in pass-1. The card text is stale (card-text-lag), not a math defect.

**Why this matters:**
Reader-facing card understates verification coverage; the dual-engine status is now satisfied but the card still says "none yet."

**Required change:**
This is a paper-side text fix routed to the user — Codex applies nothing. See `## Resolve before fix_loop` in the directive.

**Verification:**
Card line 4 should read "Mathematica audit: \StageFile{mathematica/...stage252...mathematica_audit.wl}." after user approval.

## Independent-derivation check (Mathematica)

The `.wl` is a genuine second engine for the load-bearing symbolic work, with porting-like scaffolding. Each engine independently performs the hard symbolic operations: Mathematica `Integrate[D[Vevent,{t,k}]^2,{t,0,T}]` (L120-131) vs SymPy `sp.integrate(Vdd**2,(t,0,T))` (L89-91); Mathematica `D[flat,rV]` (L98) vs SymPy `sp.diff(fl,rV)` (L66); Mathematica `Limit[...]` (L107-109) vs SymPy `sp.limit(...)` (L68-69). These are re-derived per engine, not echoed. The scaffolding (variable names Evac/Elat/Eexp, `shapeRule`/`safeRule` substitution technique, M8 phi-family substitution, and the M9 benchmark block with identical hardcoded targets/tolerances) closely parallels the `.py` outline. Judgment: independent enough — the verifiable identities are recomputed by each engine; the parallel outline reflects the stage's fixed identity sequence rather than a blind transliteration. NOT a `mathematica_transliteration` finding, but the M9 benchmark targets are shared literals in both engines (calibration pins, not re-derived).

## Engine cross-check

Both outputs agree exactly. f_vac/f_lat forms identical (sympy L8-9 ↔ wl L12-13); drift law identical (sympy L15 ↔ wl L24); I1/I2/I3 identical (sympy L22-24 ↔ wl L41-43); safe-energy reduced identical (sympy L36 ↔ wl L66, both `(Vin^2/2)(e^2-1)mu_eta(s0^2-sc^2)`); 3:1 surface identical (sympy L44 ↔ wl L85-86); benchmark numerics identical to ~1e-14 (M9 PASS lines 108-125). No `engine_disagreement`.

## Verdict justification

`findings` — exactly one finding, low-severity `paper_misalignment` (stale "Mathematica audit: none yet" card line), routed to the user. The math is sound on both engines. Attacks tried that failed: (a) **Pass-1 round-trip re-check** — the old `gamma_safe_eq` X−X tautology is resolved: the load-bearing identity is now sympy L140 / wl M7-L182 (`safe_combo/sc == Gamma3 sc^2 + Gamma5 sc^4`), a genuine algebraic identity that fails if the kernel powers are wrong; the bridge (L148/M7-L185) only documents the safe-edge reduction and is connective. (b) **subs-fired check** — confirmed both `safe_combo` substitutions actually fire: sympy output L36 shows `E_safe reduced` reduced to the mu_eta form (so L120 fired), and the `(G3 sc^3+G5 sc^5)/sc` Add survives as a structural node inside the Mul-with-inverse so L148's `.subs` matches. (c) **variable-independence trap** — `D[flat,rV]`/`sp.diff(fl,rV)`: `fl` genuinely depends on `rV`, derivative nonzero (output shows nonzero numerator), no zero-derivative trap. (d) **hardcoded benchmark** — M9/section-6 numerics pin computed values against the documented Session-IV targets; `E_v=0.25*0.01033460` checked against `0.00258365` is a calibration-consistency pin (the notes explicitly label it "a calibration consistency check, not a theorem"), acceptable. (e) **symbol assumptions** — all Gammas nonnegative, integration vars positive; positivity is justified by the physical setup (`Gamma_n^vac, Gamma_n^lat >= 0`, `s,T,Vin>0`). I read the card, notes, and appendix; the script's verified claims match the paper's seven deliverables exactly.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `f_lat = (Gamma3lat+Gamma5lat r_V)/(Gamma3+Gamma5 r_V)` | py L59 / wl L13 / out L9 | card eq L40-41, notes L166-169, appendix eq L323-324 | MATCH |
| `df_lat/dr_V` drift form | py L67 / wl L24 / out L15 | card eq L46-48, notes L181-184 | MATCH |
| endpoints `Gamma3lat/Gamma3`, `Gamma5lat/Gamma5` | py L76-77 / wl L31-32 / out L16-17 | notes L197-199 | MATCH |
| `I2 = Vin^2 s^3(e^{2sT}-1)/2`, `I3=s^2 I2`, `r_V=s^2` | py L104-105 / wl L43 / out L23-25 | card eq L61-65, notes L227-240 | MATCH |
| `gamma^eq = Gamma3 s^2 + Gamma5 s^4` | py L94 / wl L143-144 / out L28 | card eq L70-72, notes L257-264 | MATCH |
| `E_exp,min^safe = (Vin^2/2)(e^2-1)mu_eta(s0^2-sc^2)` | py L139 / wl L171 / out L36 | card eq L74/L323-331(appendix), notes L322-325 | MATCH |
| `gamma_eff,safe^eq = mu_eta(s0^2-sc^2)/sc` | py L122,148 / wl L186 / out L37 | notes L360-364 | MATCH |
| 3:1 surface `Gamma3lat+Gamma5lat sc^2 = 3(Gamma3vac+Gamma5vac sc^2)` | py L155,172 / wl L193 / out L44 | card eq L81-84, notes L389-395 | MATCH |
| phi-family `f_lat=phi`, `phi=3/4` | py L179,191 / wl L230,235 / out L45 | card eq L410-418, notes L409-419 | MATCH |
| `sc ≈ 0.5489386551`, `sc^2 ≈ 0.3013336471` | py L223-224 / wl L263 / out L50-51 | notes L441-443 | MATCH |
| `gamma_eff,safe^eq ≈ 87.26925235` | py L225 / wl L264 / out L52 | notes L455 | MATCH |
| `gamma_vac/lat,safe ≈ 21.81731309 / 65.45193926` | py L226-227 / wl L265-266 / out L53-54 | notes L461-463 | MATCH |
| `Vin_match ≈ 8.21771260e-3` | py L229 / wl L268 / out L56 | notes L480 | MATCH |
| `E_vac/E_lat ≈ 0.00258365 / 0.00775095` | py L230-231 / wl L269-270 / out L57-58 | notes L485-488 | MATCH |
| (card verification line "Mathematica audit: none yet") | wl exists, M1-M9 PASS | card L4 | MISMATCH → F1 |

INTERNAL (no finding): `I1_exp` velocity integral, `safe_combo`, `Et_safe_expected`/`Et_safe_reduced`, `safePrefactorNum=153.0353549...`, residual/PASS flags, tolerances, section banners.

reconciliation: complete; 14 deliverable values checked, 0 numeric-value mismatches; 1 card-text-lag mismatch (verification line) folded into F1.

## Self-test notes

Checked: (1) variable-independence — `sp.diff(fl,rV)`/`D[flat,rV]` depend genuinely on rV (nonzero output), no identically-zero-derivative trap. (2) round-trip re-check (pass-1 flag) — confirmed the load-bearing rate identity (py L140 / wl M7-L182 `safe_combo/sc == Gamma3 sc^2+Gamma5 sc^4`) is non-tautological and can fail; the bridge subs (py L148) genuinely fires because the `G3 sc^3+G5 sc^5` Add survives as a structural node under division by sc. (3) subs-fired — sympy output L36 confirms the M6 `safe_combo` substitution reduced to the mu_eta form. (4) trivial-case — fractions sum to 1 and benchmark numerics reproduce the documented Session-IV split exactly. Conclusion: pass-1 tautology is resolved; only a low-severity card-text-lag remains.
