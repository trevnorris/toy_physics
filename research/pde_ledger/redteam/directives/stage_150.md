---
unit_id: 150
batch: IV.5
created_at: 2026-05-29T00:00:00Z
supersedes: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-29T22:46:12-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 150 (rewrite encoding Codex review of 2026-05-29)

## What this rewrite does

The 2026-05-27 directive (now superseded) prescribed a SOURCE-level de-tautology of the
slope: replace `Sq = simplify(diff(Tq,x).subs(x,0))` with the hand-derived closed form
`Sq = Aq*k - Cq*Pi` (SymPy) / `sQ = aq*k - cq*p` (Mathematica). **That source fix is
already applied and committed (3e2b5c0)** — see py:37-39 and wl:37-39, both already in the
hand-derived form. It is correct and **must not be re-prescribed** (no-op). The downstream
assertion `T_q'(0)-S_q == 0` (py:47 / wl:46) is no longer tautological and passes.

The Codex review (`redteam/codex_reviews/stage_150.md`, findings_count = 1) found exactly
**one remaining issue, R1 — a DISPLAY/corroboration finding on the regenerated
TRANSCRIPTS, not on the source**. The directive originally required the regenerated
transcripts to DISPLAY a compact slope `S_q(Pi) = Aq*k - Cq*Pi` / `aq*k - cq*p` as visible
evidence that the old expanded-derivative display was removed. Instead both transcripts
still print a fully substituted rational at line 5, because `Aq`/`Cq` (py:33-34) and
`aq`/`cq` (wl:33-34) are concrete expressions, so `sp.pprint(Sq)` (py:42) and `fmt[sQ]`
(wl:41) auto-expand. The math is right; only the printed corroboration is missing.

Apply the single finding below. After applying, append an `## Applied: F1` block under it
with: `files_changed`, `summary` (one sentence), and `deviation` (or "none"). If the
required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with
a question instead.

Do NOT introduce features/refactors/stylistic changes beyond the named edits. Do NOT touch
paper.tex, notes/, or any prose document. After editing, RUN the affected scripts
(`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they
exit 0 with all in-file checks passing, all within `timeout 600` (a timeout is a failure —
reformulate, never raise the cap). Getting the scripts to run cleanly is your job; the
orchestrator independently re-runs afterward.

NOTE: there is NO `sp.integrate` / `Integrate` / `NIntegrate` anywhere in either Stage 150
script — there is no hang or quadrature-precision risk in this fix.

---

## F1 — insufficient_verification (DISPLAY corroboration of the compact slope)

**Targets (LIVE anchors — old directive cited line 37 / line 45; those drifted by +2 after
3e2b5c0, the assert is now py:47/wl:46 and the assignment py:39/wl:39):**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py` — display block py:41-42.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl` — display line wl:41.

**Already applied (3e2b5c0), no-op:** the source slope is already the hand-derived closed
form `Sq = Aq*k - Cq*Pi` (py:39) / `sQ = aq*k - cq*p` (wl:39). Do NOT re-DERIVE the slope by
a different route and do NOT change its VALUE. (The display fix in F1 below restructures the
*assignment* via free placeholder symbols so the printed form can be compact, but the
resulting `Sq`/`sQ` object stays mathematically identical — see the `.subs` step. Keep the
explanatory comments at py:37-38 / wl:37-38.)

**Issue (review R1):** the regenerated transcripts must DISPLAY the compact symbolic slope
`S_q(Pi) = Aq*k - Cq*Pi` (SymPy) and `S_q(Pi) = aq*k - cq*p` (Mathematica), optionally
followed by the expanded value. Currently `transcript line 5` prints the fully substituted
rational instead (confirmed: py output line 5 is the expanded fraction; wl output line 5 is
`-(p^2/((1 - E^(-p))*(-p^2 + Pi^2/4))) + ...`). The auto-expansion happens because the
slope is built from the concrete `Aq`/`Cq` (`aq`/`cq`) definitions.

**Required change — consult-Q6 approach (b): print from FREE placeholder symbols, then
`.subs` the concrete definitions for ALL load-bearing checks.** Build the printed compact
form from free coefficient symbols so the displayed `Aq*k - Cq*Pi` (`aq*k - cq*p`) is
*provably the real slope* (it is the same symbolic object, with the concrete coefficient
expressions substituted in for the assertion). Do NOT print a decoupled hardcoded display
string — that would be a fabricated-display anti-pattern (the printed text would no longer
be guaranteed equal to what the assertion uses).

### (1) SymPy — `scripts/...stage150...sympy_audit.py`

Introduce two free placeholder symbols for the coefficients, define the slope COMPACTLY
from them, print that compact form, then substitute the concrete `Aq`/`Cq` definitions to
recover the load-bearing `Sq`.

Before (py:39 assignment, py:41-42 display):
```python
Sq = Aq*k - Cq*Pi

print("S_q(Pi) =")
sp.pprint(Sq)
```

After (replace py:39 + py:41-42; keep the comment at py:37-38 above):
```python
# Build the slope compactly from FREE coefficient symbols so the PRINTED form is provably
# the real slope; then .subs the concrete Aq, Cq definitions for the load-bearing checks.
Aq_s, Cq_s = sp.symbols("Aq Cq")
Sq_symbolic = Aq_s*k - Cq_s*Pi
Sq = Sq_symbolic.subs({Aq_s: Aq, Cq_s: Cq})

print("S_q(Pi) =", Sq_symbolic)        # compact: Aq*k - Cq*Pi
print("S_q(Pi) [expanded] =")
sp.pprint(Sq)
```

`Sq` (the substituted object) remains exactly the load-bearing slope used at py:47; only the
DISPLAY now shows the compact `Aq*k - Cq*Pi` first. (`k = sp.pi/2`, so the printed form is
`Aq*pi/2 - Cq*Pi` — still the compact coefficient form, NOT the expanded rational. If a
literally `Aq*k - Cq*Pi` token is preferred over `Aq*pi/2`, leave `k` symbolic in the print
by using a free `k_s` symbol too; either is acceptable, the requirement is that `Aq`,`Cq`
appear as symbols, not as their expanded definitions.)

### (2) Mathematica — `mathematica/...stage150...mathematica_audit.wl`

Mirror with free placeholder symbols `aqS`, `cqS`; print the compact form, then substitute
the concrete `aq`/`cq` definitions for `sQ`.

Before (wl:39 assignment, wl:41 display):
```mathematica
sQ = aq*k - cq*p;

Print["S_q(Pi) = ", fmt[sQ]];
```

After (replace wl:39 + wl:41; keep the comment at wl:37-38 above):
```mathematica
(* Build the slope compactly from FREE coefficient symbols so the PRINTED form is provably
   the real slope; then substitute the concrete aq, cq definitions for the load-bearing checks. *)
sQsymbolic = aqS*k - cqS*p;
sQ = sQsymbolic /. {aqS -> aq, cqS -> cq};

Print["S_q(Pi) = ", fmt[sQsymbolic]];        (* compact: aq*k - cq*p *)
Print["S_q(Pi) [expanded] = ", fmt[sQ]];
```

`sQ` (the substituted object) stays exactly the load-bearing slope used at wl:46. `k = Pi/2`
(wl:31), so the compact print shows `aqS*(Pi/2) - cqS*p` form with the COEFFICIENT symbols
preserved; if a literal `aq*k - cq*p` token is preferred, also carry `k` as a free symbol
`kS` in `sQsymbolic`. Either is acceptable — the requirement is `aq`,`cq` appear as symbols.
(Mathematica-idiom reminder: no `*)` inside comment text; the comments above are prose-safe.)

**HARD guards (state explicitly — do NOT violate):**
1. The printed compact form MUST be the SAME symbolic object used in the load-bearing
   assertion (built via free symbols then `.subs`/`/.` to the concrete coefficients), NOT a
   decoupled hardcoded string. A divergence between what is printed and what is asserted is
   itself a finding.
2. The assertion `T_q'(0) - S_q == 0` (py:47 / wl:46) MUST remain UNCHANGED and still PASS.
   `Sq` (py) and `sQ` (wl) must continue to equal the concrete hand-derived slope.
3. Do NOT alter the curvature checks: `R(0)` (py:51 / wl:49), `R'(0)` (py:52 / wl:50),
   `R''(0)` computation + `R''(0) - target` (py:55-59 / wl:52-55). These are load-bearing
   and unaffected by a display edit.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 150` and
`redteam exec-mathematica 150`. Both scripts must:
- (a) exit 0;
- (b) at `transcript line 5` show the COMPACT coefficient form — `Aq*k - Cq*Pi` (i.e.
  `Aq*pi/2 - Cq*Pi`) in the SymPy transcript and `aq*k - cq*p` (i.e. `aqS*(Pi/2) - cqS*p` /
  `aq*k - cq*p`) in the Mathematica transcript — with `Aq`/`Cq` (`aq`/`cq`) appearing as
  symbols, optionally followed by an expanded-value line; NOT the fully substituted rational;
- (c) still report `T_q'(0)-S_q = 0` (SymPy) and `T_q'(0)-S_q = 0` + `PASS: T_q'(0)-S_q`
  (Mathematica);
- (d) still report `R''(0) - target = 0` (SymPy) and `R''(0) - target = 0` +
  `PASS: R''(0) - target` (Mathematica), with `R(0)`, `R'(0)` also still zero/PASS.

**Claim manifest (what this finding does and does not certify):**
- CERTIFIES (after fix): the regenerated transcripts visibly corroborate that the slope
  display is the compact closed form `Aq*k - Cq*Pi` / `aq*k - cq*p`, and that printed form is
  provably the same object entering the `T_q'(0)-S_q` assertion (no fabricated display).
- DOES NOT change: the source slope (already correct, 3e2b5c0), the load-bearing assertion,
  or any curvature check. This is a display/corroboration strengthening only.

**Self-test (must hold after the edit):**
- The printed compact `S_q(Pi)` string contains the symbols `Aq` and `Cq` (`aq` and `cq`),
  NOT `exp`/`sinh`/`cosh`/`E^`/`Sech` (the signature of the expanded rational).
- Deleting the `.subs`/`/.` step (so `Sq`/`sQ` stays the free-symbol form) would make
  `T_q'(0)-S_q == 0` FAIL — proving the printed form and the asserted object are linked, not
  decoupled. (Do not actually delete it; this is the conceptual can-fail check.)

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl`
- summary: Added free coefficient placeholder slope expressions for the compact display, then substituted the concrete coefficient definitions for the existing load-bearing checks.
- deviation: none

---

## RESOLVED (consult batch 7)

Per `redteam/codex_reviews/_consult_batch7.md` **Q6** (150-R1, compact display — CONCUR):

> Source slope `Sq=Aq*k-Cq*Pi` already committed; only the DISPLAY expands (Aq/Cq
> concretized). Print a compact placeholder form (build from free `Aq,Cq,k`), THEN `.subs`
> the concrete definitions for the load-bearing `T_q'(0)-S_q==0` assert, so the printed form
> is provably the real slope. **Approach (b).**

This directive encodes exactly that: free placeholder coefficient symbols → compact print →
substitute concrete definitions for the load-bearing checks. The decoupled-hardcoded-string
route (approach a) was explicitly rejected as a fabricated-display anti-pattern. No user
escalation; this is a how-it's-displayed fix resolved Claude+Codex.

---

## Orchestrator note (NOT a script finding — do not fix in this directive)

`notes/stages/review/stage_150_review.md` is MISLABELED: its body is Stage 031 content, not
Stage 150. This is a prose/notes mislabeling for the orchestrator to handle separately; it is
outside the scripts-only fix loop and must NOT be edited here.
