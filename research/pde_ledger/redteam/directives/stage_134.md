---
unit_id: 134
batch: IV.4
created_at: 2026-05-27T00:00:00Z
rewritten_at: 2026-05-29T00:00:00Z
rewrite_basis: redteam/codex_reviews/stage_134.md (authoritative; supersedes the 2026-05-27 tainted-fix prescription)
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-29T22:40:29Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 134 (REWRITTEN against codex_review)

This directive was rewritten to encode the authoritative `redteam/codex_reviews/stage_134.md`.
The original 2026-05-27 directive PRE-DATED that review and prescribed a "canonical gain line"
check that the review later flagged as still-tautological. The original F1/F2/F4 are reconciled
as `tainted-applied` (the checks they introduced PASSED the review and stay), F3 is RESOLVED
(see block below), and the single live finding is R1 — the tautological gain-line assertion.

Apply the live finding below. After applying, append an `## Applied: R1` block under it with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the lines named.
Do NOT introduce ANY new numeric literal. Do NOT generate a target by evaluating the script's
own kernel (`S(Pi, pi/2)` / `sKernel[p, Pi/2]`) at runtime — that re-introduces the tautology.

Do NOT touch paper.tex, notes/, or any prose document. The red-team only modifies scripts.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for
Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts
to run cleanly is your job; the orchestrator independently re-runs afterward.

---

## R1 — tautological_check (the only live finding)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:81-91`
(and the print-only mirror in `mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:66-78`)

**Issue (from codex_review R1, `..._sympy_audit.py:82`):**
The "Check 4: canonical gain line" block asserts
`intercept = sp.N(Pi_star, 30)` equals the literal `sp.Float("1.50882951349316", 30)` — but
`Pi_star` was itself *defined* as `sp.Float("1.50882951349316")` at line 39, so this is an
X−X comparison of a literal against itself. The companion assert `slope = -S_star` vs
`-0.658075937605428` only restates Check 3 (`S_q(Pi_*) ≈ 0.658075937605428`, lines 76-79) with
a sign flip. The block CANNOT detect a wrong gain-line derivation or an unsupported gain
selection; it only confirms constants already inserted into the script. The review's prescribed
fix: "Either remove this as a substantive check, or compare an independently derived / stored
gain-line expression or outlet-consistent gain pair against the fixed-point residual. If outlet
consistency belongs to Stage 135 as the current card says, Stage 134 should not claim this as an
independent gain-line verification."

**Resolution chosen (consistent with F3 below):**
Outlet consistency and the gain pair `(M_s, M_q)` were moved to Stage 135 (and susceptibility
closure to Stage 137) by the IV.4 paper-card downgrade — see `## RESOLVED: F3`. Stage 134
therefore must NOT claim an independent gain-line verification. NEUTRALIZE the tautological
assertion: delete the self-referential intercept/slope asserts and replace them with a
provenance comment that the gain-line / outlet-consistency is verified at Stage 135 — not a
self-referential numeric equality. The fixed-point line `Pi = M_s + M_q*S_q(Pi)` is still
printed for the transcript (it is symbolic, not an assertion), and the genuinely substantive
numeric anchors (Checks 1, 2, 3) stay intact.

**Required change (SymPy — exact replacement):**

Delete lines 81-91 (the entire "Check 4" block, from the comment
`# Check 4: canonical gain line ...` through the final
`print("OK: gain line coefficients match notes ...")`), i.e.:

```python
# Check 4: canonical gain line M_s = Pi_star - S_q(Pi_star)*M_q.
intercept = sp.N(Pi_star, 30)
slope = sp.N(-S_star, 30)
intercept_target = sp.Float("1.50882951349316", 30)
slope_target = sp.Float("-0.658075937605428", 30)
assert abs(sp.N(intercept - intercept_target, 30)) < sp.Float("1e-12"), \
    f"gain line intercept mismatch: got {intercept}, want {intercept_target}"
assert abs(sp.N(slope - slope_target, 30)) < sp.Float("1e-12"), \
    f"gain line slope mismatch: got {slope}, want {slope_target}"
print("OK: gain line coefficients match notes (intercept 1.50882951349316, "
      "slope -0.658075937605428)")
```

and replace it with EXACTLY this provenance comment (no assertion, no new literal):

```python
# Note (no in-stage gain-line assertion): the canonical gain line
#   M_s = Pi_* - S_q(Pi_*)*M_q
# is printed symbolically above for the transcript only. The intercept is the
# imported literal Pi_* (owned by stage 131; see stage 131 note) and the slope
# is -S_q(Pi_*), already validated against the notes value in Check 3 above.
# Re-asserting intercept == Pi_* and slope == -S_q(Pi_*) here would be an X-X
# tautology (it compares constants already inserted into the script), so it is
# intentionally omitted. The substantive deliverable that uses this gain line --
# outlet consistency of the gain pair (M_s, M_q) -- is verified downstream at
# Stage 135 (outlet-consistent mouth closure); susceptibility closure at Stage 137.
```

Notes for Codex:
- Lines 1-79 are UNCHANGED. In particular keep Check 1 (shell limit, lines 47-51), Check 2
  (S_q at Pi=1/2,1,2 vs mpmath literals, lines 53-72), and Check 3 (S_q(Pi_*) vs notes value
  0.658075937605428, lines 74-79) exactly as they are — the review marked all three PASS.
- The symbolic `print("Canonical gain line: ...")` / `sp.pprint(...)` at the top of the file
  (lines 42-43) is a print, not an assertion — leave it as-is.
- After deletion, `Pi_star` (line 39) and `S_star` (line 40) remain referenced by Check 3, so
  do NOT delete those definitions.

**Anti-tautology guard:** The replacement contains NO assertion at all for the gain line, so it
cannot be an X−X check. The only retained numeric assertions are Checks 1-3, each of which
compares the script kernel against an *external* source of truth (exact `1`; mpmath literals;
the stage-134 notes value 0.658075937605428) — none compares a literal against the same literal
used to define it. The gain-line deliverable is explicitly deferred to Stage 135, matching the
paper card.

**Required change (Mathematica — mirror, exact replacement):**

The `.wl` mirror at lines 66-78 only PRINTS the gain line (it has no `expectZero`/`pass`/`fail`
assertion), so it is not a hard runtime tautology — but it still presents a self-referential
"matches the Stage 134 note" claim and should be brought into alignment with the SymPy side.
Replace lines 66-78:

```mathematica
piStar = SetPrecision[1.50882951349316, 30];
sStar = N[sQ /. p -> piStar, 30];
gainLine = N[piStar - sStar*Mq, 30];

Print["S_q(Pi_star) = ", sStar];
Print["Canonical gain line Ms = Pi_star - S_q(Pi_star) M_q"];
Print[gainLine];

Print[""];
Print["RESULT:"];
Print["  Family-1 reduces the coupled mouth law to p = Ms + Mq S_q(p)."];
Print["  The static shell lane gives S_shell = 1 exactly."];
Print["  Evaluated at Pi_star, the canonical compensation line matches the Stage 134 note."];
```

with:

```mathematica
piStar = SetPrecision[1.50882951349316, 30];   (* imported literal; owned by stage 131 note *)
sStar = N[sQ /. p -> piStar, 30];
gainLine = N[piStar - sStar*Mq, 30];

Print["S_q(Pi_star) = ", sStar];
Print["Canonical gain line Ms = Pi_star - S_q(Pi_star) M_q (printed only; not asserted)"];
Print[gainLine];

(* No in-stage gain-line assertion: intercept is the imported literal Pi_*
   (stage 131) and the slope is -S_q(Pi_*); re-asserting them here would be an
   X-X tautology. Outlet consistency of (Ms, Mq) is verified at Stage 135;
   susceptibility closure at Stage 137. *)

Print[""];
Print["RESULT:"];
Print["  Family-1 reduces the coupled mouth law to p = Ms + Mq S_q(p)."];
Print["  The static shell lane gives S_shell = 1 exactly."];
Print["  S_q at p=1/2,1,2 and at Pi_star match independent literal targets (see checks above)."];
Print["  Outlet consistency of the gain pair (Ms, Mq) is verified downstream at Stage 135."];
```

Notes for Codex (.wl):
- Lines 1-62 are UNCHANGED — keep `expectZero["static shell channel", ...]` (line 44) and the
  three `expectClose[...]` numeric checks (lines 57-62) exactly. The review marked the shell and
  S_q numeric checks PASS.
- `sStar` is still computed and printed (transcript continuity); only the *claim* that it
  "matches the Stage 134 note" is removed, since there is no assertion backing it here and the
  SymPy side carries the real Check-3 assertion.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 134` and `redteam exec-mathematica 134`.
Expect both to exit 0. SymPy transcript: the `OK: gain line coefficients match notes ...` line is
GONE; `OK: S_shell = 1`, `OK: S_q matches independent numeric targets at Pi=1/2, 1, 2`, and
`OK: S_q(Pi_star) matches notes value 0.658075937605428` still appear. Mathematica transcript:
`PASS: static shell channel`, `PASS: S_q at p=1/2`, `PASS: S_q at p=1`, `PASS: S_q at p=2` still
appear; the RESULT line now references downstream Stage 135 for outlet consistency.

## Applied: R1

- files_changed:
  - `scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl`
- summary: Removed the self-referential gain-line assertion/claim and replaced it with provenance text deferring outlet consistency to Stage 135.
- deviation: Mathematica comment spelling uses `PiStar` instead of `Pi_*` where `Pi_*)` would terminate a `.wl` comment.

---

## RESOLVED: F3 (paper_misalignment) — already done in batch IV.4 Cluster C, via direction (b)

The original F3 held for user resolution on whether the paper card's two checklist items
("gain pair `(M_s, M_q)` against outlet consistency" and "self-matched susceptibility closure")
belong to stage 134 or downstream. The user chose direction (b): they belong to downstream
stages. This was APPLIED in an earlier batch (IV.4 Cluster C): the two `\item` lines were
DOWNGRADED to carry-forward citations of the owning stages.

Confirmed in `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_134.tex:21-25`,
the Checks block now reads:

```latex
\stagefield{Checks}{\begin{verificationchecklist}
\item Outlet consistency of the gain pair \((M_s,M_q)\) is verified at Stage~\ref{stage:135} (outlet-consistent mouth closure); carried forward here.
\item Self-matched susceptibility closure is verified at Stage~\ref{stage:137} (explicit core-to-mouth gain map); carried forward here.
\item Check numerical fixed points are recorded as numerically located, not closed-form constants.
\end{verificationchecklist}}
```

Consequences for the scripts (enforced by R1 above):
- Stage 134 must NOT add an outlet-consistency check or a susceptibility-closure check.
- Stage 134 must NOT claim an independent gain-line verification.
No paper-side or script-side action is required for F3 beyond the R1 neutralization. Do NOT
re-add the downgraded items.

---

## RECONCILED: F1 / F2 / F4 (tainted-applied — keep, no further change)

These were applied by the 2026-05-27 directive and the codex_review confirmed the resulting
checks are corroborated by the saved outputs and are NOT runtime tautologies:
- **F1 (SymPy substance):** Check 1 (shell limit `S_shell - 1 == 0`, lines 47-51), Check 2
  (S_q at Pi=1/2,1,2 vs mpmath literals, lines 53-72), Check 3 (S_q(Pi_*) vs 0.658075937605428,
  lines 74-79) — all PASS. Keep unchanged.
- **F2 (Mathematica de-tautologization):** the tautological `expectZero["specialized D/N kernel",
  sQ - sQExpected]` was deleted and replaced with the three `expectClose` numeric checks
  (lines 57-62) plus the retained shell `expectZero` (line 44) — all PASS. Keep unchanged.
- **F4 (transliteration):** subsumed by F1/F2 — both engines now compare `S_q` at independent
  numeric points against pasted literals sourced from a separate mpmath run, operationally
  neutralizing the shared-closed-form concern. No additional change.

Mark these `tainted-applied`: reconciled and retained. The only live edit is R1.

---

## Numeric literals retained in the script + owning source (verifier must confirm, not assume)

| literal | where used | owning source the verifier must check it against |
|---|---|---|
| `1.50882951349316` (Pi_*) | sympy lines 39, 84(removed); wl lines 66 | `notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md:8` (carried at `notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md:72`; original source-shape point Stage 233) |
| `0.658075937605428` (S_q(Pi_*)) | sympy line 76 (Check 3, retained); 85(removed) | `notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md:86` (this stage's own notes deliverable) |
| `0.608336415687717065435990381419` (S_q(1/2)) | sympy line 61; wl line 58 | mpmath run recorded in `redteam/resolutions/batch_IV4_paper_alignment.md` (M3) — NOT the script kernel |
| `0.633127670034487546375729566676` (S_q(1)) | sympy line 62; wl line 60 | same mpmath run, `batch_IV4_paper_alignment.md` (M3) |
| `0.681366857005321783286541952613` (S_q(2)) | sympy line 63; wl line 62 | same mpmath run, `batch_IV4_paper_alignment.md` (M3) |
| `-0.658075937605428` (slope target) | sympy line 85 — DELETED by R1 | (was the tautological half of Check 4; removed, not retained) |
| `1.50882951349316` (intercept target) | sympy line 84 — DELETED by R1 | (was the tautological half of Check 4; removed, not retained) |

The verifier MUST open the cited owner files and confirm each retained literal matches before
signing off. Do NOT regenerate any of these by running the kernel.
