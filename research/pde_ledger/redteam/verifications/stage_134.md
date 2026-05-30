---
unit_id: 134
batch: IV.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 134 (REMEDIATION verify)

This supersedes the stale 2026-05-27 verification (which pre-dated the recovery review
and wrongly recorded the gain-line `OK:` line at py:82 as a PASS — the very check R1
later flagged as tautological). The authoritative finding set is
`redteam/codex_reviews/stage_134.md` (verdict FINDINGS, one live finding R1 —
tautological_check). F1/F2/F4 are `tainted-applied` (their checks PASSED the review and
are retained unchanged); F3 is RESOLVED via the IV.4 paper-card downgrade (no script
action). The only live edit verified here is R1.

## Per-finding outcomes

### R1 — tautological_check (canonical gain-line assertion)

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:81-90`:
  the entire "Check 4: canonical gain line" block is DELETED. Formerly it asserted
  `intercept = sp.N(Pi_star,30)` against the literal `1.50882951349316` that *defined*
  `Pi_star` at line 39 (an X−X comparison), plus a `slope = -S_star` vs
  `-0.658075937605428` that merely restated Check 3 with a sign flip. In its place is a
  no-assertion provenance comment: the intercept is the imported literal Π_* (owned by
  stage 131), the slope is −S_q(Π_*) already validated in Check 3, re-asserting them
  would be an X−X tautology so it is intentionally omitted, and the substantive
  deliverable (outlet consistency of the gain pair) is verified downstream at Stage 135
  (susceptibility closure at Stage 137). There is NO `assert` anywhere in this block.
- The symbolic `print("Canonical gain line: ...")` + `sp.pprint(...)` at lines 42-43
  remains (a print, not an assertion) — correct per directive.
- Mathematica `mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:66-77`:
  the `.wl` block (which never had an `expectZero`/`pass`/`fail`, so was not a hard
  runtime tautology) is relabeled `"Canonical gain line Ms = Pi_star - S_q(Pi_star) M_q
  (printed only; not asserted)"` (line 71); the prior "matches the Stage 134 note" claim
  in the RESULT block is replaced with the downstream-Stage-135 outlet-consistency
  statement (line 84). A no-assertion provenance comment was added at lines 74-77.

**Assessment:**
The fix REMOVES the tautological check rather than dressing it up as a fake
"independent" gain-line check — exactly what the codex_review prescribed once
outlet/gain-line consistency was downgraded to a Stage-135 carry-forward. The gain line
is not "moved" to another self-referential assert; it is dropped entirely and the
deliverable is explicitly deferred downstream. No new numeric literal was fabricated;
the two deleted targets (`1.50882951349316`, `-0.658075937605428`) are gone and nothing
replaces them.

The three KEPT substantive checks are intact and non-tautological — each compares the
script kernel against an *external* source of truth (so none is X−X):
- Check 1 (py:47-51 / wl:44): shell limit `S_shell - 1 == 0` (external truth = exact `1`).
  Exec: `S_shell - 1 = 0` → `OK: S_shell = 1` / `PASS: static shell channel`.
- Check 2 (py:53-72 / wl:57-62): S_q(1/2), S_q(1), S_q(2) vs independent mpmath literals
  `0.608336…`, `0.633127…`, `0.681366…`. Exec: SymPy diffs ~3.9e-31 / 0 / 4.9e-31;
  Mathematica PASS×3.
- Check 3 (py:74-79): S_q(Π_*) vs notes value `0.658075937605428`. Exec:
  `OK: S_q(Pi_star) matches notes value 0.658075937605428`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `S_shell - 1 = 0` / `OK: S_shell = 1`
- `S_q(1/2) = 0.608... (target 0.608..., diff 3.94...E-31)` … `OK: S_q matches independent numeric targets at Pi=1/2, 1, 2`
- `OK: S_q(Pi_star) matches notes value 0.658075937605428`
- The former `OK: gain line coefficients match notes ...` line is ABSENT, as the directive requires.

**Mathematica:** exit=0. Notable lines:
- `PASS: static shell channel`
- `PASS: S_q at p=1/2` / `PASS: S_q at p=1` / `PASS: S_q at p=2` (diffs ~0)
- `Canonical gain line Ms = Pi_star - S_q(Pi_star) M_q (printed only; not asserted)`
- RESULT now reads `Outlet consistency of the gain pair (Ms, Mq) is verified downstream at Stage 135.`

**Output freshness:** confirmed. Saved outputs
`scripts/output/...sympy_audit.txt` and `mathematica/output/...mathematica_audit.txt`
both have mtime 2026-05-29 16:49:02, newer than the edited scripts (.py 16:41:15,
.wl 16:41:55). The saved `.txt` tails match the exec logs (no gain-line `OK`/assert
line; downstream-Stage-135 RESULT line present).

## Literal provenance (confirmed against owner files, not recomputed)

| literal | owner confirmed |
|---|---|
| `1.50882951349316` (Π_*) | `notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md:8` (carried at stage 134 notes:72) |
| `0.658075937605428` (S_q(Π_*)) | `notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md:86` |
| `0.608336415687717…` (S_q(1/2)) | M3 mpmath run, `redteam/resolutions/batch_IV4_paper_alignment.md:75-80` |
| `0.633127670034487…` (S_q(1)) | same M3 run |
| `0.681366857005321…` (S_q(2)) | same M3 run |

The two deleted Check-4 targets are NOT retained anywhere. No new literal was introduced.

## Material-change assessment

`material_change`: false.

A redundant, self-referential (X−X) assertion was removed and replaced with a
non-executing provenance comment; the `.wl` print was relabeled and its RESULT prose
re-pointed downstream. No derived result, no kernel, no retained numeric value, and no
emitted symbolic expression changed (the gain line was already only printed, not
asserted). The three substantive checks and their values are byte-identical to the prior
state. Downstream units therefore have no derived-result dependency on this edit.

## Side observations (non-blocking)

- Codex's one recorded deviation — the Mathematica comment uses `PiStar`/`piStar`
  instead of `Pi_*` — is benign and in fact necessary: the token `Pi_*)` would close a
  `(* ... *)` block comment prematurely. This is the documented `Pi_*)` comment-terminator
  pitfall; the substitution is correct and changes no logic.
- The previous (stale) verification file on disk recorded findings 4/4 resolved and cited
  the now-removed gain-line `OK:` line as a PASS. That artifact is overwritten by this
  remediation verdict; the 4-finding framing was the tainted 2026-05-27 directive's, not
  the authoritative codex_review's.
- `notes/` was read only to confirm literal provenance (explicitly required by the
  remediation directive's literal-ownership table); no prose-derivation re-audit was done.

## Verdict justification

R1 is fully resolved: the tautological gain-line assertion was removed — not moved, not
replaced with another self-referential check — a no-assertion provenance comment was put
in its place, and the gain-line/outlet deliverable was correctly deferred to Stage 135 in
line with the F3 paper-card downgrade. No new numeric literal was fabricated; all retained
literals trace to their owners (Π_*←stage 131, S_q(Π_*)←stage 134 notes, the three S_q
targets←the M3 mpmath run). The three retained checks (shell limit, three S_q mpmath
spot-checks, S_q(Π_*) vs notes) are intact, non-tautological, and pass in both engines
(sympy_exit=0, mathematica_exit=0). Saved outputs were refreshed post-fix. material_change
is false. Verdict: verified.
