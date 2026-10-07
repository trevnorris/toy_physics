# Light-leakage inventory v4: scoped-check dispositions (orchestrator, 2026-10-07)

**Artifact:** `light_leakage_scoping.md` v4. `sha256sum -c scoping_v4.sha256` → `light_leakage_scoping.md: OK`.
It was made by a fresh Claude author from v3 (frozen `light_leakage_scoping_v3.md`, `85f1d09c…`). The diff
touches lines 26, 27, 41, 71, 86, 92, 94 and 132.

**Legs (identical prompt `light_leakage_scoping_review_prompt_r4.md`):**
- Grok `scoping_review_r4_grok.txt`: **CLEAR**, no findings.
- Codex `scoping_review_r4_codex.txt`: **NEEDS REVISION**, one finding.

Both reported before adjudication. Both legs confirm that the F1 routing (to C3, with the gaps named), the
narrowed O7 and the F2 domain are faithful to the sources.

## Finding (Codex): the added amplitude-order rationale excludes too much

```
$ sed -n '41p' _scratch/light_leakage/light_leakage_scoping.md   (excerpt)
… Both point to C3, the generic nonlinear program, whatever the initiation … or WB gain class: on held profiles,
a coupling of transverse light to another wave perturbation enters as a product of wave amplitudes, beyond the
first-order truncation (this inventory's reading of SPa:194–198). …
$ sed -n '26p' … | grep -o '…first order in wave amplitude…'
… and, on this inventory's reading, a coupling of transverse light to another wave perturbation lies beyond the
records' "first order in wave amplitude" truncation (SPa:194–198; taxonomy note below)
$ sed -n '57p' …   (E1, unchanged)
| E1 | **Equation-level route:** nonuniform transverse↔`{θ,e_W,u_L}` operator, including gradient-thickness leg …
$ sed -n '123,126p' research/pde_ledger_v3/directives/S11c_decisions.md
… the transverse↔thickness coupling is `O(εη)`; the leakage probability/rate `O(ε²η²)`.
$ sed -n '366,368p' research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
… that energy is bilinear in the perturbations, so its first variation is linear in them …
```

**Verified.** The new clause says *any* coupling of transverse light to another wave perturbation is beyond
the first-order truncation. The records disagree: they put the transverse↔thickness coupling at first order in
the wave amplitude, `O(εη)` (DEC:125). The inventory's own E1 route (line 57) is exactly that first-order
interbranch operator. So the clause contradicts an unchanged route, and the bracket could misread E1's linear
conversion as nonlinear C3 work. That changes a route classification, which is inside the physics filter.

The v4 author introduced the clause unprompted, as "this inventory's reading", so this repair bred a defect in
the material it changed. The C3 routing itself is not in question: both legs found it sound.

**What must be true.** The rows 3–4 routing to C3 does not rest on any claim that couplings between wave
branches as such lie beyond first order. Nothing in rows 3–4 or the taxonomy note contradicts E1's first-order
interbranch operator (DEC:125). No v3 vertex is invented.

**Outcome under the user's rule:** not clear. One verified finding remains, so the decision goes to the user.

## User decision (2026-10-07)

"Go with your recommendation." That is:
1. A different fresh Claude author removes the over-broad claim. v4 is preserved as `light_leakage_scoping_v4.md`
   (`214d1fe6…`).
2. Codex and Grok check only that change (Claude-written → Codex + Grok), with an identical prompt.
3. The table is accepted if nothing outstanding changes what the bracket would compute or claim.

**Scope of the repair.** Change only what the finding requires, at v4 lines 26 and 41, plus any cell that would
otherwise contradict the change. Every other cell stays byte-identical. The directive still binds: no
calculations, no new physics, no verdicts; record disagreements quoted, not resolved; unopened sources marked
UNVERIFIED; cite only records and published sources.

## v5 scoped check and acceptance (orchestrator, 2026-10-07)

**Artifact:** v5, sha256 `17293320…` (`scoping_v5.sha256`; `sha256sum -c` → OK). The fresh author changed lines
26, 27 and 41 of v4 (`diff … | grep '^[0-9]'` → `26,27c26,27`, `41c41`). The report is
`scoping_repair5_claude_author_report.md`.

**Legs (identical prompt `light_leakage_scoping_review_prompt_r5.md`):**
- Codex `scoping_review_r5_codex.txt`: **CLEAR**, no findings;
- Grok `scoping_review_r5_grok.txt`: **CLEAR**, no findings.

Both legs confirm four things:
- the removed claim is absent;
- the note now affirms E1's `O(εη)` coupling (DEC:125; SC:48);
- the "homogeneous" qualifier matches S11:74–76;
- rows 3–4 stay NOT ADDRESSED with no vertex asserted.

**Author flags, adjudicated:**
- *No record line puts spontaneous Raman/Brillouin in the nonlinear program.* True. But the C3 pointer is a
  destination for a NOT ADDRESSED gap, not a mechanism claim. The note says "This choice does not decide whether
  a nonlinear or time-dependent interaction exists", and both legs read it as a pointer. The bracket would
  compute nothing from it beyond listing the gap. **Not outstanding.**
- *(v3 author)* c1 has no "Changed after close" note, though b:39 says the inertia repair propagated. This is a
  records-housekeeping item for the S11c step records, not for this inventory. **Not outstanding here.** It is
  carried to the S11c review bucket.

**Outcome under the user's rule: ACCEPTED (v5).** Nothing outstanding changes what the bracket would compute or
claim. Round history: 19 → 16 → 8 → 2 → 1 → 0 findings.
