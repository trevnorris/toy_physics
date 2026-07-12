# Design-review this computation directive BEFORE execution. READ-ONLY. Do NOT write or run code.

Read:
- `software/em_charge_attribute/directive_emergent_em_construction.md` (the directive under review — "Construction B")
- `docs/em_charge_attribute_requirements.md` (context: §2 requirements A/B/D, the reframe)

Context: A 3-AI panel judged the "EM = decoupled charge-attribute" reframe `conditional-needs-X` (not a no-go). Construction B is the prior test that the softness question is downstream of: does a local constraint on the constituents' charge-attribute actually produce emergent electromagnetism — a deconfined compact-U(1) with charge = the constraint DEFECT (not a postulated U(1) charge) and the DUAL SIGN (like-charges repel, like-currents attract) — with the model's ±w throat mapping onto the emergent charge. Grounded in quantum-spin-ice / quantum-dimer emergent-U(1) physics.

Assess, using real condensed-matter / gauge-theory knowledge:
1. **Is the DUAL SIGN the right decisive able-to-fail test** for "is it EM vs a scalar," and is the directive's demand that both signs come from ONE emergent propagator correct? Any subtlety (e.g. does the emergent compact-U(1) photon in spin ice actually give a Coulomb `1/r²` between emergent charges AND the current-current attraction, or are there caveats — magnetostatics vs full electrodynamics, static vs dynamic)?
2. **Is the anti-circularity mechanism sound** — is "charge = constraint defect (integer, emergent)" genuinely non-circular, and is the FAIL_CHARGE_POSTULATED check able to catch a construction that smuggles in the U(1)? Where is the real line between emergent and assumed here?
3. **Are the controls right** (no-constraint → no gauge; known-answer reproduces QSI/dimer Coulomb phase over a finite region)? Any vacuous/tautological one, or a missing control that would catch a wrong construction?
4. **Any physics error or hidden assumption** — e.g. does mapping a bare Z₂ ±w label onto an emergent integer charge actually work, or is there an obstruction? Is 3+1D handled right? Is the emergent photon's mode count / Gauss law correctly required?
5. **Is it tractable as a semi-analytic + small-numeric construction**, or does verifying the deconfined phase / dual sign secretly need a full lattice Monte Carlo? What is honestly analytic vs what must be cited from literature vs what needs numerics?
6. The single most important change before Codex executes it.

Output: terse, per-item 1–6, then a one-line verdict: `DIRECTIVE_READY` or `DIRECTIVE_NEEDS_FIXES` (list them). Your entire final message is the deliverable.
