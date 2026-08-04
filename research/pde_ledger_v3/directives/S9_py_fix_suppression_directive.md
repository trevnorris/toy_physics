# S9 SymPy engine — ⛔ URGENT: undo cross-package payload suppression

⛔ **One targeted fix. Change nothing else.**

**File (absolute):**
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py`

Verify:
`cd /var/projects/toy_physics/research/pde_ledger_v3 && timeout 600 python3 scripts/S9_light_requires_shear_sympy_audit.py`

⛔⛔ **Do not read any `.wl`, the registry, or any other engine.** ⛔ Do not commit.

---

## The defect

The script currently suppresses a tag whose payload is identical to one already emitted. That suppression
is being applied **across control packages**, and it is destroying the most important signal in the
output.

⚠ **Measured now:** `DIM_COEFFICIENTS` is emitted for `MAIN` but is **absent for X1, X2, X3, X4 and X5**,
because those controls produce the same coefficient dimensions as `MAIN`.

⛔⛔ **THIS INVERTS THE MEANING OF THE OUTPUT.** An automated consumer pairs each main tag with its
counterpart under every control and asks *"did this value move?"* Under that question:

- **a tag present with an identical value** = **INVARIANT** — a real, interpretable result;
- **a tag that is absent** = **indistinguishable from "no computation ever produced it"** — which is the
  exact defect this entire rebuild exists to detect.

⇒ ⭐ **Suppressing an identical payload across packages deletes the evidence that the quantity was
computed at all**, and it does so precisely for the quantities that are correctly invariant.

⚠ **The instruction that caused this was mine and it was under-specified.** It meant *"do not emit the
same object twice under two names **within one package**"* — e.g. one tag that is a verbatim alias of
another. ⛔ It did not mean *"deduplicate across controls."*

## The fix

1. ⭐ **Remove the cross-package suppression entirely.** **Every control package emits the full tag set,
   always**, whatever its values are. ⛔ A payload identical to another package's is emitted anyway —
   that identity **is** the finding.
2. ⭐ Keep de-duplication only in its intended, narrow sense: ⛔ **within a single package**, do not emit
   one computed object under two different tag names. If such an alias exists, delete the alias — ⛔ do
   not delete the original, and ⛔ do not make emission conditional on a value.
3. ⛔ **Emission must never be conditional on a payload's value.** Whether a tag appears may depend only
   on which package and which quantity it belongs to, ⛔ never on what the value turned out to be.
   ⚠ A value-dependent emission means the *presence* of a tag carries hidden information, and a consumer
   cannot tell a missing computation from a repeated one.

---

## Report back — under 15 lines

1. Exit status, runtime, total tag count, confirmation all tag **names** are unique and there is no
   untagged output.
2. Confirmation that **every** control package now emits the same tag set as `MAIN`, and the literal
   `DIM_COEFFICIENTS` line for `MAIN` and for each of `X1`–`X5`.
3. What you removed, by line number.
4. Anything that surprised you.

⚠ Note that tag **names** must remain unique while **payloads** may now legitimately repeat. ⭐ That is
the intended state.
