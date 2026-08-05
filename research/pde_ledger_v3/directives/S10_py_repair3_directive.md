# S10 SymPy engine — repair round 3: one dead control

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py`
**Its specification:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
— ⭐ re-read **§4 (the structural rule)** and **Q7** in full before you touch anything.

⛔ **Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file.**
⛔⛔ **Do not read the sibling Mathematica engine** (`mathematica/S10_brane_mode_spectrum_*.wl`) or its
output under `mathematica/out/`, and ⛔ do not read `research/pde_ledger_v3/steps/` or
`research/pde_ledger_v3/paper/`. ⚠ The two engines exist so they can **disagree**; reading the sibling
destroys that and is worse than leaving the bug in.

⚠ Re-run after the fix and confirm it completes under **10 minutes** and exits `0`.

---

## ⛔⛔ R1 — BLOCKING, and it is the only substantive item: Q7's control does not re-enter at the action

⭐ **Measured in your own committed output**, `scripts/out/S10_brane_mode_spectrum_sympy_audit.out`:

```
PY_S10_MAIN_D3_Q7_STIFFNESS          ┐
PY_S10_XFORM_FULLGRAD_D3_Q7_STIFFNESS│  all six payloads are
PY_S10_XFORM_DIVONLY_D3_Q7_STIFFNESS ├─ BYTE-IDENTICAL
PY_S10_XFORM_SIGNFLIP_D3_Q7_STIFFNESS│
PY_S10_XFORM_ANISO_D3_Q7_STIFFNESS   │
PY_S10_XCOEF_SCALE_D3_Q7_STIFFNESS   ┘
PY_S10_*_D3_Q7_DIFFERENCE: 0          ← in every one of the six
```

⛔⛔ **Two of those packages exist precisely to change the stiffness FORM.** A control that changes the
form of the stiffness density cannot leave the stiffness density unchanged. ⇒ **Q7 is comparing an
expression against itself, and its zero difference is the control not being applied — not a result.**

⭐ **§4 is the rule this violates: every control must re-enter at the ACTION.** Q7 must read the
**package's own** stiffness density — the one built from that package's action — not a curl expression
reconstructed independently of which package is running.

⭐ **Fix:** make Q7's compared object derive from the same per-package action object the rest of the
package derives from. ⛔ Do not special-case the packages by name, and ⛔ do not hand-write a second
expression per package: if you find yourself typing a stiffness density, you are re-introducing the same
defect one level down. ⭐ The object must arrive from the action.

⚠ **The zero-difference packages are not evidence you can lean on.** Some packages may legitimately
produce a zero difference after the fix and some may not; ⛔ **that is the output, and it is not yours to
predict, arrange, or check against anything.** Compute it and print it.

### ⭐ How to prove to yourself the fix is real — ⛔ do this, do not skip it

Change the stiffness **form** in one package's action — a form change, ⛔ not a coefficient rescale — and
confirm **that package's** `Q7_STIFFNESS` payload moves while another package's does **not**.
⚠ A change that moves all six identically, or none, means you have not fixed it. ⭐ Then revert your
ablation.

---

## ⚠ R2 — one sign test remains undecided under the joint assumptions

`XFORM_ANISO` at `D3` and `D4` emits `undecided_under_joint_assumptions` for a root's `Q3_SIGN`.
⭐ Look at what the joint assumption set actually contains for that package, and at what would have to be
assumed about the package's own scale factor for the sign to be decidable.

⭐ **Emit what you find as tags** — the assumption set in force, and the specific sub-expression whose
sign could not be settled. ⛔ **Do not add an assumption to make it decide.** ⚠ If the spec does not
supply the fact needed, the honest output is the undecided marker **plus** the reason, and that is a
finding I want to see, ⛔ not a thing to fix.

---

## ⭐ What NOT to change — confirmed live by ablation, ⛔ leave alone

`N3` (rank of `M_r` stacked on `kᵀ`), `N7`'s two independent routines, the separate construction of the
two matrix routes and their one-sided corruption behaviour, per-package re-entry at the action for
**everything other than Q7**, Q5's scaling, the dimension tree walk, the real-locus solve, and ⭐ **the
registry allowlist handling** — patching locus validation off before `load_registry` and restoring it
after is **correct**. ⛔ Do not "simplify" it away.

⚠ Your tag **names** are also not to change. They are wired into a comparison configuration.

---

## Report back — ⛔ under 20 lines

1. `R1`: fixed / partially / not, with line numbers, and **what the ablation in the R1 box did** — which
   package's payload moved and which did not. ⛔ Do not report the payload values themselves.
2. `R2`: what you emitted.
3. New tag count, wall-clock, exit code.
4. ⛔ **Do not report what any value came out to be**, and ⛔ do not say whether anything "worked".
5. ⭐ Anything in this directive you believe is **wrong**. ⭐ This is wanted.
