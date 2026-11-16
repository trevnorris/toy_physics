### F3. Exact one-sided mutations of the **copies** (each must FAIL or change a load-bearing number)

**M1 — confinement reaction sign** in copy of `mouth_confinement_projection.py`:

```text
Old: reaction = sp.simplify(work - r**2 * density * vr * mu)
New: reaction = sp.simplify(work + r**2 * density * vr * mu)
```

Expect: `ArithmeticError: Inconsistent angular projection or endpoint work`.

**M2 — crossflow odd divisor** in copy of `held_crossflow_budget.py`:

```text
Old: odd = (pos - neg) / (2 * speed)
New: odd = (pos - neg) / (2 * speed**2)
```

Expect: printed `per_speed_*` change by factor \(1/|U|\); `assert_allclose(neg,-pos)` still passes (mutation is post-antisymmetry). Demonstrates that the reported per-speed normalization is a distinct load-bearing step.

**M3 — transmission translated identity** in copy of `mouth_transmission_work.py`:

```text
Old: if current_jump != 0 or translated != sp.diff(energy, y) or receiving.subs(lam, 1) != 0:
New: if current_jump != 0 or translated != drive or receiving.subs(lam, 1) != 0:
```

Expect: `ArithmeticError: Interface momentum balance failed` (drive retains \(\partial_n\) pieces).

**M4 — unequal force sign** in copy of `unequal_mouth_boundary_work.py`:

```text
Old: force = -np.sum(dv * difference * 6 * axial * base(0)**2 / radius**2)
New: force = +np.sum(dv * difference * 6 * axial * base(0)**2 / radius**2)
```

Expect: `assert_allclose(f, fd, ...)` fails on every mouth.

**M5 — input-field control: corrupt aperture weight in crossflow path**
In a copy used only for a local synthetic check (or patch `signed_readings` after weight build):

```text
Old: weight = np.maximum(1 - x*x - y*y - z*z, 0)**3
New: weight = np.maximum(1 - x*x - y*y - z*z, 0)**2
```

Expect: `assert_allclose(current, 0.25*raw/...)` fails (intake rebuild disagrees with saved current).

**M6 — IBP stress support control** (add to witness or a tiny copy harness after loading one crossflow face):

```python
# After computing slope and n on the mouth plane from a saved packet:
stress_parts = dx**3 * np.sum(slope * dn_dx) / 4
stress_intr = dx**3 * np.sum(-n * d_dx_slope) / 4
# Compact-support equality:
assert abs(stress_parts - stress_intr) < 1e-10 * max(1, abs(stress_intr))
# One-sided: extend slope outside mouth artificially
slope2 = slope + 0.01 * (~mouth)
stress_parts2 = dx**3 * np.sum(slope2 * dn_dx) / 4
stress_intr2 = dx**3 * np.sum(-n * np.gradient(slope2, dx, axis=0)) / 4
assert abs(stress_parts2 - stress_intr2) > 10 * abs(stress_parts - stress_intr)
```

Expect: first assert pass, second assert pass (equality breaks when support leaks).
