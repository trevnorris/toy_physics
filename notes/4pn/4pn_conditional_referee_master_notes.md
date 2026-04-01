# 4PN conditional referee master — result ledger

## What this step accomplished

This step freezes the sharpest **full conservative 4PN theorem statement** currently available inside the declared closure hierarchy.

It does not try to solve the full moving-throat PDE. Instead it replays the three decisive already-solved stages and then proves the interface theorem tying them together:

1. the **2.5PN master audit** already reduces the universal dissipative bridge to the canonically normalized STF quadrupole branch, with one narrow normalization gap left open;
2. the **local 4PN referee master** already closes the entire local instantaneous 4PN sector exactly;
3. the **4PN hereditary bridge audit** shows that the hereditary coefficient is not a new free datum, but is fixed by the same quadrupole normalization.

So the full conservative 4PN problem is now conditionally reduced to the **same** passive/outgoing quadrupole normalization problem already isolated on the 2.5PN side.

---

## Stage-audit status now frozen

The master replay confirms:

- `2_5pn_master_session_sympy_audit.py` passes,
- `4pn_local_referee_master_sympy_audit.py` passes,
- `4pn_tail_hereditary_bridge_audit.py` passes.

That means the following three legs are simultaneously stable:

- the narrowed **quadrupole normalization gap** from 2.5PN,
- the exact **local 4PN** generic-frame ordinary representative,
- the exact **4PN hereditary bridge coefficient relation**.

---

## Exact interface theorem

Let the canonically normalized STF quadrupole odd coefficient be

\[
\gamma_{\rm quad}^{\rm eff}.
\]

Then the unique compatible conservative 4PN hereditary coefficient is

\[
\boxed{
C_{\rm tail}=rac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
}
\]

So if the 2.5PN normalization target holds,

\[
\gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
\]

then automatically

\[
\boxed{
C_{\rm tail}=\frac{G^2M}{5c^8}.
}
\]

This is the exact GR conservative 4PN tail coefficient.

---

## Equivalent branch forms

The same result may be written in either of the two branch languages already isolated in the 2.5PN program.

### Geometric normalization form

\[
\gamma_{\rm quad}^{\rm eff}=\mathcal N_Q\frac{a^5}{27c_s^5},
\qquad
\mathcal N_Q^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
\]

### Canonical invariant-pair form

\[
\gamma_{\rm quad}^{\rm eff}
=
9\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]

On the Burke–Thorne target branch,

\[
\overline K_0^{\rm target}=\frac{64G}{45c^5}\Omega_Q^5,
\qquad
\overline K_2^{\rm target}=\frac{16G}{45c^5}\Omega_Q^3,
\]

and substitution gives

\[
\gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
\qquad
C_{\rm tail}=\frac{G^2M}{5c^8}.
\]

So **no new independent normalization datum opens at 4PN**.

---

## Strongest honest theorem statement available now

Within the declared closure hierarchy:

\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}},
\]

with

- `L_4PN^local` already assembled exactly by the local referee chain,
- `L_4PN^tail` fixed by the same STF quadrupole normalization that controls the 2.5PN Burke–Thorne channel.

Therefore the remaining gap between the present hierarchy and a fully unconditional conservative 4PN theorem is **not** a new 4PN-specific normalization constant. It is exactly the same narrow passive/outgoing quadrupole normalization problem already isolated by the 2.5PN program.

---

## Best next move after this step

The next sharp target is no longer local 4PN algebra. That part is done.

The next real theorem gate is the one already exposed by the 2.5PN notes:

- derive the passive/outgoing quadrupole normalization from the moving-throat side,
- or derive enough of the higher conservative grouped real `P2` bundle to fix the canonical pair `(\overline K_0,\overline K_2)` on the natural branch.

Once that is closed, both

- the normalized 2.5PN odd quadrupole coefficient, and
- the normalized 4PN conservative hereditary coefficient

close simultaneously.
