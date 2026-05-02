# Projected Maxwell Push — Bundle Bridge Master Notes

This bundle contains two linked continuations of the projected-Maxwell / near-throat work.

## 1. Grouped `P2` / full-bundle bridge

File:
- `step_09_projected_maxwell_p2_bridge_notes.md`

This note shows how near-throat projection-first Maxwell corrections feed the moving-throat bundle objects

\[
D_0,\ D_2,\ D_4,\ N_0,\ N_2,\ N_4
\]

and therefore the derived quantities

\[
u_2,\ u_4,\ P_0,\ P_2,\ P_4,\ \Xi_1=\frac{P_1}{P_0}.
\]

Main takeaways:

- mouth-local projection can move the bundle at first order \(O(\ell)\),
- the isotropic compatibility surface is first sensitive to \(z_2,z_4,n_0\),
- \(z_0\) shifts \(P_0\) directly but cancels out of the fixed-target Stage-18 compatibility relation,
- weak-axisymmetric projected-Maxwell slippage feeds \(\Xi_1\) with the correct grouped signature.

In the updated SymPy audit, the \(z_0\) cancellation is now shown in two
equivalent ways:

- by solving the perturbed one-pole and normalization equations for the shifted
  wall stiffness \(K\) and subtracting those two surfaces,
- and by matching that difference to the directly eliminated compatibility
  surface
  \[
  \frac{N_0+\varepsilon n_0}{P_{0,\mathrm{target}}}
  -
  \frac{3(S+\varepsilon z_2)^2}{T+\varepsilon z_4},
  \]
  which contains no \(z_0\) slot at all.

So the disappearance of \(z_0\) is now an exact elimination result, not a
parallel-typing artifact.

The script now also checks the stronger ratio-form target
\[
P_{0,\mathrm{target}}(\varepsilon)=\frac{N_0+\varepsilon n_0}{D_{0,\mathrm{target}}},
\]
with an explicit transported denominator slot \(D_{0,\mathrm{target}}\). In
that form the normalization \(K\)-surface becomes
\[
K=B_0+Z_0+\varepsilon z_0+D_{0,\mathrm{target}},
\]
so both \(z_0\) and \(n_0\) cancel from the resulting compatibility transport,
leaving only the geometric pair \(z_2,z_4\).

## 2. Primitive one-port bridge

File:
- `step_10_projected_maxwell_stage4_primitive_bridge_notes.md`

This note goes one level more microscopic and starts from the Stage-4/Stage-5 one-port primitive data

\[
Q,\ S_2,\ H,\ \Delta,\ P,\ G_W
\]

with mouth-local slippages

\[
q_1,\ s_1,\ h_1,\ d_1,\ p_1,\ g_1.
\]

It computes the induced

\[
z_0,\ z_2,\ z_4,\ n_0,\ n_2,\ n_4
\]

exactly.

Main takeaways:

- \(z_0\) depends only on \(q_1,d_1\),
- \(n_0\) depends only on \(p_1,d_1\),
- the static prefactor slope \(\Xi_1\) is first controlled by \((q_1,d_1,p_1)\),
- \(s_1,h_1\) first enter through the one-pole / compatibility transport,
- \(g_1\) first enters through \(n_2,n_4\), i.e. the constant-prefactor transport.

## 3. Scripts

General bridge:
- `step_08_projected_maxwell_push_bundle_master_sympy.py`
- `step_09_projected_maxwell_p2_bridge_sympy.py`

Primitive one-port bridge:
- `step_10_projected_maxwell_stage4_primitive_bridge_sympy.py`

## 4. Practical next move

The clean next derivation step is to choose one concrete throat-local projected ansatz that produces explicit primitive slippages

\[
(q_1,\ s_1,\ h_1,\ d_1,\ p_1,\ g_1),
\]

preferably in the structurally clean \(H=Z\) gauge-localization scheme, and test whether it moves

\[
\delta\mathcal C,\ \delta P_2,\ \delta P_4,\ \Xi_1
\]

in the right direction on the actual moving-throat branch.
