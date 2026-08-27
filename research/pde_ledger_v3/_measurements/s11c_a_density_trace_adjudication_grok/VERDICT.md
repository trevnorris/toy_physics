# Adjudication — S11c-a T-e shifted density trace (LAB_HELD, FACE_MINUS, DELTA_W)

## Commands

```bash
python3 research/pde_ledger_v3/_measurements/s11c_a_density_trace_adjudication_grok/derive_density_bg_normal_deriv.py
python3 research/pde_ledger_v3/_measurements/s11c_a_density_trace_adjudication_grok/compare_wl_extra_to_forced_zero.py
```

Literal stdout: `derive_density_bg_normal_deriv.stdout`, `compare_wl_extra_to_forced_zero.stdout`.

## Spec quotes (governing)

- §2b (lines 214–222): both density representatives are maps of the in-plane anchor through `W_bg(y)`.
- §2d (lines 250–251): `𝔅⁰` contains `ρ_4D,bg⁰`.
- §3c (lines 375–392):
  - law `δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰(x,h_s⁰)`;
  - "Every background face value or normal derivative appearing in this law is obtained by differentiating a member of the supplied background state `𝔅⁰` (§2d); none may be introduced as a free premise";
  - "the supplied density background depends on the in-plane anchor, not on `w`".

## Verdict

(a) Correct first-order §3c shifted density trace at the disputed key:
`delta_rho_4D_face_minus - (W_0/2) eta_bg w1_profile delta_rho_4D_face_minus_dw`

(b) `∂_w ρ_4D,bg⁰ = 0` for both §2b members of `𝔅⁰` (script stdout).

(c) PY is §3c-faithful. WL introduces undefined `rhoBulkBackground(x,w,t)` (no `:=` in the WL source) and carries live `∂_w`/`∂²_w` of it, violating the free-premise ban.

(d) Difference is a spec-forced zero, not a genuine missing physical term.

(e) No ambiguity: §3c states the density background is not `w`-dependent and requires grounding in `𝔅⁰`.
