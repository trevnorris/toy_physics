#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional


PI = math.pi


@dataclass
class CoherentBranchState:
    chi0: float
    delta_U: float
    Z_W: float
    epsilon_W: float
    epsilon_eta: float
    Lambda: float
    zeta: Optional[float] = None
    chi_Q: Optional[float] = None
    N_Q: Optional[float] = None
    Pi_tr: Optional[float] = None
    mhat0: float = 1.0


@dataclass
class WeakAxisymmetricDrifts:
    dln_chi0: float
    dln_delta_U: float
    dln_Z_W: float
    dln_epsilon_W: float
    dln_epsilon_eta: float
    dln_Lambda: float


@dataclass
class EvaluationResult:
    epsilon_split: float
    R_tr: float
    R_target: float
    M_mix: float
    M_tr: Optional[float]
    C_mix: float
    rho_alpha: Optional[float]
    support_regime: Optional[str]
    dln_R_tr: float
    dln_R_target: float
    dln_epsilon_eta: float
    tracking_codrift: float
    outgoing_residual_NQ: Optional[float]
    outgoing_residual_chiQ: Optional[float]
    source_map_consistency_residual: Optional[float]


def support_enhancement(zeta: float, epsilon_split: float) -> float:
    denom = 1.0 - zeta * epsilon_split
    if abs(denom) < 1e-15:
        raise ZeroDivisionError('support enhancement denominator is too small')
    return 1.0 + zeta * (1.0 - epsilon_split) / denom


def epsilon_split(state: CoherentBranchState) -> float:
    return state.epsilon_W * (1.0 - (2.0 / 11.0) * state.delta_U / (1.0 + state.delta_U))


def compute_R_tr(state: CoherentBranchState) -> float:
    return (1.0 + state.chi0 / (1.0 + state.delta_U)) / (1.0 + state.chi0)


def compute_R_target(state: CoherentBranchState, eps: float) -> float:
    denom = state.Z_W * (1.0 + state.chi0) ** 2
    return state.Lambda * (1.0 - state.epsilon_eta) * (1.0 - eps) ** 2 / denom


def compute_M_mix(state: CoherentBranchState, eps: float) -> float:
    denom = PI ** 2 * (1.0 - state.epsilon_eta) * (1.0 - eps)
    return 8.0 * state.Z_W * (1.0 + state.chi0) ** 2 / denom


def compute_C_mix(state: CoherentBranchState, eps: float) -> float:
    return 8.0 * state.Lambda * (1.0 - eps) / (PI ** 2)


def classify_support_regime(rho_alpha: float) -> str:
    if rho_alpha <= 1.0:
        return 'mixed-only'
    if rho_alpha <= 2.0:
        return 'lowest-symmetric-twin'
    return 'non-twin'


def evaluate(state: CoherentBranchState, drifts: WeakAxisymmetricDrifts) -> EvaluationResult:
    eps = epsilon_split(state)
    r_tr = compute_R_tr(state)
    r_target = compute_R_target(state, eps)
    m_mix = compute_M_mix(state, eps)
    c_mix = compute_C_mix(state, eps)

    m_tr = None
    if state.zeta is not None:
        m_tr = m_mix * support_enhancement(state.zeta, eps)

    rho_alpha = None
    support_regime = None
    if state.Pi_tr is not None:
        rho_alpha = state.Pi_tr / c_mix
        support_regime = classify_support_regime(rho_alpha)

    dln_eps = drifts.dln_epsilon_W - (2.0 * state.delta_U / ((1.0 + state.delta_U) * (11.0 + 9.0 * state.delta_U))) * drifts.dln_delta_U
    dln_r_tr = -(
        state.chi0 * state.delta_U
        / ((1.0 + state.chi0) * (1.0 + state.delta_U) * (1.0 + state.chi0 + state.delta_U))
    ) * ((1.0 + state.delta_U) * drifts.dln_chi0 + (1.0 + state.chi0) * drifts.dln_delta_U)
    dln_r_target = (
        drifts.dln_Lambda
        - drifts.dln_Z_W
        - (state.epsilon_eta / (1.0 - state.epsilon_eta)) * drifts.dln_epsilon_eta
        - (2.0 * state.chi0 / (1.0 + state.chi0)) * drifts.dln_chi0
        - (2.0 * eps / (1.0 - eps)) * dln_eps
    )
    tracking_codrift = (1.0 + state.delta_U) * drifts.dln_chi0 + (1.0 + state.chi0) * drifts.dln_delta_U

    outgoing_residual_nq = None
    if state.N_Q is not None:
        outgoing_residual_nq = state.N_Q - 1.0

    outgoing_residual_chiq = None
    if state.chi_Q is not None:
        outgoing_residual_chiq = state.chi_Q - 1.0

    source_map_consistency = None
    if state.N_Q is not None and state.chi_Q is not None:
        source_map_consistency = state.N_Q - 1.0 / state.chi_Q

    return EvaluationResult(
        epsilon_split=eps,
        R_tr=r_tr,
        R_target=r_target,
        M_mix=m_mix,
        M_tr=m_tr,
        C_mix=c_mix,
        rho_alpha=rho_alpha,
        support_regime=support_regime,
        dln_R_tr=dln_r_tr,
        dln_R_target=dln_r_target,
        dln_epsilon_eta=drifts.dln_epsilon_eta,
        tracking_codrift=tracking_codrift,
        outgoing_residual_NQ=outgoing_residual_nq,
        outgoing_residual_chiQ=outgoing_residual_chiq,
        source_map_consistency_residual=source_map_consistency,
    )


def template_payload() -> dict:
    return {
        'state': asdict(CoherentBranchState(
            chi0=0.0,
            delta_U=0.0,
            Z_W=0.0,
            epsilon_W=0.0,
            epsilon_eta=0.0,
            Lambda=0.0,
            zeta=0.0,
            chi_Q=1.0,
            N_Q=1.0,
            Pi_tr=1.0,
            mhat0=1.0,
        )),
        'drifts': asdict(WeakAxisymmetricDrifts(
            dln_chi0=0.0,
            dln_delta_U=0.0,
            dln_Z_W=0.0,
            dln_epsilon_W=0.0,
            dln_epsilon_eta=0.0,
            dln_Lambda=0.0,
        )),
    }


def load_payload(path: Path) -> tuple[CoherentBranchState, WeakAxisymmetricDrifts]:
    data = json.loads(path.read_text())
    state = CoherentBranchState(**data['state'])
    drifts = WeakAxisymmetricDrifts(**data['drifts'])
    return state, drifts


def main() -> None:
    parser = argparse.ArgumentParser(description='Evaluate the reduced 5PN finish packet on an actual branch dataset.')
    parser.add_argument('--input', type=Path, help='JSON file containing "state" and "drifts" blocks')
    parser.add_argument('--template-json', type=Path, help='Write a template JSON payload and exit')
    parser.add_argument('--demo', action='store_true', help='Run an internal demo payload')
    args = parser.parse_args()

    if args.template_json is not None:
        args.template_json.write_text(json.dumps(template_payload(), indent=2))
        print(f'wrote template JSON to {args.template_json}')
        return

    if args.demo:
        demo_state = CoherentBranchState(
            chi0=1.5,
            delta_U=2.0/3.0,
            Z_W=0.2,
            epsilon_W=0.05,
            epsilon_eta=0.1,
            Lambda=0.3,
            zeta=1.0,
            chi_Q=1.0,
            N_Q=1.0,
            Pi_tr=1.3333333333333333,
            mhat0=1.0,
        )
        demo_drifts = WeakAxisymmetricDrifts(
            dln_chi0=0.0,
            dln_delta_U=0.0,
            dln_Z_W=0.0,
            dln_epsilon_W=0.0,
            dln_epsilon_eta=0.0,
            dln_Lambda=0.0,
        )
        result = evaluate(demo_state, demo_drifts)
        print(json.dumps(asdict(result), indent=2))
        return

    if args.input is None:
        parser.error('provide --input, --template-json, or --demo')

    state, drifts = load_payload(args.input)
    result = evaluate(state, drifts)
    print(json.dumps(asdict(result), indent=2))


if __name__ == '__main__':
    main()
