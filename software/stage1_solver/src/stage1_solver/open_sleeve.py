"""Steady neutral Stage-1 carrier in an open sleeve joined to a finite receiver.

The feed phase adjusts to a specified inward number rate. The far receiver
fixes the phase reference. Both baths fix the carrier amplitude. The sleeve
and receiver walls reflect the carrier. The solve has no body mass input.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import spsolve

from .config import TensorGridSpec


@dataclass
class OpenSleeve:
    grid: TensorGridSpec
    mouth_cells: int
    sleeve_cells: int
    hbar: float = 1.0
    constituent_mass: float = 1.0
    eos_k: float = 0.5
    bath_density: float = 1.0

    def __post_init__(self):
        g = self.grid
        if not (0 < self.mouth_cells < g.nr and 0 < self.sleeve_cells < g.nw):
            raise ValueError("Mouth edge and sleeve exit must be interior grid faces")
        if not (g.r_min == 0 and g.w_min == 0 and g.r_max > 0 and g.w_max > 0):
            raise ValueError("The domain must start at the axis and entrance")
        if not all(np.isfinite(x) and x > 0 for x in
                   (g.r_max, g.w_max, self.hbar, self.constituent_mass, self.eos_k, self.bath_density)):
            raise ValueError("Physical scales must be positive and finite")
        self.dr = g.r_max / g.nr
        self.dw = g.w_max / g.nw
        radii = np.arange(g.nr + 1) * self.dr
        shells = 4 * np.pi / 3 * np.diff(radii**3)
        # Only occupied cells enter the graph. Missing links implement no
        # normal derivative on the sleeve, shoulder and outer radial wall.
        r, w = np.indices((g.nr, g.nw))
        mask = (r < self.mouth_cells) | (w >= self.sleeve_cells)
        ids = np.full(mask.shape, -1, dtype=int)
        ids[mask] = np.arange(mask.sum())
        self.r_index, self.w_index = r[mask], w[mask]
        self.volume = shells[self.r_index] * self.dw
        rr, ww = np.nonzero(mask[:-1, :] & mask[1:, :])
        radial_links = 4 * np.pi * radii[rr + 1]**2 * self.dw / self.dr
        rz, wz = np.nonzero(mask[:, :-1] & mask[:, 1:])
        self.left = np.r_[ids[rr, ww], ids[rz, wz]]
        self.right = np.r_[ids[rr + 1, ww], ids[rz, wz + 1]]
        self.conductance = np.r_[radial_links, shells[rz] / self.dw]
        self.inlet = np.where(self.w_index == 0, 2 * shells[self.r_index] / self.dw, 0.0)
        self.outlet = np.where(self.w_index == g.nw - 1, 2 * shells[self.r_index] / self.dw, 0.0)
        size = len(self.volume)
        off_diag = coo_matrix((-self.conductance, (self.left, self.right)), shape=(size, size)).tocsr()
        off_diag = off_diag + off_diag.T
        diagonal = -np.asarray(off_diag.sum(axis=1)).ravel() + self.inlet + self.outlet
        self.kinetic = self.hbar**2 / (2 * self.constituent_mass)
        self.operator = self.kinetic * (off_diag + diags(diagonal))

    def solve(self, rate):
        if not np.isfinite(rate) or rate < 0:
            raise ValueError("Specify a finite inward rate or zero for the control")
        bath = np.sqrt(self.bath_density)
        h_coeff = 5 * self.eos_k / 4
        mu = h_coeff * self.bath_density**4
        scale = self.hbar / self.constituent_mass
        x = np.full_like(self.volume, bath)
        y = np.zeros_like(x)
        phase = 0.0

        def evaluate(x, y, phase):
            psi = x + 1j * y
            feed = bath * np.exp(1j * phase)
            n = x*x + y*y
            local = self.volume * (h_coeff * n**4 - mu)
            forcing = self.kinetic * (self.inlet * feed + self.outlet * bath)
            field = self.operator @ psi + local * psi - forcing
            inlet_rate = scale * np.dot(self.inlet, (feed.conjugate() * psi).imag)
            residual = np.r_[field.real, field.imag, inlet_rate - rate]
            return residual, local, feed

        for step in range(25):
            residual, local, feed = evaluate(x, y, phase)
            error = np.max(np.abs(residual))
            if error < 2e-13:
                return x + 1j * y, phase, error
            slope = self.volume * 8 * h_coeff * (x*x + y*y)**3
            base = self.operator + diags(local)
            xy = diags(slope*x*y)
            phase_column = -1j * self.kinetic * self.inlet * feed
            flux_phase = -scale * np.dot(self.inlet, feed.real*x + feed.imag*y)
            jac = bmat([
                [base + diags(slope*x*x), xy, csr_matrix(phase_column.real[:, None])],
                [xy, base + diags(slope*y*y), csr_matrix(phase_column.imag[:, None])],
                [csr_matrix(-scale*feed.imag*self.inlet), csr_matrix(scale*feed.real*self.inlet),
                 csr_matrix([[flux_phase]])],
            ], format="csc")
            delta = spsolve(jac, -residual)
            dx, dy = np.split(delta[:-1], 2)
            for power in range(16):
                factor = 0.5**power
                trial = x + factor*dx, y + factor*dy, phase + factor*delta[-1]
                if np.max(np.abs(evaluate(*trial)[0])) < error:
                    x, y, phase = trial
                    break
            else:
                raise ArithmeticError(f"Newton step {step} did not reduce residual {error:g}")
        raise ArithmeticError("Stationary carrier failed to converge")

    def budget(self, psi, phase):
        # Read currents directly on each link from the solved complex field.
        scale = self.hbar / self.constituent_mass
        bath = np.sqrt(self.bath_density)
        inlet = scale * self.inlet * (bath*np.exp(-1j*phase)*psi).imag
        exterior = scale * self.outlet * (psi.conj()*bath).imag
        link = scale * self.conductance * (psi[self.left].conj()*psi[self.right]).imag
        imbalance = inlet - exterior
        np.add.at(imbalance, self.left, -link)
        np.add.at(imbalance, self.right, link)
        sleeve = self.w_index < self.sleeve_cells
        into_receiver = link[sleeve[self.left] & ~sleeve[self.right]].sum()
        return (inlet.sum(), into_receiver, exterior.sum(),
                np.max(np.abs(imbalance / self.volume)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rate", type=float, required=True, help="inward number per unit time")
    parser.add_argument("--output", type=Path, required=True, help="new directory for solved fields")
    parser.add_argument("--levels", type=int, nargs="+", default=[8, 16, 32],
                        help="cells per unit mouth radius and per sleeve length")
    args = parser.parse_args()
    if any(level < 2 for level in args.levels) or len(args.levels) != len(set(args.levels)):
        parser.error("each resolution must be distinct and at least two")
    if not np.isfinite(args.rate) or args.rate < 0:
        parser.error("rate must be finite and nonnegative")
    args.output.mkdir(parents=True, exist_ok=False)
    for level in args.levels:
        grid = TensorGridSpec(r_max=2.0, nr=2*level, w_min=0.0, w_max=3.7, nw=2*level)
        solver = OpenSleeve(grid, mouth_cells=level, sleeve_cells=level)
        psi, phase, residual = solver.solve(args.rate)
        entry, transfer, out, local_error = solver.budget(psi, phase)
        print(f"level={level} prescribed={args.rate:.12g} inlet={entry:.12g} "
              f"into_receiver={transfer:.12g} far_exit={out:.12g} feed_phase={phase:.12g} "
              f"min_density={min(abs(psi)**2):.12g} max_cell_rate_per_volume={local_error:.3e} "
              f"max_integrated_equation_residual={residual:.3e}", flush=True)
        rates = np.array([entry, transfer, out])
        if not np.all(np.isfinite(rates)) or np.max(abs(rates - args.rate)) > 1e-10 or local_error > 1e-8:
            raise ArithmeticError("Prescribed intake or internal flux balance failed")
        np.savez_compressed(args.output / f"level_{level}.npz", psi=psi, feed_phase=phase,
                            r_index=solver.r_index, w_index=solver.w_index, dr=solver.dr, dw=solver.dw,
                            mouth_cells=solver.mouth_cells, sleeve_cells=solver.sleeve_cells,
                            bath_density=solver.bath_density, hbar=solver.hbar,
                            constituent_mass=solver.constituent_mass, eos_k=solver.eos_k, rate=args.rate)


if __name__ == "__main__":
    main()
