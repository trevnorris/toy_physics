"""Static electric-source readout of the neutral carrier under a common test charge.

Reflect the whole inlet/receiver preparation in w. The prescribed positive
number rate is the same in both local frames. No electromagnetic backreaction.
"""

from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import breadth_first_order
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from stage1_solver.config import TensorGridSpec
from stage1_solver.open_sleeve import OpenSleeve


def transit_readout(solver, psi, phase, weights):
    """Expected weighted residence to absorption using the actual directed link fluxes."""
    inventory = solver.volume * abs(psi)**2
    # hbar/m=1, bath amplitude=1 in the experiments below.
    incoming = solver.inlet * (np.exp(-1j * phase) * psi).imag
    escaping = -solver.outlet * psi.imag
    signed = solver.conductance * (psi[solver.left].conj() * psi[solver.right]).imag
    start = np.where(signed > 0, solver.left, solver.right)
    end = np.where(signed > 0, solver.right, solver.left)
    links = coo_matrix((abs(signed), (start, end)), shape=(psi.size, psi.size)).tocsr()
    links.eliminate_zeros()
    if np.any(incoming < 0) or np.any(escaping < 0):
        raise ArithmeticError("Boundary backflow would need separate entry and exit tags")
    # Attach one common reservoir vertex. Every measured cell must both receive
    # tags from the feed and have a directed path out, including near the wall.
    origin = psi.size
    rows, cols = links.nonzero()
    feed_cells = np.flatnonzero(incoming)
    exit_cells = np.flatnonzero(escaping)
    rows = np.r_[rows, np.full(feed_cells.size, origin), exit_cells]
    cols = np.r_[cols, feed_cells, np.full(exit_cells.size, origin)]
    graph = coo_matrix((np.ones(rows.size), (rows, cols)), shape=(origin+1, origin+1)).tocsr()
    for direction in (graph, graph.T.tocsr()):
        reached = breadth_first_order(direction, origin, directed=True, return_predecessors=False)
        if reached.size != origin + 1:
            raise ArithmeticError("Unresolved stagnant or disconnected carrier cells")

    outgoing = np.asarray(links.sum(axis=1)).ravel() + escaping
    # Backward residence equation: outgoing_i T_i - sum_j F_ij T_j = N_i f_i.
    times = spsolve(diags(outgoing, format="csr") - links, inventory[:, None] * weights)
    return incoming @ times


def run(radius, cells, rate=0.25):
    grid = TensorGridSpec(r_max=3.0, nr=3*cells, w_min=0.0, w_max=3.0, nw=3*cells)
    solver = OpenSleeve(grid, mouth_cells=round(radius*cells), sleeve_cells=cells)
    psi, phase, residual = solver.solve(rate)
    entry, into_receiver, exit_rate, continuity = solver.budget(psi, phase)
    # The common test charge and Gaussian localization length are one. Fix
    # the laboratory origin at the inlet. Orientation is reversed by reflection.
    w = (solver.w_index + 0.5) * solver.dw
    weights = np.column_stack((np.ones_like(w), 2*w, -2*w))
    inventory = solver.volume * abs(psi)**2
    raw = inventory @ weights
    transit = transit_readout(solver, psi, phase, weights)
    error = raw - transit
    excess = (solver.volume * (abs(psi)**2 - solver.bath_density)) @ weights
    print(f"R={radius:g} cells_per_length={cells} inlet={entry:.12g} cut={into_receiver:.12g} "
          f"exit={exit_rate:.12g} GNLS_residual={residual:.3e} local_continuity={continuity:.3e}")
    print(f"  flat_carrier={raw[0]:.12g} flat_transit={transit[0]:.12g} "
          f"mean_residence={transit[0]/entry:.12g} flat_residual={error[0]:+.5e}")
    print(f"  excess_over_rest_flat={excess[0]:+.12g} "
          f"excess_over_rest_depth_plus={excess[1]:+.12g} excess_over_rest_depth_minus={excess[2]:+.12g}")
    print(f"  depth_residual_plus={error[1]:+.5e} depth_residual_minus={error[2]:+.5e}", flush=True)
    if not np.allclose([entry, into_receiver, exit_rate], rate, atol=1e-10, rtol=0):
        raise ArithmeticError("Unmatched intake or outlet")
    if continuity > 1e-8 or not np.allclose(transit, raw, atol=1e-9, rtol=1e-10):
        raise ArithmeticError("Continuity or transported inventory failed")


if __name__ == "__main__":
    for resolution in (8, 16):
        for mouth_radius in (1.0, 2.0):
            run(mouth_radius, resolution)
