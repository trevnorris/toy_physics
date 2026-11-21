"""Check the endpoint response by integrating the nonlinear signed shear drives."""
import argparse

import numpy as np
from scipy.integrate import solve_ivp

from neutral_patch_sites import emit
from neutral_patch_structural_endpoints import Cell, preparation


def main(side, seed):
    cell = Cell(side,6.)
    rest, _, _, _ = cell.settle(preparation(cell,seed))
    times = [.5, 1., 1.5, 2., 4.]

    def response(amplitude, pulse):
        state = rest.ravel()
        readings = []
        for begin,end,strain,outputs in [(0.,1.,amplitude,[.5,1.]),
                (1.,4.,0. if pulse else amplitude,[1.5,2.,4.])]:
            # Convert the unsheared cell-coordinate force to a lab force and
            # its velocity back to cell coordinates for Cartesian unit mobility.
            mobility = np.array([[1+strain**2,-strain],[-strain,1.]])
            def drift(t, flat):
                _,g,_ = cell.read(flat.reshape(-1,2),strain)
                return (-g@mobility).ravel()
            ivp = solve_ivp(drift,(begin,end),state,t_eval=outputs,method='DOP853',
                            rtol=1e-11,atol=5e-13)
            if not ivp.success:
                raise ArithmeticError(ivp.message)
            for t,y in zip(outputs,ivp.y.T,strict=True):
                reading_strain = 0. if pulse and t==1. else strain
                readings.append(cell.read(y.reshape(-1,2),reading_strain)[2])
            state = ivp.y[:,-1]
        return np.array(readings)

    # Both trials have identical nonlinear trajectories through the first step
    # interval; the step is subsequently removed only in the pulse preparation.
    for pulse in [False,True]:
        amplitude = .0001
        plus = response(amplitude,pulse)
        minus = response(-amplitude,pulse)
        emit(test='nonlinear_endpoint_drive',side=side,seed=seed,pulse=pulse,amplitude=amplitude,
             times=times,positive=plus.tolist(),negative=minus.tolist(),
             incremental_stress_per_strain=((plus-minus)/(2*amplitude)).tolist())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('side',type=int,choices=[4,6])
    parser.add_argument('seed',type=int)
    args=parser.parse_args()
    main(args.side,args.seed)
