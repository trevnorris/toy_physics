"""Export the three archived pressure-prepared receivers into the inertial probe format."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


def load(path):
    return json.loads(path.read_bytes())


def export(out):
    batch = Path('_scratch/s11c/thermal87_batch02')
    manifest = load(Path('research/pde_ledger_v3/_measurements/mouth_bulk_temperature_steps.json'))
    out.mkdir(exist_ok=True,parents=True)
    for history in range(3):
        setup = {}
        for phase in ('prepare','undriven'):
            name = f'r4h{history}_{phase}'
            folder = batch/name
            raw = (folder/'stdout').read_bytes()
            receipt = load(folder/'outcome.json')
            assert receipt['childOutcome']==load(folder/'child-outcome.json')
            assert receipt['childOutcome']['guardReason'] is None and receipt['exitCode']==receipt['childOutcome']['exitCode']==0
            assert load(folder/'limit-validation.json')['verified'] and load(folder/'duration-validation.json')['verified']
            assert sha256(raw).hexdigest()==manifest['runs'][name]['stdout_sha256']
            setup[phase], *frames = map(json.loads,raw.splitlines())
            assert frames[-1]['frame']==48 and setup[phase]['temperature']==1.024
        path = batch/f'r4h{history}_prepare/native/end.npz'
        digest = sha256(path.read_bytes()).hexdigest()
        assert setup['undriven']['source_sha256']==digest
        with np.load(path,allow_pickle=False) as state:
            x, box, axis = (state[k] for k in ('pos','box','axis'))
            assert int(state['history'])==history and int(state['rung'])==4 and str(state['phase'])=='prepare'
            np.testing.assert_allclose(axis,setup['prepare']['axis'])
            np.testing.assert_allclose(axis,axis.T,atol=1.e-14)
            np.testing.assert_allclose([np.trace(axis),np.sum(axis*axis)],[0.,.5],atol=1.e-12)
            output = out/f'history{history}.txt'
            with output.open('x') as stream:
                print(len(x),history,setup['prepare']['temperature'],setup['prepare']['volume_mobility'],setup['prepare']['cutoff'],digest,file=stream)
                for array in (box, np.linalg.inv(box), axis, x):
                    np.savetxt(stream,array,fmt='%.17g')
        print(json.dumps(dict(history=history,archive=str(path),source_sha256=digest,file=str(output),export_sha256=sha256(output.read_bytes()).hexdigest())))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    export(parser.parse_args().output)
