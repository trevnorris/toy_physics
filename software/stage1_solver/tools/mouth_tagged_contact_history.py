"""Export and plot the signed readouts of the fine tagged-junction experiment."""
import argparse
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mouth_order_wave_pulse import digest, emit


def main(folder, history, figure):
    outcome = json.loads((folder/'outcome.json').read_text())
    limits = json.loads((folder/'limit-validation.json').read_text())
    if (not limits['verified'] or int(limits['actual']['memory.swap.max']) or
            not outcome['limitsVerified'] or outcome['exitCode'] or outcome['childOutcome']['exitCode']):
        raise ValueError('Unverified input')
    rows = map(json.loads, (folder/'stdout').read_text().splitlines())
    record = next(r for r in rows if r['case'] == 'n48_line128')
    raw = Path(record['raw_path'])
    if digest(raw) != record['raw_sha256']:
        raise ValueError('Raw contact history differs from the recorded run')
    with np.load(raw) as archive:
        columns = archive['columns'].tolist()
        data = archive['samples']
    fields = dict(zip(columns, data.T, strict=True))
    metadata = {key:record[key] for key in ('case','wall_site','ordered_site','lobes','registry',
                                            'offsets','probe_slips','cut','initial_energy',
                                            'events','probes','fields','raw_sha256')}
    metadata.update(columns=columns, live_stiffness='chi/2',
                    samples='Step midpoints in model units. Depth is the mean of step endpoints. '
                            'Forces are on the ordered endpoint or downstream part of the cut. '
                            'Null zero_force_slip denotes two disconnected attachments.')
    with gzip.GzipFile(filename=str(history), mode='wb',mtime=0) as stream:
        stream.write((json.dumps(metadata, allow_nan=False)+'\n').encode())
        for row in data:
            values = [None if np.isnan(value) else float(value) for value in row]
            stream.write((json.dumps(values, separators=(',', ':'), allow_nan=False)+'\n').encode())

    time, step = fields['time'], fields['dt']
    prestrain = max(record['offsets'])
    line = np.cumsum(step*fields['force_on_ordered']*fields['ordered_end_velocity'])
    cut = np.cumsum(step*fields['cut_force']*fields['cut_velocity'])
    np.testing.assert_allclose([line[-1], cut[-1]], [record['probes'][-1][key]
                               for key in ('endpoint_work', 'cut_work')], atol=1e-18, rtol=0)
    fig, (a, b) = plt.subplots(2, 1, figsize=(7.4, 5.2), sharex=True, layout='constrained')
    a.plot(time, fields['zero_force_slip']/prestrain, color='#9e4218', label='Zero-force position')
    a.plot(time, fields['ordered_end_slip']/prestrain, color='#145d72', label='Ordered-end position')
    a.set(ylabel='Sideways position / ε', title='One tagged junction during conversion')
    scale = 100/record['initial_energy']
    b.plot(time+step/2, line*scale, color='#145d72', label='Work at ordered end')
    b.plot(time+step/2, cut*scale, color='#9e4218', ls='--', label='Work across cut 0.25 farther in')
    b.set(xlabel='Time since entry (model units)', ylabel='Work / entry energy (%)')
    for ax in (a, b):
        for event in record['events'][::2]:
            ax.axvline(event['time'], color='0.49', lw=.85, ls=':')
        ax.set_xlim(0, time[-1]+step[-1]/2)
        ax.grid(axis='y', alpha=.18)
        ax.spines[['top','right']].set_visible(False)
        ax.legend(frameon=False)
    b.text(record['events'][0]['time']-.04, .96, 'First loss', rotation=90, ha='right',
           va='top', transform=b.get_xaxis_transform(), fontsize=8)
    b.text(record['events'][2]['time']+.04, .05, 'Breakup', rotation=90, ha='left',
           va='bottom', transform=b.get_xaxis_transform(), fontsize=8)
    fig.savefig(figure, dpi=160)
    plt.close(fig)
    emit(rows=len(data), history=str(history), history_sha256=digest(history),
         figure=str(figure), figure_sha256=digest(figure))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ('folder', 'history', 'figure'):
        parser.add_argument(field, type=Path)
    args = parser.parse_args()
    main(args.folder, args.history, args.figure)
