"""Plot the settled midplanes of the two original period-24 releases."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from neutral_patch_sites import emit


def profile(path):
    with np.load(path) as state:
        pos, flat, box = state['position'], state['flat'], state['box']
    lower, upper = np.split(pos[:, 2], 2)
    mid = (lower+upper)/2
    mid -= mid.mean()
    half_columns = 2*round(np.sqrt(mid.size))
    index = np.rint(flat[:len(mid), 0]/box[0]*half_columns).astype(int) % half_columns
    height = np.bincount(index, mid)/np.bincount(index)
    spread = np.max(abs(mid-height[index]))
    height = np.roll(height, half_columns//2-np.argmax(height))
    return np.linspace(-.5, .5, half_columns+1), np.r_[height, height[0]], spread


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('first', type=Path)
    parser.add_argument('second', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    fig, ax = plt.subplots(figsize=(7.1, 3.7), layout='constrained')
    ax.axhline(0, linewidth=.8, color='#8a8c8c')
    for path, label, color in [(args.first, r'$f$', '#005987'),
                                (args.second, r'$f^2$', '#ab4d0b')]:
        x, y, spread = profile(path)
        ax.plot(x, y, color=color, linewidth=2.4, label=label, marker='.', markersize=3)
        emit(state=str(path), peak_to_trough=np.ptp(y), max_deviation_from_column_mean=spread)
    ax.set(xlim=(-.5,.5), xlabel='Position along x relative to the crest (periods)',
           ylabel='Midplane height relative to the mean\n(initial layer-gap units)')
    ax.set_title('The settled wrinkles at the same projected area', loc='left', weight='bold')
    ax.spines[['right','top']].set_visible(False)
    ax.legend(title='Interaction weight', frameon=False, loc='upper right')
    fig.savefig(args.output, dpi=180)
