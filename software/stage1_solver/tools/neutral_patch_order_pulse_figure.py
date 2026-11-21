"""Show the strain and directly recorded gauge tractions on the free-backed run."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import digest
from neutral_patch_conversion import clock


def main():
    source = Path('_scratch/s11c/pulse_probe01/B_central_fine.npz')
    archive = Path('research/pde_ledger_v3/_measurements/neutral_patch_shear_pulse.json')
    if digest(source) != json.loads(archive.read_text())['gauge_file_sha256'][str(source)]:
        raise ValueError('Changed free-backed gauge data')
    t, _, order, _, gates = clock(WallTrace(48, .5), 2)
    with np.load(source) as data:
        np.testing.assert_array_equal(data['t'], t)
        i = data['probe_cases'].tolist().index('pulse')
        tau = data['tau'][:, i]
        speed = data['u'][:, i]
        gap = data['z'][0] - data['z'][1]
    strain = cumulative_trapezoid((speed[:, 0]-speed[:, 1])/gap, t, initial=0.)
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.2), sharex=True)
    fig.subplots_adjust(left=.15, right=.97, top=.78, bottom=.12, hspace=.24)
    fig.suptitle('One mouth with a free backing', x=.15, y=.975, ha='left', fontsize=15)
    axes[0].plot(t, strain*1e5, color='#2f5c99', linewidth=2.2)
    axes[0].set_ylabel('Mean strain\n(10⁻⁵)')
    for i, (label, color, style) in enumerate((('Near gauge', '#A84135', '-'), ('Deeper gauge', '#24796A', '--'))):
        axes[1].plot(t, tau[:, i]*1e5, color=color, linestyle=style, linewidth=2, label=label)
    axes[1].set_ylabel('Measured traction\n(10⁻⁵ patch units)')
    axes[1].set_xlabel('Time along the converting parcel')
    axes[1].legend(loc='upper right', frameon=False)
    order_axis = axes[0].secondary_xaxis('top')
    order_axis.set_xticks(t[gates], [f'{v:.2f}' for v in order[gates]])
    order_axis.set_xlabel('Local shear order')
    for ax in axes:
        ax.axhline(0, linewidth=.9, color='#4f5866', zorder=0)
        ax.axvspan(0, .45*t[-1], color='#e2eaf1', zorder=-1)
        ax.set_xlim(t[0], t[-1])
        ax.grid(axis='y', color='#d5dce2', lw=.6)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].text(.025, .88, 'Pulse applied', transform=axes[0].transAxes, color='#374955')
    result = Path('software/stage1_solver/reports/figures/neutral_patch_order_pulse.png')
    fig.savefig(result, dpi=160)
    plt.close(fig)
    print(json.dumps(dict(test='direct_exit_readings', positive_drive=True, file=str(source),
                          input_sha256=digest(source), exit_order=order[-1], exit_time=t[-1],
                          near_traction=tau[-1, 0], deep_traction=tau[-1, 1],
                          mean_strain_trapezoidal=strain[-1], figure=str(result),
                          figure_sha256=digest(result)), allow_nan=False))


if __name__ == '__main__':
    main()
