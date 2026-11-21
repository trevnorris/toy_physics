"""Predict a short strain pulse from fixed-order hold forces, then score it.

The prediction command reads only the held steps. Frequencies, impulse weights,
and force scales are never fitted to the validation drive. Strain is the change
in relative mean displacement in gap units; both planes move symmetrically.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve


LAWS = ('central', 'layered')
AMPLITUDES = [0., .001, -.001, .0005, -.0005]
PULSE_FRACTION = .25


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def rms(a):
    return np.sqrt(np.mean(a*a, axis=-1))


def apply_hold(kernel, displacement):
    """Integrate the measured hold kernel against changes of imposed strain."""
    change = np.diff(displacement)
    middle = (kernel[..., 1:] + kernel[..., :-1])/2
    shape = (1,)*(middle.ndim-1) + change.shape
    convolution = fftconvolve(middle, change.reshape(shape), axes=-1)
    return np.concatenate([np.zeros_like(kernel[..., :1]),
                           convolution[..., :len(change)]], axis=-1)


def kinematics(time):
    length = time[-1]*PULSE_FRACTION
    strain = np.where(time <= length, np.sin(np.pi*time/length)**2, 0.)
    rate = np.where(time <= length, (np.pi/length)*np.sin(2*np.pi*time/length), 0.)
    return strain, rate


def read_material(path, phase):
    with np.load(path, allow_pickle=False) as data:
        force = data['force']
        time = data['time']
        orders, parents, seeds = (data[k] for k in ('orders','parent_index','seeds'))
        amplitudes = data['amplitudes']
        np.testing.assert_allclose(amplitudes, AMPLITUDES, rtol=0, atol=0)
        np.testing.assert_array_equal(parents, np.repeat(np.arange(4), 2))
        if len(set(seeds)) != len(parents) or not np.all(np.isfinite(force)):
            raise ValueError(f'Invalid thermal histories in {path}')
        if force.shape != (2, len(orders), len(parents), len(amplitudes), len(time)):
            raise ValueError(f'Unexpected force dimensions in {path}')
        np.testing.assert_allclose(np.diff(time), time[-1]/(len(time)-1), rtol=1e-10)
        np.testing.assert_allclose(time[0], 0., atol=1e-14)
        shape = np.ones_like(time) if phase == 'hold' else kinematics(time)[0]
        np.testing.assert_allclose(data['strain'], amplitudes[:, None]*shape, atol=1e-13)
        if phase == 'pulse':
            np.testing.assert_allclose(data['bath'], 6*amplitudes[:,None]*kinematics(time)[1], atol=1e-12)
        pair = (force[:,:,:,1] - force[:,:,:,2])/(2*amplitudes[1])
        half = (force[:,:,:,3] - force[:,:,:,4])/(2*amplitudes[3])
        even = (force[:,:,:,1]+force[:,:,:,2]-2*force[:,:,:,0])/(2*amplitudes[1])
    return time, orders, parents, seeds, pair, half, even


def predict(folder, frozen):
    saved = {}
    pins = {}
    for law in LAWS:
        path = folder/f'hold_{law}.npz'
        time, orders, parents, seeds, pair, half, even = read_material(path, 'hold')
        if saved:
            for key, value in [('time',time), ('orders',orders), ('parent_index',parents), ('seeds',seeds)]:
                np.testing.assert_array_equal(saved[key], value)
        shape, _ = kinematics(time)
        # Constant and exponentially relaxing material give analytic solutions
        # for this exact waveform. The first also exercises the recovery tail.
        np.testing.assert_allclose(apply_hold(np.ones_like(time), shape), shape, atol=2e-13)
        rate = 2/time[-1]
        omega = 2*np.pi/(PULSE_FRACTION*time[-1])
        end = np.minimum(time, PULSE_FRACTION*time[-1])
        exact = omega/2 * ((rate*np.sin(omega*end)-omega*np.cos(omega*end) +
                            omega*np.exp(-rate*end))/(rate**2+omega**2)) * np.exp(-rate*(time-end))
        control_error = np.max(abs(apply_hold(np.exp(-rate*time), shape)-exact))
        if control_error > 3e-6:
            raise ArithmeticError('Analytic relaxation check failed')
        pins[law] = sha(path)
        saved.update(time=time, orders=orders, parent_index=parents, seeds=seeds, strain=shape)
        saved[law+'_held'] = pair
        saved[law+'_predicted'] = apply_hold(pair, shape)
        end_window = time >= .9*time[-1]
        for i, order in enumerate(orders):
            k = pair[-1, i].mean(axis=0)
            retained = k/k[0]
            # Report a first crossing only if it persists for 5% of the transit.
            width = max(1, int(np.ceil(.05*(len(time)-1))))
            crossings = np.flatnonzero(np.convolve((retained<=.5).astype(int),
                                        np.ones(width,dtype=int), mode='valid') == width)
            last_ratio = pair[-1,i][:,end_window].mean(axis=1)/pair[-1,i,:,0]
            emit(test='strain_hold', law=law, order=order, initial_force_per_unit_strain=k[0],
                 retained_at_tenth= float(np.interp(.1*time[-1],time,retained)),
                 retained_at_quarter=float(np.interp(.25*time[-1],time,retained)),
                 retained_last_tenth=float(np.mean(retained[end_window])),
                 individual_tail_ratio_range=[float(last_ratio.min()),float(last_ratio.max())],
                 first_sustained_half_time=float(time[crossings[0]]) if crossings.size else None,
                 half_drive_relative_rms=rms((pair-half)[-1,i].mean(axis=0))/rms(k),
                 even_relative_rms=rms(even[-1,i].mean(axis=0))/rms(k),
                 refinement_relative_rms=rms(k-pair[0,i].mean(axis=0))/rms(k))
        emit(test='pulse_forecast_control',law=law, exponential_control_max_abs_error=control_error)
    saved['hold_pins'] = json.dumps(pins)
    with frozen.open('xb') as f:
        np.savez_compressed(f, **saved)
    emit(test='forecast_frozen',source_sha256=sha(Path(__file__)),sha256=sha(frozen),hold_inputs=pins,
         pulse_duration=PULSE_FRACTION*time[-1],hold_duration=time[-1],
         diagnostic='10% relative RMS of signed internal force over full record; peak and tail reported separately')


def score(folder, frozen, expected, figure):
    if sha(frozen) != expected:
        raise ValueError('Frozen forecast digest changed')
    with np.load(frozen, allow_pickle=False) as forecasts:
        pins = json.loads(str(forecasts['hold_pins']))
        plot_rows = []
        for law in LAWS:
            if sha(folder/f'hold_{law}.npz') != pins[law]:
                raise ValueError('Held-strain reference changed')
            path = folder/f'pulse_{law}.npz'
            t, orders, parents, seeds, obs, half, even = read_material(path, 'pulse')
            for key, value in [('time',t), ('orders',orders), ('parent_index',parents)]:
                np.testing.assert_array_equal(forecasts[key], value)
            if set(seeds) & set(forecasts['seeds']):
                raise ValueError('Validation reuses a future bath from calibration')
            emit(test='pulse_source',law=law,sha256=sha(path),forecast_sha256=expected)
            shape, rate = kinematics(t)
            drive_on = t <= PULSE_FRACTION*t[-1]
            for i, order in enumerate(orders):
                per_history = forecasts[law+'_predicted'][-1,i]
                held = forecasts[law+'_held'][-1,i].mean(axis=0)
                prediction = per_history.mean(axis=0)
                measured = obs[-1,i].mean(axis=0)
                residual = prediction-measured
                scale = rms(measured)
                loo = []
                for j in np.unique(parents):
                    pred = per_history[parents!=j].mean(axis=0)
                    direct = obs[-1,i,parents!=j].mean(axis=0)
                    loo.append(float(rms(pred-direct)/rms(direct)))
                emit(test='pulse_score',law=law,order=order,relative_rms=float(rms(residual)/scale),
                     direct_peak=float(np.max(measured)),predicted_peak=float(np.max(prediction)),
                     late_rms_rel_to_whole_force=float(rms(residual[~drive_on])/scale),
                     late_observed_rms_rel_to_whole_force=float(rms(measured[~drive_on])/scale),
                     no_relaxation_relative_rms=float(rms(held[0]*shape-measured)/scale),
                     half_drive_relative_rms=float(rms(half[-1,i].mean(axis=0)-measured)/scale),
                     even_relative_rms=float(rms(even[-1,i].mean(axis=0))/scale),
                     direct_step_refinement=float(rms(measured-obs[0,i].mean(axis=0))/scale),
                     forecast_step_refinement=float(rms(prediction-forecasts[law+'_predicted'][0,i].mean(axis=0))/scale),
                     leave_one_parent_out_range=[min(loo),max(loo)],
                     independent_bath_split_relative_rms=float(rms(obs[-1,i,::2].mean(axis=0)-obs[-1,i,1::2].mean(axis=0))/scale),
                     total_with_direct_drag_relative_rms=float(rms(residual)/rms(measured+6*rate)))
                plot_rows.append((law, order, t, held, prediction, measured))
    draw(figure, plot_rows)
    emit(test='figure',path=str(figure),sha256=sha(figure))


def draw(path, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {0.95:'#244d86', .5:'#b64d18', .1:'#187d63', .01:'#844887'}
    fig, axs = plt.subplots(2, 2, figsize=(9.2, 5.7), sharex=True)
    for law, order, t, held, pred, measured in rows:
        j = LAWS.index(law)
        color = colors[round(float(order),2)]
        axs[0,j].plot(t,held/held[0],lw=1.4,label=f'{order:.2f}',color=color)
        # Show the half-order force on its own signed scale.
        if np.isclose(order,.5):
            axs[1,j].plot(t, measured, color='#9c462a', lw=1.7,label='Independent pulse')
            axs[1,j].plot(t, pred, color='#164e72', lw=1.7,ls='--',label='From held strain')
            axs[1,j].axvspan(0,.25*t[-1],color='#e6ecf0',zorder=-2)
        axs[0,j].set_title(law.capitalize()+' bonds')
    for ax in axs.flat:
        ax.axhline(0,color='#717780',lw=.65,zorder=-1)
        ax.grid(axis='y',alpha=.25)
        ax.spines[['top','right']].set_visible(False)
        ax.set_xlim(0,rows[0][2][-1])
    axs[0,0].legend(title='Shear order',ncol=4,frameon=False,fontsize=8)
    axs[1,1].legend(frameon=False,fontsize=9)
    axs[0,0].set_ylabel('Held force / initial force')
    axs[1,0].set_ylabel('Force / peak strain\nat half order')
    for ax in axs[1]:
        ax.set_xlabel('Time (carrier units)')
    fig.tight_layout()
    fig.savefig(path,dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['predict','score'])
    parser.add_argument('folder', type=Path)
    parser.add_argument('frozen', type=Path)
    parser.add_argument('--expected')
    parser.add_argument('--figure',type=Path)
    args = parser.parse_args()
    if args.phase=='predict':
        predict(args.folder,args.frozen)
    else:
        score(args.folder,args.frozen,args.expected,args.figure)
