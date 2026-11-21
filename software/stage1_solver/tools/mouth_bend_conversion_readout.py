"""Compare matched-intake conversion savings, leading order response and contact costs."""
import json
from pathlib import Path

from mouth_flared_receiver import emit

ROOT=Path('_scratch/s11c')


def rows(name,test):
    folder=ROOT/name
    guard=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if (guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Incomplete or unverified input {folder}')
    parsed=[json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    return [r for r in parsed if r.get('test')==test]


def main():
    by_run={}
    for label in ('n34_02','n48plus_02','n48minus_02','modes4_accepted01','weak48_01'):
        data=rows('bend_conversion_'+label,'matched_intake')
        by_run[label]={r['bend']:r for r in data}
        for r in data:
            if (abs(r['feed']-.25)>2e-9 or abs(r['bulk']-.25)>2e-9 or
                    r['order_correction']>3e-12 or r['field_residual']>3e-11):
                raise ArithmeticError(f'Carrier or order unresolved in {label}')
    positive=by_run['n48plus_02']
    for b in (.005,.01):
        a=positive[b]
        m=by_run['n48minus_02'][-b]
        coarse=by_run['n34_02'][b]
        extra=by_run['modes4_accepted01'][b]
        emit(test='finite_comparison',bend=b,baseline=positive[0.]['passive_loss'],
             remaining=a['passive_loss'],saved=a['saved_power'],mirror_difference=a['saved_power']-m['saved_power'],
             change_on_mesh_refinement=a['saved_power']-coarse['saved_power'],
             coarse_change_on_angular_refinement=extra['saved_power']-coarse['saved_power'],
             conditional_force_on_flow_x=a['force_on_carrier_xw'][0],
             conditional_force_on_flow_x_mirror=m['force_on_carrier_xw'][0])
    weak=by_run['weak48_01']
    sample={b:r['saved_power']/b**2 for b,r in (positive|weak).items() if b>0}
    # Removing the first even finite-amplitude error yields independent
    # quadratic-coefficient estimates from successive pairs of amplitudes.
    estimates=[(4*sample[b/2]-sample[b])/3 for b in (.01,.005,.0025)]
    emit(test='weak_limit',savings_per_bend_squared=sample,extrapolations=estimates)
    if max(abs(x/estimates[-1]-1) for x in estimates)>.01:
        raise ArithmeticError('Inconsistent weak limit')
    contact={r['inlet']:r for r in rows('contact_order_softening_m48sub4_01',None)
             if r.get('bend')==.01}
    for row in rows('bend_conversion_gate48fine_01','order_feedback'):
        name=row['inlet']
        other=contact['relaxed' if name=='traction_free' else name]
        gate=sum(row['derivative_loss_per_bend_squared'])
        saving=estimates[-1]-gate
        renewal=other['renewal_loss']/.01**2
        spring=other['spring_energy_in']/.01**2
        emit(test='conditional_budget',inlet=name,bare_saving_per_bend_squared=estimates[-1],
             extra_loss_from_strain_per_bend_squared=gate,net_saving_per_bend_squared=saving,
             renewal_per_bend_squared=renewal,saving_over_renewal=saving/renewal,
             fraction_if_all_renewal=renewal/saving,fraction_after_incoming_spring_credit=(renewal-spring)/saving,
             power_leaving_softened_springs_at_b01=other['softening_release'],
             material_slip_work_at_b01=other['material_slip_work'])
    fine=[r for r in rows('bend_conversion_cuts_01','order_energy_cuts') if r['mesh']==48]
    for r in fine:
        if r['bend'] in (0.,.01):
            emit(test='material_cut_summary',bend=r['bend'],
                 distant_order_power_out=sum(v['sum'] for v in r['baths'].values()),
                 middle_enter_leave=[r['segments']['(1.36, 1.52)'][side]['sum'] for side in ('enter','leave')],
                 raw_sidewall_power_error=r['segments']['(1.173, 2.067)']['wall']['sum'])


if __name__=='__main__':
    main()
