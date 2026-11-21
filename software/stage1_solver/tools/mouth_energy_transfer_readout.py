"""Summarize stationary signed power, supply modulation and independent reaction tests."""
import json
from pathlib import Path

from mouth_supplied_contact_work import verified
from mouth_order_wave_pulse import emit

ROOT=Path('_scratch/s11c')


def rows(name,kind):
    return [r for r in verified(ROOT/('contact_supply_'+name+'_01')) if r['test']==kind]


def keyed(records,key):
    return {r[key]:r for r in records}


def main():
    for n in (34,48):
        ref=keyed(rows('steps'+str(n),'supplied_order'),'multiple')
        all_contact=rows('stepcontact'+str(n),'contact')
        for inlet in ('straight','traction_free'):
            contacts=keyed([r for r in all_contact if r['inlet']==inlet],'multiple')
            for delta in (.01,.02):
                lo,hi=round(.5-delta,2),round(.5+delta,2)
                dh=ref[hi]['head_power']-ref[lo]['head_power']
                heat_slope=(ref[hi]['passive_loss']-ref[lo]['passive_loss'])/dh
                work_slope=(ref[hi]['bare_force_work_on_carrier']-ref[lo]['bare_force_work_on_carrier'])/dh
                if abs(heat_slope+work_slope-1)>.001:
                    raise ArithmeticError('Modulation of source power not accounted for')
                low,high=contacts[lo],contacts[hi]
                emit(test='symmetric_supply',mesh=n,inlet=inlet,supply_radius_in_flat_head_units=delta,
                     gain_into_conversion_heat=heat_slope,gain_into_carrier_work=work_slope,
                     supply_power_change=dh,order_work_at_b01=[low['order_work_on_contacts']*1e-4,high['order_work_on_contacts']*1e-4],
                     incremental_gate_power_per_extra_supply_at_b01=1e-4*(high['order_work_on_contacts']-low['order_work_on_contacts'])/dh,
                     segment_force_on_material_at_b01=1e-2*contacts[.5]['segment_force_on_material_x'],
                     segment_force_on_wall_at_b01=1e-2*contacts[.5]['segment_force_on_wall_x'])
            for r in contacts.values():
                if r['order_work_positive']>1e-13 or r['order_work_on_contacts']>=0:
                    raise ArithmeticError('Passive work sign changed')
        fold=rows('branch'+str(n),'branch')
        peak=max(fold,key=lambda r:r['head'])
        after=next(r for r in fold if r['position']==1.26)
        emit(test='neutral_return',mesh=n,top_sampled_head=peak['head'],position=peak['position'],
             supply=peak['head_power'],loss=peak['loss'],unpaid=peak['unpaid_loss'],
             returned_head=after['head'],negative_variation_curvature=after['free_variation_curvature'])
    coupon=rows('coupon','fixed_slip_coupon')
    h=max(r['head'] for r in coupon)
    for beta in (0,1):
        pos=next(r for r in coupon if r['prestrain']==.001 and r['head']==h and r['beta_input']==beta)
        neg=next(r for r in coupon if r['prestrain']==-.001 and r['head']==h and r['beta_input']==beta)
        half=next(r for r in coupon if r['prestrain']==.0005 and r['head']==h and r['beta_input']==beta)
        emit(test='source_attribution',head=h,beta=beta,supply=pos['supply'],conversion_heat=pos['conversion_heat'],
             stress_work_exported=pos['flow_cut_stress_work'],net_delivered=pos['net'],renewal=pos['renewal'],
             spring_energy_flux_gain=pos['increase_in_carried_spring_power'],closure=pos['supply_closure'],
             force_on_material=pos['force_on_material'],force_on_fixture=pos['force_on_fixture'],
             net_at_opposite_prestrain=neg['net'],net_at_half_prestrain=half['net'])
        if not (0<pos['flow_cut_stress_work'] and abs(pos['net']-neg['net'])<1e-13
                and abs(neg['force_on_material']+pos['force_on_material'])<2e-13):
            raise ArithmeticError('Positive head port or strain reversal failed')
    trial=[]
    for name in ('stroke48','strong48','lower48','margin48s3','margin48s4','margin48s5','lowresponse48'):
        trial.extend(rows(name,'stroke'))
    stable=[r for r in trial if r['min_decay_rate']>0]
    if any(r['net_per_renewal']>=1 for r in stable):
        raise ArithmeticError('Stable full-power case: reassess interpretation')
    emit(test='scanned_cases',positive_decay_rows=len(stable),total_rows=len(trial),
         largest_net_renewal_ratio_in_decaying_rows=max(r['net_per_renewal'] for r in stable),
         all_sampled_max_ratio=max(r['net_per_renewal'] for r in trial))
    margin=[r for r in trial if r['multiple']==.05 and r['beta']==5.425]
    emit(test='positive_small_signal',mesh_subdivisions=[r['sub'] for r in margin],
         net_over_renewal=[r['net_per_renewal'] for r in margin],min_decay=[r['min_decay_rate'] for r in margin])
    sweep=keyed([r for r in trial if r['sub']==4 and r['beta']==5.425],'multiple')
    for lo,hi in ((.04,.06),(.045,.055)):
        emit(test='active_supply_control',lo=lo,hi=hi,
             net_per_renewal=[sweep[lo]['net_per_renewal'],sweep[hi]['net_per_renewal']],
             net_power_per_bend_squared=[sweep[lo]['net_conversion_power'],sweep[hi]['net_conversion_power']])


if __name__=='__main__':
    main()
