"""Builder controls for the scoped comparator. Execute under s11c_guarded_run.py."""
import unittest

import sympy as s
import S15_reading_comparator as c
from S15_reading_comparator import algebra as m


class ComparatorChecks(unittest.TestCase):
    def test_roster_and_signed_surface_sources(self):
        roster=c.required_names()
        self.assertEqual(len(roster),336)
        sums=[name for name in roster if c.inputs.sum_operands(name)]
        self.assertEqual(len(sums),27)
        self.assertEqual(sum(len(c.inputs.sum_operands(q)) for q in sums),144)
        src=c.inputs.sum_operands('both_surfaces/sum_crossing/component_1')
        self.assertIn('inner/pressure/component_1',src)
        self.assertIn('outer/transport_radius/component_1',src)
        self.assertEqual(c.object_key('outer/carried/component_2/labels_{drain1, drain2}'),
                         'outer/carried/component_2/labels_drain1+drain2')

    def test_live_drain_slot_and_evaluation_point(self):
        py=m.parse('Subs(Derivative(rho_1(_a,t),_a),_a,r_in(t))*u_2(R(t),t)','py')
        wl='Derivative[1,0][rr[1]][rin[t],t]*uu[2][rout[t],t]'
        self.assertEqual(m.exact_residual(py,m.parse(wl,'wl')),0)
        for changed in (wl.replace('[1,0]','[0,1]'),wl.replace('rin[t]','rout[t]'),
                        wl.replace('rr[1]','rr[2]'),'-('+wl+')'):
            other=m.parse(changed,'wl')
            self.assertNotEqual(m.exact_residual(py,other),0)
            ctx=m.Sample(m.FORMS[1],'generic',0,(1,2,3,4),16)
            self.assertNotAlmostEqual(ctx.value(py),ctx.value(other))

    def test_dimension_empty_and_nonzero_sign(self):
        self.assertEqual(c.fraction_rows('[]','py'),c.fraction_rows('{}','wl'))
        self.assertNotEqual(c.fraction_rows('[]','py'),c.fraction_rows('{{0,0,0}}','wl'))
        self.assertNotEqual(c.fraction_rows('[(1,4,-2)]','py'),c.fraction_rows('{{1,4,2}}','wl'))
        inner=m.parse('Integral(-r_in(t)**2*sin(theta),(theta,0,pi))','py')
        outer=m.parse('Integral(R(t)**2*sin(theta),(theta,0,pi))','py')
        self.assertNotEqual(m.exact_residual(inner+outer,outer-inner),0)

    def test_one_sided_branch_change(self):
        py='Piecewise((u_1(r_in(t),t),Ne(a_1,-1)), (u_2(R(t),t),True))'
        wl='Piecewise[{{uu[2][rout[t],t],aa[1]==-1}},uu[1][rin[t],t]]'
        loci=m.branch_loci("['Ne(a_1,-1)','True']",'py')
        self.assertEqual(loci,m.branch_loci('{aa[1]==-1}','wl'))
        self.assertNotEqual(loci,m.branch_loci('{aa[2]==-1}','wl'))
        py=m.parse(py,'py')
        good=m.parse(wl,'wl')
        bad=m.parse(wl.replace('uu[2][rout[t],t]','(-uu[2][rout[t],t])'),'wl')
        for case,exp in m.exponent_samples(0,loci):
            ctx=m.Sample(m.FORMS[0],m.STRATA[0],0,exp,16)
            self.assertEqual(ctx.value(py),ctx.value(good))
            self.assertEqual(ctx.value(py)==ctx.value(bad),case=='generic')

    def test_held_zero_and_live_velocity_change(self):
        self.assertEqual(m.parse('Inactive[Integrate][0,{theta,0,Pi}]','wl'),0)
        ctx=m.Sample(m.FORMS[0],'X1_velocity_zero',0,(1,2,3,4),16)
        self.assertEqual(ctx.value(m.parse('Derivative[1][xx[1,1]][t]','wl')),0)
        self.assertNotEqual(ctx.value(m.parse('Derivative[2][xx[1,1]][t]','wl')),0)
        change=m.parse('centreVelocityChange[1][t]','wl')
        self.assertEqual(change,m.parse('X_1_1_velocity_change(t)','py'))
        self.assertNotEqual(ctx.value(change),0)


if __name__=='__main__':
    unittest.main(verbosity=2)
