"""Configuration isolation and scientific-body preservation checks."""
import ast
from pathlib import Path
import unittest
from S11c_d_numerical_radiating_balance import matrix_group_setting

class Tests(unittest.TestCase):
    def test_only_single_leg_order_changes_without_mutating_source(self):
        setting={'outerOrder':24,'singleLegOuterOrder':48,'innerOrders':(8,8),'momentumBound':4.}
        single=matrix_group_setting(setting,('k',));self.assertEqual(single['outerOrder'],48)
        self.assertEqual(single['innerOrders'],(8,8));self.assertEqual(single['momentumBound'],4.)
        self.assertEqual(setting['outerOrder'],24)
        for variables in ((),('k','q'),('k','q','p')):
            self.assertEqual(matrix_group_setting(setting,variables),setting)
    def test_finest_order_preserves_nested_refinement(self):
        setting={'outerOrder':32,'singleLegOuterOrder':64,'innerOrders':(12,12),'regulator':.05}
        self.assertEqual(matrix_group_setting(setting,('k',))['outerOrder'],64)
        self.assertEqual(matrix_group_setting(setting,('k','p')),setting)
    def test_no_completed_integration_or_producer_calls(self):
        source=Path(__file__).with_name('S11c_d_numerical_radiating_balance.py').read_text();tree=ast.parse(source)
        funcs={n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,(ast.Name,ast.Attribute))}
        self.assertFalse(funcs & {'uniform_gaussian','middle_check','bind_sources','endpoint_orders','radical_inventory','source_branch_joins','quad_vec','alarm'})
        self.assertTrue({'matrix_group','matrices','solve','control_report'} <= funcs)

if __name__=='__main__':unittest.main()
