"""Roundoff-panel and immutable continuation tests; no scientific objects."""
import json
import math
from pathlib import Path
import tempfile
import unittest
import S11c_d_numerical_radiating_integration_v3 as integration
import S11c_d_numerical_radiating_end_maps_v2 as storage
from S11c_d_numerical_radiating_saved_operations import SavedOperationsJournal

class Tests(unittest.TestCase):
    def test_actual_synthetic_panel_collapse_and_correction(self):
        points={-4.,-.2,.2,4.}
        for c in (-.31,.11,.41):
            points.add(c);d=.01
            while d<2*(8+abs(c)):
                points.update(v for v in (c-d,c+d) if -4<v<4);d*=2
        def mapped_gap(lo,hi):
            if lo>=-.2 and hi<=.2:return math.asin(hi/.2)-math.asin(lo/.2)
            a,b=sorted(math.acosh(max(1,abs(v)/.2)) for v in (lo,hi));return b-a
        old=sorted(points);self.assertTrue(any(mapped_gap(a,b)==0 for a,b in zip(old[:-1],old[1:])))
        new=integration.stable_panels(points,(-4.,-.2,.2,4.));self.assertTrue(all(mapped_gap(a,b)>0 for a,b in zip(new[:-1],new[1:])))
        self.assertEqual(new[0],-4.);self.assertEqual(new[-1],4.);self.assertTrue({-.2,.2}<=set(new));self.assertAlmostEqual(sum(b-a for a,b in zip(new[:-1],new[1:])),8.,places=14)
    def test_distinct_cuts_preserved_and_critical_endpoints_win(self):
        points=[-4.,-.2,0.,math.nextafter(.2,0),.2,math.nextafter(.2,1),.20001,4.]
        new=integration.stable_panels(points,(-4.,-.2,.2,4.));self.assertEqual(new,[-4.,-.2,0.,.2,.20001,4.])
    def test_completed_operations_not_replayed_pending_input_joined(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);old=root/'old';old.mkdir();j=storage.Journal(old)
            j.call('complete-a',{'x':1},lambda:{'answer':2});j.call('complete-b',{'x':2},lambda:3)
            try:j.call('unfinished',{'x':3},lambda:storage.require(False,'fixture'))
            except ValueError:pass
            j.blob('non-journal-evidence',{'saved':4});j.store.close()
            files={str(p.relative_to(old)):{'sha256':storage.digest(p),'bytes':p.stat().st_size} for p in old.rglob('*') if p.is_file()}
            new=root/'new';new.mkdir();k=SavedOperationsJournal(new,{'path':str(old),'files':files,'completeOperations':2})
            def fail():raise AssertionError('replayed')
            self.assertEqual(k.call('complete-a',{'x':1},fail),{'answer':2});self.assertEqual(k.call('complete-b',{'x':2},fail),3)
            self.assertEqual(k.call('unfinished',{'x':3},lambda:5),5);self.assertEqual(k.restored,2);self.assertIsNone(k.pending)
            k.store.close();k.saved_store.close()

if __name__=='__main__':unittest.main()
