"""Storage regression and source-only scope checks; no scientific imports."""
import ast
from pathlib import Path
import pickle
import tempfile
import unittest
import S11c_d_numerical_radiating_uniform_refinement as worker
import S11c_d_numerical_radiating_end_maps_v2 as storage
from S11c_d_numerical_radiating_blob_store import BlobStore


class Tests(unittest.TestCase):
    def test_saved_return_is_copied_exactly_and_labelled_without_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=BlobStore(root/'old.sqlite',create=True)
            raw=pickle.dumps({'record':(1,2,3)},protocol=5);ref=source.put('return',raw)
            dest=root/'new';dest.mkdir();journal=storage.Journal(dest)
            result=worker.restored_blob(journal,source,ref,'saved',{'operation':'original-complete'})
            self.assertEqual(result,{'record':(1,2,3)})
            self.assertEqual(journal.store.get(journal.refs['restore/saved']),raw)
            import json
            record=json.loads((dest/'operation-index.jsonl').read_text())
            self.assertFalse(record['functionCalled']);self.assertEqual(record['status'],'RESTORED_PRIOR_COMPLETE_RETURN')
            self.assertEqual(journal.restored,1);journal.store.close();source.close()

    def test_corrupt_return_receipt_refused_before_decode(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=BlobStore(root/'old.sqlite',create=True)
            ref=source.put('return',pickle.dumps(1));ref['sha256']='0'*64
            dest=root/'new';dest.mkdir();journal=storage.Journal(dest)
            with self.assertRaisesRegex(ValueError,'integrity mismatch'):
                worker.restored_blob(journal,source,ref,'saved',{})
            self.assertEqual(journal.count,0);self.assertEqual(journal.store.count(),0)
            journal.store.close();source.close()

    def test_focus_has_no_reference_middle_solve_or_source_producer_calls(self):
        tree=ast.parse(Path(worker.__file__).read_text())
        run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run')
        prohibited={'uniform_gaussian','middle_check','bind_sources','endpoint_orders','source_branch_joins','settings','solve','lstsq','quad_vec','load_packet','Poly'}
        calls={n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id for n in ast.walk(run) if isinstance(n,ast.Call) and isinstance(n.func,(ast.Name,ast.Attribute))}
        self.assertFalse(calls & prohibited)
        self.assertIn('row_integral',calls)
        self.assertNotIn('alarm',calls)


if __name__=='__main__':unittest.main()
