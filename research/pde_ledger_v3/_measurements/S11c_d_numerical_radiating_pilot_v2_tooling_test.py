"""Local startup-order/reuse regression tests, synthetic stdlib objects only."""
import ast
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import S11c_d_numerical_radiating_pilot_v2 as new
import S11c_d_numerical_radiating_end_maps_v2 as storage

M=Path(__file__).resolve().parent

class Tests(unittest.TestCase):
    def execute_prefix(self,path):
        tree=ast.parse(path.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run');prefix=[]
        for statement in run.body:
            if isinstance(statement,ast.Expr) and isinstance(statement.value,ast.Call) and ast.unparse(statement.value.func)=='J.report':break
            prefix.append(statement)
        class Helpers:
            ready=False
            def native_helpers(self,*args):self.ready=True;return SimpleNamespace()
            def context(self,*args):
                if not self.ready:raise RuntimeError('numerical dependencies not initialized')
                return 'context'
        class Journal:
            def call(self,name,args,fn):return {'bound':{}} if name=='restore/native' else {}
        env={'I':Helpers(),'F':SimpleNamespace(initialize=lambda *x:None),'J':Journal(),'np':object(),'sp':object(),'manifest':{'packets':{k:{} for k in ('native','frequency','grades')},'endMaps':{},'nativeEngine':'engine','nativeFinite':'finite'}}
        exec(compile(ast.Module(body=prefix,type_ignores=[]),str(path),'exec'),env)
        self.assertEqual(env['r'],'context')
    def test_regression_original_fails_corrected_initializes_first(self):
        with self.assertRaisesRegex(RuntimeError,'not initialized'):self.execute_prefix(M/'S11c_d_numerical_radiating_pilot.py')
        self.execute_prefix(Path(new.__file__))
    def test_four_saved_reads_restored_without_functions(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);old=root/'old';old.mkdir();j=storage.Journal(old)
            for i,n in enumerate(('frequency','grades','native','completed-end-maps')):j.call('restore/'+n,{'index':i},lambda i=i:{'saved':i})
            ref=j.blob('inputs/completed-end-maps.pickle',{'end':1});storage.save(old/'consumed-end-map-evidence.json',{'receipt':ref});j.store.close()
            files={str(p.relative_to(old)):{'sha256':storage.digest(p),'bytes':p.stat().st_size} for p in old.rglob('*') if p.is_file()};dest=root/'new';dest.mkdir();k=new.SavedStartupJournal(dest,{'path':str(old),'files':files})
            def forbidden():raise AssertionError('completed function replayed')
            for i,n in enumerate(('frequency','grades','native','completed-end-maps')):self.assertEqual(k.call('restore/'+n,{'index':i},forbidden),{'saved':i})
            self.assertEqual(k.restored,4);self.assertEqual(k.call('new-work',{},lambda:5),5)
            for rel,item in files.items():self.assertEqual(storage.digest(dest/'prior-startup'/rel),item['sha256'])
            k.store.close()
    def test_actual_read_input_mismatch_is_fatal(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);old=root/'old';old.mkdir();j=storage.Journal(old)
            for i,n in enumerate(('frequency','grades','native','completed-end-maps')):j.call('restore/'+n,{'i':i},lambda:1)
            j.blob('inputs/completed-end-maps.pickle',1);storage.save(old/'consumed-end-map-evidence.json',{});j.store.close()
            files={str(p.relative_to(old)):{'sha256':storage.digest(p),'bytes':p.stat().st_size} for p in old.rglob('*') if p.is_file()};dest=root/'new';dest.mkdir();k=new.SavedStartupJournal(dest,{'path':str(old),'files':files})
            with self.assertRaisesRegex(ValueError,'input byte identity'):k.call('restore/frequency',{'i':99},lambda:1)
            self.assertEqual(k.restored,0);self.assertTrue((dest/'active-operation.json').exists());k.store.close();k.saved_store.close()

if __name__=='__main__':unittest.main()
