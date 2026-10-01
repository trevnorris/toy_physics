"""Pure standard-library tests for continuation decisions/storage, no science."""
import json
from pathlib import Path
import pickle
import tempfile
import unittest
from unittest.mock import patch
import S11c_d_near_unity_uniform_continue as worker


class ArrayStandIn:
    def __init__(self,data,shape=(2,),dtype='f8'):self.data=data;self.shape=shape;self.dtype=dtype
    def tobytes(self):return self.data
    def __eq__(self,other):raise ValueError('Array truth conversion must not be used')


class Tests(unittest.TestCase):
    def test_certificates_fail_closed(self):
        approve=worker.certificate_acceptance
        good=(True,(),True,True,True,0,1)
        self.assertTrue(approve(*good))
        self.assertTrue(approve(True,(),True,True,True,-1,None))
        for at,value in [(0,False),(1,('unbound',)),(2,False),(3,None),(4,False),(6,None),(6,0)]:
            bad=list(good);bad[at]=value;self.assertFalse(approve(*bad))
        self.assertFalse(approve(True,(),True,True,True,None,None))

    def test_exact_array_and_container_identity(self):
        eq=worker.exact_structure
        self.assertTrue(eq({'x':[ArrayStandIn(b'ab')]},{'x':[ArrayStandIn(b'ab')]}))
        self.assertFalse(eq(ArrayStandIn(b'ab'),ArrayStandIn(b'ac')))
        self.assertFalse(eq(ArrayStandIn(b'ab'),ArrayStandIn(b'ab',shape=(1,2))))
        self.assertFalse(eq(ArrayStandIn(b'ab'),ArrayStandIn(b'ab',dtype='i8')))
        self.assertFalse(eq([1],(1,)))
        self.assertFalse(eq({'x':1},{'y':1}))

    def test_all_completed_returns_restored_without_functions(self):
        obj=object.__new__(worker.ContinuationJournal)
        obj.terminals={str(n):dict(status='COMPLETE',input={'member':'i'+str(n)},result={'member':'r'+str(n)}) for n in range(27)}
        obj.terminals['pending']=dict(status='UNRESOLVED')
        obj.prior_inputs={};obj.restored={};obj.restored_count=0;events=[];reads=[]
        def read(ref):reads.append(ref['member']);return ('saved',ref['member'])
        obj.read_prior=read;obj.event=events.append
        obj.restore_completed()
        self.assertEqual(obj.restored_count,27);self.assertEqual(len(reads),54)
        self.assertEqual(len(events),27);self.assertTrue(all(e['functionCalled'] is False for e in events))
        self.assertNotIn('pending',obj.restored)

    def test_argument_join_and_mismatch_are_persisted(self):
        with tempfile.TemporaryDirectory() as directory:
            obj=object.__new__(worker.ContinuationJournal);obj.out=Path(directory)
            obj.prior_store=worker.BlobStore(obj.out/'prior.sqlite',create=True)
            args={'schedule':{'index':2},'sign':1,'source':[1,2]}
            ref=obj.prior_store.put('input',pickle.dumps(args,protocol=5))
            obj.terminals={'point':dict(input=ref)};obj.prior_inputs={'point':args};saved=[]
            obj.emit=lambda name,value:(saved.append((name,value)) or {'name':name})
            self.assertEqual(obj.join('point',args),'EXACT_SERIALIZED_BYTES')
            with self.assertRaises(worker.original.IntegrityError):obj.join('point',args|{'sign':-1})
            evidence=json.loads((obj.out/'argument-mismatch.json').read_text())
            self.assertEqual(evidence['classification'],'INTEGRITY_FAILURE_NOT_PHYSICS')
            self.assertFalse(saved[-1][1]['exactStructuralEqual']);obj.prior_store.close()

    def test_corrupted_prior_bytes_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            store=worker.BlobStore(Path(directory)/'data.sqlite',create=True)
            ref=store.put('x',b'original');bad=ref|{'sha256':'0'*64}
            with self.assertRaises(ValueError):store.get(bad)
            self.assertEqual(store.get(ref),b'original');store.close()

    def test_actual_certificate_policy_with_standin(self):
        class Expr:
            free_symbols=();is_zero=None;is_real=True;is_positive=False;is_negative=False
            def __init__(self,n):self.n=n;self.is_positive=n>0;self.is_negative=n<0;self.is_zero=n==0
            def as_real_imag(self,deep=True):return Expr(0),Expr(self.n)
            def __sub__(self,other):return self
            def __eq__(self,other):return self.n==other
        class OriginalUnknown(Expr):
            def __init__(self,n):super().__init__(n);self.is_zero=None
        class ImaginaryUnit:
            def __mul__(self,other):return other
        class FakeSP:
            I=ImaginaryUnit()
            simplify=staticmethod(lambda x:x)
            cancel=staticmethod(lambda x:0)
        with patch.object(worker,'sp',FakeSP,create=True),patch.object(worker.original,'finite',lambda x:True):
            out=worker.exact_nonzero(OriginalUnknown(5))
            self.assertTrue(out['certified']);self.assertEqual(out['route'],'EXACT_REAL_IMAGINARY_SIGN_WITNESS')
            out=worker.exact_nonzero(OriginalUnknown(0));self.assertFalse(out['certified'])


if __name__=='__main__':unittest.main()
