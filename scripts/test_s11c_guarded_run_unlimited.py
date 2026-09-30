"""No science: exercise unlimited duration and retained host-memory protection."""
import copy
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import MagicMock, patch
import s11c_guarded_run_unlimited as guard


class UnlimitedGuard(unittest.TestCase):
    def run_child(self, low_memory=False):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);group=root/'cgroup';group.mkdir()
            for name,value in {'memory.max':str(2*1024**3),'memory.swap.max':'0','pids.max':'32'}.items():
                (group/name).write_text(value)
            spec={**guard.resource_limits(),'seconds':0,'cpu':3,'unit':'synthetic-unit','command':['synthetic-worker'],'cwd':directory}
            manifest=root/'invocation.json';manifest.write_text(json.dumps(spec))
            child=MagicMock();child.pid=123
            child.wait.side_effect=[subprocess.TimeoutExpired('synthetic-worker',2),subprocess.TimeoutExpired('synthetic-worker',2),0]
            available=2*1024**3 if low_memory else 8*1024**3
            with patch.object(guard,'group_path',return_value=group), \
                 patch.object(guard.os,'sched_setaffinity'), \
                 patch.object(guard.os,'sched_getaffinity',return_value={3}), \
                 patch.object(guard.os,'getpriority',return_value=15), \
                 patch.dict(os.environ,{name:'1' for name in guard.THREADS}), \
                 patch.object(guard,'available_memory',return_value=8*1024**3), \
                 patch.object(guard,'sample',return_value={'hostAvailableBytes':available}), \
                 patch.object(guard.subprocess,'check_output',return_value='RuntimeMaxUSec=infinity\nRestart=no\n'), \
                 patch.object(guard.subprocess,'Popen',return_value=child), \
                 patch.object(guard,'stop') as stop, \
                 patch.object(guard.time,'monotonic',side_effect=[0,1000000]):
                code=guard.child_main(manifest)
            return code,json.loads((root/'child-outcome.json').read_text()),child.wait.call_count,stop.call_count

    def test_arbitrarily_large_elapsed_time_does_not_stop_job(self):
        code,outcome,waits,stops=self.run_child()
        self.assertEqual(code,0)
        self.assertIsNone(outcome['guardReason'])
        self.assertEqual(outcome['wallSeconds'],1000000)
        self.assertEqual(waits,3)
        self.assertEqual(stops,1) # final cleanup only

    def test_host_memory_stop_still_operates(self):
        code,outcome,waits,stops=self.run_child(low_memory=True)
        self.assertEqual(code,124)
        self.assertEqual(outcome['guardReason'],'host available memory below 4 GiB for two observations')
        self.assertEqual(waits,2)
        self.assertEqual(stops,2)

    def test_existing_resource_controls_retained(self):
        spec={**guard.resource_limits(),'cpu':3}
        actual={'memory.max':str(2*1024**3),'memory.swap.max':'0','pids.max':'32','nice':15,'affinity':[3],'threads':dict.fromkeys(guard.THREADS,'1')}
        self.assertTrue(guard.limits_match(actual,spec))
        for key,value in [('memory.max','max'),('memory.swap.max','1'),('pids.max','max'),('affinity',[3,4]),('nice',0),('threads',dict.fromkeys(guard.THREADS,'2'))]:
            with self.subTest(key=key):
                changed=copy.deepcopy(actual);changed[key]=value
                self.assertFalse(guard.limits_match(changed,spec))

if __name__=='__main__':unittest.main()
