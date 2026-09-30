"""Synthetic storage tests; no scientific imports or saved payloads."""
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile

from S11c_d_numerical_radiating_blob_store import BlobStore


def main():
    checks = {}
    with tempfile.TemporaryDirectory(prefix='s11c-byte-store-') as directory:
        path = Path(directory) / 'operations.sqlite'
        writer = BlobStore(path, create=True)
        receipts = []
        for index in range(512):
            payload = b'\x00\xff' + str(index).encode() + bytes(range(256))
            receipt = writer.put('synthetic/' + str(index), payload)
            assert writer.get(receipt) == payload
            receipts.append(receipt)
        checks['all512OpaqueByteRoundtrips'] = True
        try:
            writer.put(receipts[0]['member'], b'replacement')
        except sqlite3.IntegrityError:
            pass
        else:
            raise AssertionError('Duplicate overwrite was accepted')
        assert writer.count() == 512
        assert writer.get(receipts[0]).startswith(b'\x00\xff0')
        checks['duplicateRefusedWithoutChangingPriorBlob'] = True
        writer.integrity_check()
        writer.close()

        # Abruptly end a separate process during an uncommitted write.
        # Recovery must preserve committed bytes and omit the incomplete blob.
        program = '''import sqlite3,sys,os
c=sqlite3.connect(sys.argv[1]);c.execute('BEGIN IMMEDIATE')
c.execute('INSERT INTO blobs VALUES (?,?,?,?)',('uncommitted',b'x','bad',1))
os._exit(23)
'''
        child = subprocess.run([sys.executable, '-c', program, str(path)],
                               capture_output=True, check=False)
        assert child.returncode == 23 and not child.stderr
        # Read-write open solely to permit SQLite's crash-journal recovery.
        with sqlite3.connect(path) as recovery:
            assert recovery.execute('SELECT count(*) FROM blobs').fetchone()[0] == 512
        reader = BlobStore(path)
        reader.integrity_check()
        for receipt in receipts:
            payload = reader.get(receipt)
            assert hashlib.sha256(payload).hexdigest() == receipt['sha256']
        checks['crashRollbackAndAllCommittedBlobsPreserved'] = True
        try:
            reader.put('read-only-write', b'x')
        except sqlite3.OperationalError:
            pass
        else:
            raise AssertionError('Read-only connection accepted a write')
        reader.close()
        checks['readOnlyReopenEnforced'] = True
        try:
            BlobStore(path, create=True)
        except FileExistsError:
            pass
        else:
            raise AssertionError('Existing database was overwritten')
        checks['existingFileCreateRefused'] = True

        with sqlite3.connect(path) as tamper:
            tamper.execute('UPDATE blobs SET payload=? WHERE name=?',
                           (b'changed', receipts[0]['member']))
        reader = BlobStore(path)
        try:
            reader.get(receipts[0])
        except ValueError:
            pass
        else:
            raise AssertionError('Corrupted payload was accepted')
        reader.close()
        checks['contentCorruptionDetected'] = True
    print(json.dumps({'status': 'PASS', 'checks': checks,
                      'scientificInputsOrImports': False}, indent=2))


if __name__ == '__main__':
    main()
