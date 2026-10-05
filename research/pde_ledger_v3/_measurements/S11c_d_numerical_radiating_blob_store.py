"""Immutable byte checkpoints without rescanning an archive per append.

Pure storage: no scientific imports, decoding, transformations or equations.
One SQLite connection, durable transaction per blob, unique names and SHA256.
This helper is prepared for future use; it does not resume the stopped job.
"""
import hashlib
from pathlib import Path
import sqlite3


class BlobStore:
    def __init__(self, path, *, create=False):
        self.path = Path(path).resolve()
        if create:
            # Refuse overwrite before SQLite gets a chance to open the file.
            with self.path.open('xb'):
                pass
            self.connection = sqlite3.connect(self.path)
            self.connection.execute('PRAGMA journal_mode=DELETE')
            self.connection.execute('PRAGMA synchronous=FULL')
            with self.connection:
                self.connection.execute('''CREATE TABLE blobs (
                    name TEXT PRIMARY KEY,
                    payload BLOB NOT NULL,
                    sha256 TEXT NOT NULL,
                    bytes INTEGER NOT NULL CHECK(bytes >= 0)
                )''')
        else:
            self.connection = sqlite3.connect(
                self.path.as_uri() + '?mode=ro', uri=True)

    def put(self, name, payload):
        if not isinstance(name, str) or not name:
            raise ValueError('Nonempty immutable blob name required')
        if not isinstance(payload, bytes):
            raise TypeError('Only already-encoded bytes are accepted')
        digest = hashlib.sha256(payload).hexdigest()
        # A failed/duplicate write rolls back without changing a prior blob.
        with self.connection:
            self.connection.execute(
                'INSERT INTO blobs VALUES (?, ?, ?, ?)',
                (name, payload, digest, len(payload)))
        return {'storage': 'sqlite', 'member': name,
                'sha256': digest, 'bytes': len(payload)}

    def get(self, receipt):
        row = self.connection.execute(
            'SELECT payload, sha256, bytes FROM blobs WHERE name=?',
            (receipt['member'],)).fetchone()
        if row is None:
            raise KeyError(receipt['member'])
        payload, digest, size = row
        if (digest != receipt['sha256'] or size != receipt['bytes']
                or size != len(payload)
                or hashlib.sha256(payload).hexdigest() != digest):
            raise ValueError('Blob/receipt integrity mismatch')
        return payload

    def count(self):
        return self.connection.execute('SELECT count(*) FROM blobs').fetchone()[0]

    def integrity_check(self):
        rows = self.connection.execute('PRAGMA integrity_check').fetchall()
        if rows != [('ok',)]:
            raise ValueError('SQLite structural integrity check failed: ' + repr(rows))

    def close(self):
        self.connection.close()
