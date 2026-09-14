#!/usr/bin/env python3
"""Replay unchanged pairing arithmetic with scalar checkpoints and direct errors.

The extra runtime record is diagnostic provenance. It must accompany any
downstream use of the checker packets produced by this wrapper.
"""
import hashlib
import json
from pathlib import Path
import pickle
import resource
import sys
import time
import traceback

import S11c_d_end_pairing_check as checker


def main():
    # Capture the original exception without the system crash-reporting hook.
    sys.excepthook = sys.__excepthook__
    destination = Path(sys.argv[sys.argv.index('--run-directory') + 1])
    state = destination / 'scalar-state'
    state.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    instrument = Path(__file__).resolve()
    signature = {
        'instrument': str(instrument),
        'instrumentSha256': checker.digest(instrument),
        'checkerSignatureSha256': checker.digest(destination / 'signature.json'),
        'python': sys.version,
        'sympy': checker.sp.__version__,
        'stackLimitBytes': list(resource.getrlimit(resource.RLIMIT_STACK)),
        'recursionLimit': sys.getrecursionlimit(),
        'methods': ['rational_coefficient', 'carrier_expansion'],
        'arithmeticSource': 'unchanged native methods; exact-operand cache',
    }
    pin = state / 'signature.json'
    if pin.exists():
        if json.loads(pin.read_text()) != signature:
            raise ValueError('diagnostic resume signature differs')
    else:
        checker.atomic(pin, (json.dumps(signature, indent=2) + '\n').encode())
        checker.atomic(state / instrument.name, instrument.read_bytes())

    def progress(record):
        record = {**record, 'elapsedSeconds': time.monotonic() - started,
                  'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        with (state / 'operations.jsonl').open('a') as stream:
            stream.write(json.dumps(record) + '\n')

    def wrap(name, method):
        def calculate(expression):
            body = pickle.dumps(expression, protocol=5)
            source_sha = hashlib.sha256(body).hexdigest()
            stem = name + '_' + source_sha
            path = state / (stem + '.pickle')
            pin = state / (stem + '.sha256')
            caller = sys._getframe(1)
            context = {'method': name, 'sourceSha256': source_sha,
                       'caller': [caller.f_code.co_filename, caller.f_code.co_name, caller.f_lineno]}
            if path.exists():
                if checker.digest(path) != pin.read_text().strip():
                    raise ValueError(('diagnostic scalar digest', stem))
                saved, value, known = pickle.loads(path.read_bytes())
                if saved != expression:
                    raise ValueError(('diagnostic scalar operand', stem))
                checker.engine.PHYSICAL_METADATA.dimensions.known.update(known)
                progress({**context, 'stage': 'resumed'})
                return value
            checker.atomic(state / (stem + '.input.pickle'), body)
            progress({**context, 'stage': 'started', 'inputBytes': len(body)})
            try:
                value = method(expression)
            except BaseException as error:
                detail = {**context, 'stage': 'exception', 'exceptionType': type(error).__name__,
                          'message': str(error), 'traceback': traceback.format_exc()}
                checker.atomic(state / 'exception.json', (json.dumps(detail, indent=2) + '\n').encode())
                progress({**context, 'stage': 'exception', 'exceptionType': type(error).__name__})
                raise
            checker.atomic(path, pickle.dumps((expression, value,
                checker.engine.PHYSICAL_METADATA.dimensions.known), protocol=5))
            checker.atomic(pin, (checker.digest(path) + '\n').encode())
            progress({**context, 'stage': 'saved', 'sha256': checker.digest(path)})
            return value
        return calculate

    for name in signature['methods']:
        original = getattr(checker.engine.ClosedCurrentPairing, name)
        setattr(checker.engine.ClosedCurrentPairing, name, staticmethod(wrap(name, original)))
    try:
        checker.run()
    except BaseException as error:
        detail = {'exceptionType': type(error).__name__, 'message': str(error),
                  'traceback': traceback.format_exc()}
        checker.atomic(state / 'run-exception.json', (json.dumps(detail, indent=2) + '\n').encode())
        progress({'stage': 'run_exception', 'exceptionType': type(error).__name__})
        raise
    progress({'stage': 'complete'})


if __name__ == '__main__':
    main()
