"""Lossless shared-payload encoding of the S11c-d computed tag stream.

Definitions contain ordinary srepr CAS objects. References point only to an
earlier definition in the same stream; every original tag is still emitted.
Decoding returns the exact original payload text, including float precision.
This is transcript storage, not a replacement for a physical computation.
"""
from collections import defaultdict
import hashlib
import re


DEFINITION = "Tuple(Str('s11cdSharedPayloadDefinition'), Integer("
REFERENCE = "Tuple(Str('s11cdSharedPayloadReference'), Integer("
DEFINITION_PATTERN = re.compile(re.escape(DEFINITION) + r"([0-9]+)\), ")
REFERENCE_PATTERN = re.compile(re.escape(REFERENCE) + r"([0-9]+)\)\)")


class PayloadEncoder:
    def __init__(self, minimum_bytes=192):
        self.minimum_bytes = minimum_bytes
        self.payload_ids = {}

    def encode(self, payload):
        if len(payload) < self.minimum_bytes:
            return payload
        if payload in self.payload_ids:
            return REFERENCE + str(self.payload_ids[payload]) + '))'
        index = len(self.payload_ids)
        self.payload_ids[payload] = index
        return DEFINITION + str(index) + '), ' + payload + ')'


class PayloadDecoder:
    def __init__(self):
        self.payloads = []

    def decode(self, payload):
        match = DEFINITION_PATTERN.match(payload)
        if match:
            index = int(match[1])
            if index != len(self.payloads) or not payload.endswith(')'):
                raise ValueError('invalid shared-payload definition')
            body = payload[match.end():-1]
            self.payloads.append(body)
            return body
        match = REFERENCE_PATTERN.fullmatch(payload)
        if match:
            index = int(match[1])
            if index >= len(self.payloads):
                raise ValueError('undefined shared-payload reference')
            return self.payloads[index]
        if payload.startswith((DEFINITION, REFERENCE)):
            raise ValueError('malformed shared-payload record')
        return payload


def decoded_lines(path):
    """Read both the old expanded stream and the shared-payload stream."""
    decoder = PayloadDecoder()
    with path.open() as stream:
        for line in stream:
            tag, separator, body = line.rstrip('\n').partition(': ')
            if not separator:
                raise ValueError(('invalid transcript line', tag))
            yield tag + separator + decoder.decode(body) + '\n'


def emission_index(lines):
    """Encode tag-order positions as arithmetic progressions per source line."""
    positions = defaultdict(list)
    for index, source_line in enumerate(lines.values()):
        positions[source_line].append(index)
    groups = []
    for source_line, indices in sorted(positions.items()):
        runs = []
        cursor = 0
        while cursor < len(indices):
            start = indices[cursor]
            step = indices[cursor+1]-start if cursor+1 < len(indices) else 1
            count = 1
            while cursor+count < len(indices) and indices[cursor+count] == start+count*step:
                count += 1
            runs.append((start, step, count))
            cursor += count
        groups.append((source_line, tuple(runs)))
    return {'s11cdTagOrderSourceLineProgressions':tuple(groups),
            's11cdIndexedTagCount':len(lines),
            's11cdIndexedTagOrderSha256':hashlib.sha256('\n'.join(lines).encode()).hexdigest()}


def restore_emission_index(record, tags):
    """Reconstruct and validate every tag-to-source-line assignment."""
    count = int(record['s11cdIndexedTagCount'])
    if len(tags) != count:
        raise ValueError('emission index tag-count mismatch')
    digest = hashlib.sha256('\n'.join(tags).encode()).hexdigest()
    if digest != str(record['s11cdIndexedTagOrderSha256']):
        raise ValueError('emission index tag-order mismatch')
    result = [None]*count
    for source_line, runs in record['s11cdTagOrderSourceLineProgressions']:
        for start, step, length in runs:
            start, step, length = map(int, (start, step, length))
            if start < 0 or step < 1 or length < 1:
                raise ValueError('invalid source-line progression')
            for index in range(start, start+step*length, step):
                if index >= count or result[index] is not None:
                    raise ValueError('overlapping or out-of-range source-line progression')
                result[index] = int(source_line)
    if any(line is None for line in result):
        raise ValueError('incomplete source-line index')
    return dict(zip(tags, result))


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('destination',type=Path)
    options = parser.parse_args()
    if options.source.resolve() == options.destination.resolve():
        raise ValueError('expansion requires a separate destination')
    # Exclusive creation also prevents writing through an existing annex link.
    tags = []
    with options.destination.open('x') as output:
        for line in decoded_lines(options.source):
            tag,_,body = line.rstrip('\n').partition(': ')
            if tag == 'PY_S11CD_EMISSION_LINES' and 's11cdTagOrderSourceLineProgressions' in body:
                import sympy as sp
                from sympy.core.symbol import Str
                from ledger_fold import _restore
                record = {str(k):v for k,v in _restore(body)}
                restored = restore_emission_index(record,tags)
                body = sp.srepr(sp.Tuple(*(sp.Tuple(Str(k),sp.Integer(v)) for k,v in restored.items())))
                line = tag+': '+body+'\n'
            output.write(line)
            tags.append(tag)
