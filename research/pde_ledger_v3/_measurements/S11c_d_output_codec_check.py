#!/usr/bin/env python3
"""Exact payload and source-index round trip of an existing full transcript."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import resource
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from S11c_d_output_codec import PayloadEncoder,PayloadDecoder,emission_index,restore_emission_index
from S11c_d_mixing_scattering_sympy_audit import cas,sp


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('source',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    encoder=PayloadEncoder();decoder=PayloadDecoder();tags=[]
    original=hashlib.sha256();decoded=hashlib.sha256();mismatches=[];bytes_out=0;indices=[]
    source_sha=hashlib.sha256()
    for raw in args.source.open():
        source_sha.update(raw.encode())
        tag,_,body=raw.rstrip('\n').partition(': ')
        if tag=='PY_S11CD_EMISSION_LINES':
            pairs=re.findall(r"Tuple\(Str\('([^']+)'\), Integer\(([0-9]+)\)\)",body)
            lines={name:int(line)for name,line in pairs}
            if list(lines)!=tags:raise ValueError('old source index does not match tag order')
            compact=emission_index(lines)
            restored=restore_emission_index(compact,tags)
            indices.append({'tagCount':len(tags),'differentAssignments':sum(restored[k]!=v for k,v in lines.items()),
                'oldBytes':len(raw),'newBytes':len(sp.srepr(cas(compact)))})
            body=sp.srepr(cas(compact))
        encoded=encoder.encode(body);restored=decoder.decode(encoded)
        if restored!=body:mismatches.append(tag)
        original.update(body.encode());decoded.update(restored.encode())
        bytes_out+=len(tag+': '+encoded+'\n');tags.append(tag)
    malformed=[]
    for payload in ("Tuple(Str('s11cdSharedPayloadReference'), Integer(2))",
                    "Tuple(Str('s11cdSharedPayloadDefinition'), Integer(7), Integer(1))"):
        try:PayloadDecoder().decode(payload)
        except ValueError:malformed.append(True)
        else:malformed.append(False)
    result={'source':str(args.source.resolve()),'sourceSha256':source_sha.hexdigest(),
        'sourceBytes':args.source.stat().st_size,'encodedBytes':bytes_out,'records':len(tags),
        'distinctSharedPayloads':len(encoder.payload_ids),'payloadMismatches':mismatches,
        'unencodedPayloadSha256':original.hexdigest(),'decodedPayloadSha256':decoded.hexdigest(),
        'sourceIndexChecks':indices,'malformedReferenceRejections':malformed,
        'instrumentSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'codecSha256':hashlib.sha256((Path(__file__).resolve().parents[1]/'scripts/S11c_d_output_codec.py').read_bytes()).hexdigest(),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if mismatches or any(r['differentAssignments']for r in indices)or not all(malformed):raise ValueError('codec round trip')


if __name__=='__main__':run()
