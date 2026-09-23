from pathlib import Path
import ast,hashlib,json
root=Path('/var/projects/toy_physics');m=root/'research/pde_ledger_v3/_measurements';index_path=m/'S11c_d_A9_saved_dependency_routes.json';index=json.loads(index_path.read_text())['transcripts'][1]
sha=lambda b:hashlib.sha256(b).hexdigest()
def end_expr(text,start):
 depth=0;quote=None;escape=False
 for i in range(start,len(text)):
  c=text[i]
  if quote:
   if escape:escape=False
   elif c=='\\':escape=True
   elif c==quote:quote=None
  elif c in "'\"":quote=c
  elif c in '([{':depth+=1
  elif c in ')]}':
   if depth==0:return i
   depth-=1
  elif c==',' and depth==0:return i
 raise ValueError('unterminated serialized expression')
def field(text,key):
 marker="Str('"+key+"'), ";start=text.index(marker)+len(marker);stop=end_expr(text,start);return text[start:stop],start,stop
records=[]
with Path(index['path']).open('rb') as stream:
 for record in index['records']:
  stream.seek(record['byteOffset']);raw=stream.read(record['bytesIncludingNewline']);assert sha(raw)==record['sha256'];text=raw.decode();value,start,stop=field(text,'VALUE')
  # This parses literal source syntax only: no eval, constructors or SymPy imports.
  tree=ast.parse(value,mode='eval').body
  for pair in tree.args:
   assert isinstance(pair,ast.Call) and isinstance(pair.args[0],ast.Call) and pair.args[0].func.id=='Integer'
   side=ast.literal_eval(pair.args[0].args[0]);body=ast.get_source_segment(value,pair.args[1]);ident,a,b=field(body,'IDENTIFICATIONS');idt=ast.parse(ident,mode='eval').body
   velocities=[]
   for mapping in idt.args:
    name=ast.literal_eval(mapping.args[0].args[0])
    if '_V_' in name:
     literal=ast.get_source_segment(ident,mapping.args[1]);velocities.append({'target':name,'savedValueSrepr':literal,'literalSha256':sha(literal.encode())})
   records.append({'tag':record['tag'],'face':side,'recordByteOffset':record['byteOffset'],'recordBytes':record['bytesIncludingNewline'],'recordSha256':record['sha256'],'selector':['VALUE',side,'IDENTIFICATIONS'],'velocityIdentifications':velocities,'sourceDerivedQualification':'This is the literal saved velocity identification for this face; no centre elimination or full bulk reconstruction is inferred.'})
result={'status':'LITERAL_SAVED_FACE_DRIVE_NAVIGATION_NOT_SCIENTIFIC_VALIDATION','method':'Read indexed transcript ranges and parse serialization source with Python ast only; no eval, scientific imports or physical computation.','sourceIndex':{'path':str(index_path),'sha256':sha(index_path.read_bytes())},'transcript':{'path':index['path'],'sha256':index['sha256']},'records':records,'scienceCalls':0}
p=m/'S11c_d_A9_saved_face_drives.json'
with p.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'records':len(records),'cases':len({v['tag'] for v in records}),'output':str(p),'bytes':p.stat().st_size}))
