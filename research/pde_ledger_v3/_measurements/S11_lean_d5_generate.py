#!/usr/bin/env python3
"""D5 finite quadratic completeness certificates, independent of native Q9.

SymPy selects rational necessary constraints. Lean checks the exhaustive
representation, all retained equations, reconstruction and full-group
sufficiency. This generator is not a trusted mathematical oracle.
"""
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
import sys
import re

OUT=Path(__file__).resolve().parents[1]/'lean/s11/S11D5Invariants'
PAIRS=list(combinations_with_replacement(range(25),2))
def vec(xs):return '!['+','.join(map(str,xs))+']'
def mat(xs):return '!['+','.join(vec(xs[5*i:5*i+5]) for i in range(5))+']'
def rat(x):
 x=s.Rational(x)
 return str(x.p) if x.q==1 else f'({x.p}/{x.q})'
def write(name,text):
 path=OUT/(name+'.lean')
 if '--check' in sys.argv:assert path.read_text()==text,f'Generated source drift: {path}'
 else:path.write_text(text)

def header(imports,doc,heartbeats=1600000):return imports+'\n\n/-! '+doc+' -/\nnamespace S11D5Invariants\nnoncomputable section\nset_option maxRecDepth 4096\nset_option maxHeartbeats '+str(heartbeats)+'\n\n'
end='\nend\nend S11D5Invariants\n'

def split_reconstruction(monolithic,necessary_blocks,free):
 """Keep every target and exact linear combination, splitting only proof units."""
 steps=[(int(k),target,combo) for k,target,combo in re.findall(
  r'  have h(\d+) : (.*?) := by\n    linear_combination (.*)\n',monolithic)]
 assert [k for k,_,_ in steps]==[k for k in range(325) if k not in free]
 groups=[steps[i:i+20] for i in range(0,len(steps),20)]
 equation_block={e:bi for bi,block in enumerate(necessary_blocks) for e in block}
 all_used=set()
 for bi,group in enumerate(groups):
  used={int(e) for _,_,combo in group for e in re.findall(r'\be(\d+)\b',combo)}
  all_used|=used
  required_blocks=sorted({equation_block[e] for e in used})
  part=header('\n'.join(f'import S11D5Invariants.ConstraintBlock{j}' for j in required_blocks),
   'Bounded coefficient reconstruction using unchanged exact linear combinations.',3200000)
  part+=f'theorem reconstructedBlock{bi} {{Q : Quad}} (hQ : OInvariant Q) (c : Coefficients)\n'
  part+='    (hc : ∀ G, Q G = polynomial c (coordinates G)) :\n'
  part+='    '+' ∧\n    '.join('('+target+')' for _,target,_ in group)+' := by\n'
  for j in required_blocks:
   names=','.join(f'e{e}' if e in used else '_' for e in necessary_blocks[j])
   part+=f'  obtain ⟨{names}⟩ := necessaryBlock{j} hQ c hc\n'
  part+='  refine ⟨'+','.join('?_ ' for _ in group).replace(' ','')+'⟩\n'
  for _,_,combo in group:part+=f'  · linear_combination {combo}\n'
  write(f'ReconstructionBlock{bi}',part+end)
 assert all_used==set(equation_block), 'Every retained necessary equation must remain used'
 signature=monolithic[monolithic.index('theorem invariant_polynomial'):monolithic.index('  obtain ⟨')]
 result=header('\n'.join(f'import S11D5Invariants.ReconstructionBlock{j}' for j in range(len(groups))),
  'Assemble the complete coefficient reconstruction. The full-group sufficiency\nproof is separate; the finite certificate is independent of native Q9.',3200000)
 result+=signature
 for bi,group in enumerate(groups):
  names=','.join(f'h{k}' for k,_,_ in group)
  result+=f'  obtain ⟨{names}⟩ := reconstructedBlock{bi} hQ c hc\n'
 result+=monolithic[monolithic.index('  intro G\n'):]
 write('Constraints',result+end)
 return steps,groups

coords=vec([f'G {i} {j}' for i in range(5) for j in range(5)])
frames=[mat([int(j==i) for j in range(25)]) for i in range(25)]
poly=' +\n    '.join(f'c {k} * v {i} * v {j}' for k,(i,j) in enumerate(PAIRS))
frameproof='\n'.join(f'  · change G {k//5} {k%5} = '+ ' + '.join(f"G {j//5} {j%5} * {int(j==k)}" for j in range(25))+'\n    ring' for k in range(25))
coeff=vec([f'A {i} {j}' if i==j else f'A {i} {j} + A {j} {i}' for i,j in PAIRS])
abstract_expansion=' + '.join(f'v {i} • e {i}' for i in range(25))
def right_sum(xs):
 xs=list(xs)
 return xs[0] if len(xs)==1 else xs[0]+' + ('+right_sum(xs[1:])+')'
def tail_vector(i):return right_sum(f'v {j} • e {j}' for j in range(i,25))
def pair_term(i,j):
 a=f'b (e {i}) (e {j})' if i==j else f'(b (e {i}) (e {j}) + b (e {j}) (e {i}))'
 return f'{a} * v {i} * v {j}'
bilinear_steps='''-- Only the small mixed rows are normalized by ring. The complete 325-term
-- identity is assembled with already checked rows and associativity of addition.
private theorem bilinear_add (b : Mat →ₗ[ℝ] Mat →ₗ[ℝ] ℝ) (x y : Mat) :
    b (x + y) (x + y) = b x x + (b x y + b y x) + b y y := by
  simp only [map_add, LinearMap.add_apply]
  ring

private theorem bilinear_diagonal (b : Mat →ₗ[ℝ] Mat →ₗ[ℝ] ℝ)
    (x : Mat) (r : ℝ) : b (r • x) (r • x) = b x x * r * r := by
  simp only [map_smul, LinearMap.smul_apply, smul_eq_mul]
  ring

'''
for i in range(24):
 cross=right_sum(pair_term(i,j) for j in range(i+1,25))
 row_simp='map_smul, LinearMap.smul_apply, smul_eq_mul' if i==23 else 'map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul'
 bilinear_steps+=f'''private theorem bilinear_cross{i} (b : Mat →ₗ[ℝ] Mat →ₗ[ℝ] ℝ)
    (e : Fin 25 → Mat) (v : Coordinates) :
    b (v {i} • e {i}) ({tail_vector(i+1)}) +
      b ({tail_vector(i+1)}) (v {i} • e {i}) =
      {cross} := by
  simp only [{row_simp}]
  ring

'''
for i in reversed(range(25)):
 target=right_sum(pair_term(k,j) for k,j in PAIRS if k>=i)
 proof='  exact bilinear_diagonal b (e 24) (v 24)\n' if i==24 else f'''  rw [bilinear_add, bilinear_cross{i} b e v, bilinear_diagonal,
    bilinear_tail{i+1} b e v]
  simp only [add_assoc]
'''
 bilinear_steps+=f'''private theorem bilinear_tail{i} (b : Mat →ₗ[ℝ] Mat →ₗ[ℝ] ℝ)
    (e : Fin 25 → Mat) (v : Coordinates) :
    b ({tail_vector(i)}) ({tail_vector(i)}) =
      {target} := by
{proof}
'''
text=header('import Mathlib.LinearAlgebra.QuadraticForm.Basic\nimport Mathlib.LinearAlgebra.Dimension.Constructions\nimport Mathlib.LinearAlgebra.Matrix.Trace\nimport Mathlib.Data.List.Fold\nimport Mathlib.Tactic','Exhaustive quadratic forms on all real 5x5 matrices. Generated finite\ncoordinate algebra, not native output transcription; see S11_lean_d5_generate.py.')
text+=f'''abbrev Mat := Matrix (Fin 5) (Fin 5) ℝ
abbrev Quad := QuadraticForm ℝ Mat
abbrev Coordinates := Fin 25 → ℝ
abbrev Coefficients := Fin 325 → ℝ

def coordinates (G : Mat) : Coordinates := {coords}
def decode (v : Coordinates) : Mat := {mat([f'v {i}' for i in range(25)])}

theorem coordinates_decode (v : Coordinates) : coordinates (decode v) = v := by
  ext i
  fin_cases i <;> simp [coordinates, decode]

theorem decode_coordinates (G : Mat) : decode (coordinates G) = G := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [coordinates, decode]

def coordinateMap (i : Fin 25) : Mat →ₗ[ℝ] ℝ where
  toFun G := coordinates G i
  map_add' G H := by fin_cases i <;> simp [coordinates]
  map_smul' r G := by fin_cases i <;> simp [coordinates]

def frame : Fin 25 → Mat := {vec(frames)}

theorem frame_expansion (G : Mat) :
    G = {' + '.join(f'coordinates G {i} • frame {i}' for i in range(25))} := by
  ext i j
  fin_cases i <;> fin_cases j
{frameproof}

def polynomial (c : Coefficients) (v : Coordinates) : ℝ :=
    {poly}

def monomial (i j : Fin 25) : Quad :=
  QuadraticMap.linMulLin (coordinateMap i) (coordinateMap j)

def bilinearCoefficients (A : Fin 25 → Fin 25 → ℝ) : Coefficients := {coeff}

{bilinear_steps}
-- Assemble the checked rows without a whole-polynomial ring certificate.
theorem bilinear_polynomial (b : Mat →ₗ[ℝ] Mat →ₗ[ℝ] ℝ)
    (e : Fin 25 → Mat) (v : Coordinates) :
    b ({abstract_expansion}) ({abstract_expansion}) =
      polynomial (bilinearCoefficients (fun i j => b (e i) (e j))) v := by
  change b ({abstract_expansion}) ({abstract_expansion}) =
    {' + '.join(pair_term(i,j) for i,j in PAIRS)}
  have hexp : ({abstract_expansion}) = ({tail_vector(0)}) :=
    List.foldl1_eq_foldr1 (f := (· + ·))
      (l := [{', '.join(f'v {i} • e {i}' for i in range(1,24))}])
      (a := v 0 • e 0) (b := v 24 • e 24)
  rw [hexp, bilinear_tail0 b e v]
  exact (List.foldl1_eq_foldr1 (f := (· + ·))
    (l := [{', '.join(pair_term(i,j) for i,j in PAIRS[1:-1])}])
    (a := {pair_term(0,0)}) (b := {pair_term(24,24)})).symm

theorem quadratic_representation (Q : Quad) :
    ∃ c : Coefficients, ∀ G, Q G = polynomial c (coordinates G) := by
  let b := Q.associated
  refine ⟨bilinearCoefficients (fun i j => b (frame i) (frame j)), fun G => ?_⟩
  have hb : Q G = b G G := (Q.associated_eq_self_apply ℝ G).symm
  rw [hb, congrArg (fun Z => b Z Z) (frame_expansion G)]
  exact bilinear_polynomial b frame (coordinates G)
'''
# Reset elaboration/kernel working memory at module boundaries. These pieces
# preserve the same definitions and proof text; only the final tail lemma
# becomes visible so the representation module can apply it across an import.
coordinate_text,representation_text=text.split(bilinear_steps,1)
write('CoordinateAlgebra',coordinate_text+end)
write('BilinearExpansion',header('import S11D5Invariants.CoordinateAlgebra',
 'Bounded bilinear row and tail identities for the complete quadratic representation.')+
 bilinear_steps.replace('private theorem bilinear_tail0 ', 'theorem bilinear_tail0 ')+end)
write('Quadratic',header('import S11D5Invariants.BilinearExpansion',
 'Complete quadratic representation assembled from checked coordinate and bilinear identities.')+
 representation_text+end)

# Rotation and Forms are maintained as compact handwritten modules.
x=s.symbols('x0:25');G=s.Matrix(5,5,x)
a,b=s.symbols('a b')
rotations=[]
for p,q in [(0,1),(1,2),(2,3),(3,4)]:
 R=s.eye(5);R[p,p]=a;R[q,q]=a;R[p,q]=-b;R[q,p]=b;rotations.append(R)
names=['rotationXY','rotationYZ','rotationZW','rotationWV']

# Necessary equations from five orthogonal rotations. Preproved coordinate images
# avoid proving the same matrix multiplication at every sample.
Rs=[R.subs({a:0,b:1}) for R in rotations]+[rotations[0].subs({a:s.Rational(3,5),b:s.Rational(4,5)})]
rnames=['rotationXY 0 1','rotationYZ 0 1','rotationZW 0 1','rotationWV 0 1','rotationXY (3/5) (4/5)']
rows=[];samples=[]
for ri,R in enumerate(Rs):
 assert R.T*R==s.eye(5) and R.det()==1
 for i,j in PAIRS:
  v=[0]*25;v[i]=1;v[j]=1
  w=list(R*s.Matrix(5,5,v)*R.T)
  rows.append([w[p]*w[q]-v[p]*v[q] for p,q in PAIRS]);samples.append((ri,v,w))
indices=s.Matrix(rows).T.rref()[1];A=s.Matrix([rows[i] for i in indices])
assert len(indices)==322
forms=[s.trace(G)**2,s.trace(G*G),s.trace(G*G.T)]
B=s.Matrix([[s.Poly(s.expand(f),*x).coeff_monomial(x[i]*x[j]) for i,j in PAIRS] for f in forms])
assert B.rank()==3 and A*B.T==s.zeros(322,3)
free=[PAIRS.index((0,6)),PAIRS.index((1,5)),PAIRS.index((1,1))]
T=s.zeros(3,325)
for i,f in enumerate(free):T[i,f]=s.Rational(1,2) if i<2 else 1
assert T*B.T==s.eye(3)
D=s.eye(325)-B.T*T
# Solve all targets together, once, instead of 325 repeated reductions.
W,params=A.T.gauss_jordan_solve(D.T)
assert params.rows==0
coeffs=[f'(c {free[0]}/2)',f'(c {free[1]}/2)',f'c {free[2]}']
# Split the certificate into bounded builds. Each block proves a
# conjunction of necessary equations; reconstruction below retains all of them.
text=header('import S11D5Invariants.Forms','Coordinate polynomial expansion used by the finite necessary constraints.')
args=' '.join('a'+str(i) for i in range(25))
expanded=' + '.join(f'c {k} * a{i} * a{j}' for k,(i,j) in enumerate(PAIRS))
text+=f'theorem polynomial_vec (c : Coefficients) ({args} : ℝ) :\n    polynomial c {vec(["a"+str(i) for i in range(25)])} = {expanded} := by rfl\n'
write('ConstraintPolynomial',text+end)
blocks=[list(range(i,min(i+20,len(indices)))) for i in range(0,len(indices),20)]
for bi,block in enumerate(blocks):
 text=header('import S11D5Invariants.ConstraintPolynomial','A bounded block of necessary equations from explicit orthogonal rotations.\nThe rational certificate is split only to limit compilation resources.',3200000)
 text+=f'theorem necessaryBlock{bi} {{Q : Quad}} (hQ : OInvariant Q) (c : Coefficients)\n    (hc : ∀ G, Q G = polynomial c (coordinates G)) :\n'
 equations=[]
 for e in block:
  expr=' + '.join(f'({rat(z)}) * c {k}' for k,z in enumerate(A.row(e)) if z)
  equations.append('('+expr+' = 0)')
 text+='    '+' ∧\n    '.join(equations)+' := by\n'
 text+='  refine ⟨'+','.join('?_ ' for _ in block).replace(' ','')+'⟩\n'
 for e in block:
  ri,v,w=samples[indices[e]];rname=names[ri] if ri<4 else names[0]
  text+=f'  · have raw := hQ ({rnames[ri]}) ({rname}_orthogonal (by norm_num)) (decode {vec(v)})\n'
  text+=f'    rw [hc, hc, coordinates_{rname}, coordinates_decode] at raw\n'
  text+='    simp only [decode, Matrix.cons_val] at raw\n'
  text+='    rw [polynomial_vec, polynomial_vec] at raw\n    norm_num at raw\n    linarith only [raw]\n'
 write(f'ConstraintBlock{bi}',text+end)
text=header('\n'.join(f'import S11D5Invariants.ConstraintBlock{bi}' for bi in range(len(blocks))),
 'Reconstruction from every necessary equation. The full-group sufficiency\nproof is separate; the finite certificate is independent of native Q9.',3200000)
text+='theorem invariant_polynomial {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)\n    (hc : ∀ G, Q G = polynomial c (coordinates G)) :\n    ∀ G, Q G = '+ ' +\n      '.join(f'{q} * {f}' for q,f in zip(coeffs,['G.trace ^ 2','(G*G).trace','(G*G.transpose).trace']))+' := by\n'
for bi,block in enumerate(blocks):
 text+='  obtain ⟨'+','.join(f'e{e}' for e in block)+f'⟩ := necessaryBlock{bi} hQ c hc\n'
for k in range(325):
 if k in free:continue
 target=' + '.join(f'{rat(B[r,k])} * {coeffs[r]}' for r in range(3) if B[r,k]) or '0'
 combo=' + '.join(f'({rat(z)}) * e{i}' for i,z in enumerate(W.col(k)) if z)
 text+=f'  have h{k} : c {k} = {target} := by\n    linear_combination {combo}\n'
text+='  intro G\n  rw [hc]\n  change polynomial c '+coords+' = _\n  rw [polynomial_vec]\n  simp only ['+', '.join(f'h{i}' for i in range(325) if i not in free)+']\n'
text+='  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply, sum_five]\n  ring\n'
split_reconstruction(text,blocks,free)
print('Generated complete 325-coefficient representation and 322 necessary constraints; free coordinates',free)
