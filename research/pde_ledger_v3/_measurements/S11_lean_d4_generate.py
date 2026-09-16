#!/usr/bin/env python3
"""D4 finite quadratic completeness certificates, independent of native Q9.

SymPy selects rational necessary constraints. Lean checks the exhaustive
representation, all retained equations, reconstruction and full-group
sufficiency. This generator is not a trusted mathematical oracle.
"""
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
import sys

OUT=Path(__file__).resolve().parents[1]/'lean/s11/S11D4Invariants'
PAIRS=list(combinations_with_replacement(range(16),2))
def vec(xs):return '!['+','.join(map(str,xs))+']'
def mat(xs):return '!['+','.join(vec(xs[4*i:4*i+4]) for i in range(4))+']'
def rat(x):
 x=s.Rational(x)
 return str(x.p) if x.q==1 else f'({x.p}/{x.q})'
def write(name,text):
 path=OUT/(name+'.lean')
 if '--check' in sys.argv:assert path.read_text()==text,f'Generated source drift: {path}'
 else:path.write_text(text)

def header(imports,doc,heartbeats=1600000):return imports+'\n\n/-! '+doc+' -/\nnamespace S11D4Invariants\nnoncomputable section\nset_option maxRecDepth 4096\nset_option maxHeartbeats '+str(heartbeats)+'\n\n'
end='\nend\nend S11D4Invariants\n'
coords=vec([f'G {i} {j}' for i in range(4) for j in range(4)])
frames=[mat([int(j==i) for j in range(16)]) for i in range(16)]
poly=' +\n    '.join(f'c {k} * v {i} * v {j}' for k,(i,j) in enumerate(PAIRS))
frameproof='\n'.join(f'  · change G {k//4} {k%4} = '+ ' + '.join(f"G {j//4} {j%4} * {int(j==k)}" for j in range(16))+'\n    ring' for k in range(16))
coeff=vec([f'b (frame {i}) (frame {j})' if i==j else f'b (frame {i}) (frame {j}) + b (frame {j}) (frame {i})' for i,j in PAIRS])
text=header('import Mathlib.LinearAlgebra.QuadraticForm.Basic\nimport Mathlib.LinearAlgebra.Dimension.Constructions\nimport Mathlib.LinearAlgebra.Matrix.Trace\nimport Mathlib.Tactic','Exhaustive quadratic forms on all real 4x4 matrices. Generated finite\ncoordinate algebra, not native output transcription; see S11_lean_d4_generate.py.')
text+=f'''abbrev Mat := Matrix (Fin 4) (Fin 4) ℝ
abbrev Quad := QuadraticForm ℝ Mat
abbrev Coordinates := Fin 16 → ℝ
abbrev Coefficients := Fin 136 → ℝ

def coordinates (G : Mat) : Coordinates := {coords}
def decode (v : Coordinates) : Mat := {mat([f'v {i}' for i in range(16)])}

theorem coordinates_decode (v : Coordinates) : coordinates (decode v) = v := by
  ext i
  fin_cases i <;> simp [coordinates, decode]

theorem decode_coordinates (G : Mat) : decode (coordinates G) = G := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [coordinates, decode]

def coordinateMap (i : Fin 16) : Mat →ₗ[ℝ] ℝ where
  toFun G := coordinates G i
  map_add' G H := by fin_cases i <;> simp [coordinates]
  map_smul' r G := by fin_cases i <;> simp [coordinates]

def frame : Fin 16 → Mat := {vec(frames)}

theorem frame_expansion (G : Mat) :
    G = {' + '.join(f'coordinates G {i} • frame {i}' for i in range(16))} := by
  ext i j
  fin_cases i <;> fin_cases j
{frameproof}

def polynomial (c : Coefficients) (v : Coordinates) : ℝ :=
    {poly}

def monomial (i j : Fin 16) : Quad :=
  QuadraticMap.linMulLin (coordinateMap i) (coordinateMap j)

theorem quadratic_representation (Q : Quad) :
    ∃ c : Coefficients, ∀ G, Q G = polynomial c (coordinates G) := by
  let b := Q.associated
  let c : Coefficients := {coeff}
  refine ⟨c, fun G => ?_⟩
  have hb : Q G = b G G := (Q.associated_eq_self_apply ℝ G).symm
  rw [hb, congrArg (fun Z => b Z Z) (frame_expansion G)]
  simp only [polynomial, c, Matrix.cons_val]
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul]
  ring
'''
write('Quadratic',text+end)

# Compact determinant expression and three actual plane rotations.
x=s.symbols('x0:16');G=s.Matrix(4,4,x)
a,b=s.symbols('a b')
rotations=[]
for p,q in [(0,1),(1,2),(2,3)]:
 R=s.eye(4);R[p,p]=a;R[q,q]=a;R[p,q]=-b;R[q,p]=b;rotations.append(R)
names=['rotationXY','rotationYZ','rotationZW']
def leanexpr(f,prefix='G'):
 t=str(s.expand(f)).replace('**','^')
 import re
 return re.sub(r'\bx(\d+)\b',lambda m:f'({prefix} {int(m[1])//4} {int(m[1])%4})',t)
text=header('import S11D4Invariants.Quadratic','Full group predicates, elementary rotations and trace identities.\nGenerated coordinate identities are checked in Lean.')
text+='''def conjugate (R G : Mat) : Mat := R * G * R.transpose
def Orthogonal (R : Mat) : Prop := R.transpose * R = 1
def Proper (R : Mat) : Prop := Orthogonal R ∧ R.det = 1

def SOInvariant (Q : Quad) : Prop := ∀ R, Proper R → ∀ G, Q (conjugate R G) = Q G
def OInvariant (Q : Quad) : Prop := ∀ R, Orthogonal R → ∀ G, Q (conjugate R G) = Q G

def reflection : Mat := ![![-1,0,0,0],![0,1,0,0],![0,0,1,0],![0,0,0,1]]
def ReflectionOdd (Q : Quad) : Prop := ∀ G, Q (conjugate reflection G) = -Q G

'''
text+=f'''theorem det_four (R : Mat) : R.det = {leanexpr(G.det(),'R')} := by
  rw [Matrix.det_succ_row_zero]
  simp [Fin.sum_univ_four, Matrix.det_fin_three, Matrix.submatrix_apply, Fin.succAbove,
    Fin.succ, Fin.castSucc]
  ring

'''
for name,R in zip(names,rotations):
 text+=f'''def {name} (a b : ℝ) : Mat := {mat(list(R))}

theorem {name}_proper {{a b : ℝ}} (h : a^2+b^2=1) : Proper ({name} a b) := by
  constructor
  · unfold Orthogonal
    ext i j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_four]
    fin_cases i <;> fin_cases j <;> simp [{name}] <;> nlinarith
  · rw [det_four]
    simp [{name}]
    nlinarith

theorem coordinates_{name} (a b : ℝ) (G : Mat) :
    coordinates (conjugate ({name} a b) G) = {vec([leanexpr(v) for v in R*G*R.T])} := by
  ext i
  fin_cases i <;>
    simp only [coordinates, conjugate, Matrix.mul_apply,
      Matrix.transpose_apply, Fin.sum_univ_four] <;>
    norm_num [{name}, Matrix.cons_val_two, Matrix.cons_val_three] <;> ring

'''
text+='''theorem reflection_orthogonal : Orthogonal reflection := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_four]
  fin_cases i <;> fin_cases j <;> norm_num [reflection, Matrix.cons_val_two, Matrix.cons_val_three]

theorem reflection_det : reflection.det = -1 := by
  rw [det_four]
  norm_num [reflection, Matrix.cons_val_two, Matrix.cons_val_three]

theorem trace_conjugate {R : Mat} (hR : Orthogonal R) (G : Mat) :
    (conjugate R G).trace = G.trace := by
  unfold conjugate
  rw [Matrix.trace_mul_cycle, hR, Matrix.one_mul]

theorem conjugate_mul {R : Mat} (hR : Orthogonal R) (G H : Mat) :
    conjugate R G * conjugate R H = conjugate R (G*H) := by
  unfold conjugate
  calc R*G*R.transpose*(R*H*R.transpose) = R*G*(R.transpose*R)*H*R.transpose := by
         simp only [Matrix.mul_assoc]
       _ = R*(G*H)*R.transpose := by
         rw [hR, Matrix.mul_one]
         simp only [Matrix.mul_assoc]

theorem conjugate_transpose (R G : Mat) :
    (conjugate R G).transpose = conjugate R G.transpose := by
  simp only [conjugate, Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.mul_assoc]
'''
write('Rotation',text+end)

# Necessary equations from four proper rotations. Preproved coordinate images
# avoid proving the same matrix multiplication at every sample.
Rs=[R.subs({a:0,b:1}) for R in rotations]+[rotations[0].subs({a:s.Rational(3,5),b:s.Rational(4,5)})]
rnames=['rotationXY 0 1','rotationYZ 0 1','rotationZW 0 1','rotationXY (3/5) (4/5)']
rows=[];samples=[]
for ri,R in enumerate(Rs):
 assert R.T*R==s.eye(4) and R.det()==1
 for i,j in PAIRS:
  v=[0]*16;v[i]=1;v[j]=1
  w=list(R*s.Matrix(4,4,v)*R.T)
  rows.append([w[p]*w[q]-v[p]*v[q] for p,q in PAIRS]);samples.append((ri,v,w))
indices=s.Matrix(rows).T.rref()[1];A=s.Matrix([rows[i] for i in indices])
assert len(indices)==132
P=(G[0,1]-G[1,0])*(G[2,3]-G[3,2])-(G[0,2]-G[2,0])*(G[1,3]-G[3,1])+(G[0,3]-G[3,0])*(G[1,2]-G[2,1])
forms=[s.trace(G)**2,s.trace(G*G),s.trace(G*G.T),P]
B=s.Matrix([[s.Poly(s.expand(f),*x).coeff_monomial(x[i]*x[j]) for i,j in PAIRS] for f in forms])
assert B.rank()==4 and A*B.T==s.zeros(132,4)
free=[PAIRS.index((0,5)),PAIRS.index((1,4)),PAIRS.index((1,1)),PAIRS.index((1,11))]
T=s.zeros(4,136)
for i,f in enumerate(free):T[i,f]=s.Rational(1,2) if i<2 else 1
assert T*B.T==s.eye(4)
D=s.eye(136)-B.T*T
# Solve all targets together, once, instead of 136 repeated reductions.
W,params=A.T.gauss_jordan_solve(D.T)
assert params.rows==0
coeffs=[f'(c {free[0]}/2)',f'(c {free[1]}/2)',f'c {free[2]}',f'c {free[3]}']
# Split the unchanged certificate into bounded builds. Each block proves a
# conjunction of necessary equations; reconstruction below retains all of them.
text=header('import S11D4Invariants.Forms','Coordinate polynomial expansion used by the finite necessary constraints.')
args=' '.join('a'+str(i) for i in range(16))
expanded=' + '.join(f'c {k} * a{i} * a{j}' for k,(i,j) in enumerate(PAIRS))
text+=f'theorem polynomial_vec (c : Coefficients) ({args} : ℝ) :\n    polynomial c {vec(["a"+str(i) for i in range(16)])} = {expanded} := by rfl\n'
write('ConstraintPolynomial',text+end)
blocks=[list(range(i,min(i+33,len(indices)))) for i in range(0,len(indices),33)]
for bi,block in enumerate(blocks):
 text=header('import S11D4Invariants.ConstraintPolynomial','A bounded block of necessary equations from explicit proper rotations.\nThis is the same rational certificate, split only to limit compilation resources.',3200000)
 text+=f'theorem necessaryBlock{bi} {{Q : Quad}} (hQ : SOInvariant Q) (c : Coefficients)\n    (hc : ∀ G, Q G = polynomial c (coordinates G)) :\n'
 equations=[]
 for e in block:
  expr=' + '.join(f'({rat(z)}) * c {k}' for k,z in enumerate(A.row(e)) if z)
  equations.append('('+expr+' = 0)')
 text+='    '+' ∧\n    '.join(equations)+' := by\n'
 text+='  refine ⟨'+','.join('?_ ' for _ in block).replace(' ','')+'⟩\n'
 for e in block:
  ri,v,w=samples[indices[e]];rname=names[ri] if ri<3 else names[0]
  text+=f'  · have raw := hQ ({rnames[ri]}) ({rname}_proper (by norm_num)) (decode {vec(v)})\n'
  text+=f'    rw [hc, hc, coordinates_{rname}, coordinates_decode] at raw\n'
  text+='    simp only [decode, Matrix.cons_val] at raw\n'
  text+='    rw [polynomial_vec, polynomial_vec] at raw\n    norm_num at raw\n    linarith only [raw]\n'
 write(f'ConstraintBlock{bi}',text+end)
text=header('\n'.join(f'import S11D4Invariants.ConstraintBlock{bi}' for bi in range(len(blocks))),
 'Reconstruction from every necessary equation. The full-group sufficiency\nproof is separate; the finite certificate is independent of native Q9.',3200000)
text+='theorem invariant_polynomial {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)\n    (hc : ∀ G, Q G = polynomial c (coordinates G)) :\n    ∀ G, Q G = '+ ' +\n      '.join(f'{q} * {f}' for q,f in zip(coeffs,['G.trace ^ 2','(G*G).trace','(G*G.transpose).trace','orientation G']))+' := by\n'
for bi,block in enumerate(blocks):
 text+='  obtain ⟨'+','.join(f'e{e}' for e in block)+f'⟩ := necessaryBlock{bi} hQ c hc\n'
for k in range(136):
 if k in free:continue
 target=' + '.join(f'{rat(B[r,k])} * {coeffs[r]}' for r in range(4) if B[r,k]) or '0'
 combo=' + '.join(f'({rat(z)}) * e{i}' for i,z in enumerate(W.col(k)) if z)
 text+=f'  have h{k} : c {k} = {target} := by\n    linear_combination {combo}\n'
text+='  intro G\n  rw [hc]\n  change polynomial c '+coords+' = _\n  rw [polynomial_vec]\n  simp only ['+', '.join(f'h{i}' for i in range(136) if i not in free)+']\n'
text+='  simp only [orientation, Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_four]\n  ring\n'
write('Constraints',text+end)
print('Generated complete 136-coefficient representation and 132 necessary constraints; free coordinates',free)
