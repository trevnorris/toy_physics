#!/usr/bin/env python3
"""Finite-coordinate proof text for J1, independent of native CAS audit outputs.

Generate the exhaustive quadratic representation and a rational linear
certificate from three explicit proper rotations. All generated mathematical
steps must be kernel checked; SymPy selects certificates but is not trusted.
No production source, export, or previous formalization is read or changed.
"""
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
import sys

OUT=Path(__file__).resolve().parents[1]/'lean/s11/S11D3Invariants'
PAIRS=list(combinations_with_replacement(range(9),2))
def vec(xs): return '!['+','.join(str(x) for x in xs)+']'
def mat(xs): return '!['+','.join(vec(xs[3*i:3*i+3]) for i in range(3))+']'
def rat(x):
 x=s.Rational(x)
 return str(x.p) if x.q==1 else f'({x.p}/{x.q})'
def write(name,text):
 path=OUT/(name+'.lean')
 if '--check' in sys.argv: assert path.read_text()==text, f'Generated source drift: {path}'
 else: path.write_text(text)

header='''import Mathlib.LinearAlgebra.QuadraticForm.Basic
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic

/-! Exhaustive quadratic forms on all real 3x3 matrices. Generated finite
coordinate algebra, not CAS-output transcription; see S11_lean_d3_generate.py. -/
namespace S11D3Invariants
noncomputable section
set_option maxRecDepth 2048
set_option maxHeartbeats 1600000

abbrev Mat := Matrix (Fin 3) (Fin 3) ℝ
abbrev Quad := QuadraticForm ℝ Mat
abbrev Coordinates := Fin 9 → ℝ
abbrev Coefficients := Fin 45 → ℝ

'''
coords=vec([f'G {i} {j}' for i in range(3) for j in range(3)])
decode=mat([f'v {i}' for i in range(9)])
frames=[mat([int(j==i) for j in range(9)]) for i in range(9)]
poly=' +\n    '.join(f'c {k} * v {i} * v {j}' for k,(i,j) in enumerate(PAIRS))
frameproof='\n'.join(f'  · change G {k//3} {k%3} = '+ ' + '.join(f"G {j//3} {j%3} * {int(j==k)}" for j in range(9))+'\n    ring' for k in range(9))
coeff=vec([f'b (frame {i}) (frame {j})' if i==j else f'b (frame {i}) (frame {j}) + b (frame {j}) (frame {i})' for i,j in PAIRS])
write('Quadratic',header+f'''def coordinates (G : Mat) : Coordinates := {coords}
def decode (v : Coordinates) : Mat := {decode}

theorem coordinates_decode (v : Coordinates) : coordinates (decode v) = v := by
  ext i
  fin_cases i <;> simp [coordinates, decode]

theorem decode_coordinates (G : Mat) : decode (coordinates G) = G := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [coordinates, decode]

def coordinateMap (i : Fin 9) : Mat →ₗ[ℝ] ℝ where
  toFun G := coordinates G i
  map_add' G H := by fin_cases i <;> simp [coordinates]
  map_smul' r G := by fin_cases i <;> simp [coordinates]

def frame : Fin 9 → Mat := {vec(frames)}

theorem frame_expansion (G : Mat) :
    G = {' + '.join(f'coordinates G {i} • frame {i}' for i in range(9))} := by
  ext i j
  fin_cases i <;> fin_cases j
{frameproof}

def polynomial (c : Coefficients) (v : Coordinates) : ℝ :=
    {poly}

def monomial (i j : Fin 9) : Quad :=
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

end
end S11D3Invariants
''')

# For each explicit SO rotation, evaluate invariance at coordinate basis vectors
# and sums of two distinct basis vectors. Select independent necessary equations.
rotations=[s.Matrix([[0,-1,0],[1,0,0],[0,0,1]]),s.Matrix([[1,0,0],[0,0,-1],[0,1,0]]),s.Matrix([[s.Rational(3,5),s.Rational(-4,5),0],[s.Rational(4,5),s.Rational(3,5),0],[0,0,1]])]
rotation_names=['rotationXY 0 1','rotationYZ 0 1','rotationXY (3/5) (4/5)']
proper_names=['rotationXY_proper','rotationYZ_proper','rotationXY_proper']
rows=[]; samples=[]
for ri,R in enumerate(rotations):
 assert R.T*R==s.eye(3) and R.det()==1
 for i,j in PAIRS:
  v=[0]*9;v[i]=1;v[j]=1
  w=list(R*s.Matrix(3,3,v)*R.T)
  row=[w[a]*w[b]-v[a]*v[b] for a,b in PAIRS]
  rows.append(row); samples.append((ri,v))
indices=s.Matrix(rows).T.rref()[1]; A=s.Matrix([rows[i] for i in indices])
assert A.rank()==42
# Target parameter extraction: coefficients of g00*g11, g01*g10 and g01^2.
x=s.symbols('x0:9');G=s.Matrix(3,3,x)
forms=[s.trace(G)**2,s.trace(G*G),s.trace(G*G.T)]
B=s.Matrix([[s.Poly(s.expand(f),*x).coeff_monomial(x[i]*x[j]) for i,j in PAIRS] for f in forms])
assert B.rank()==3 and A*B.T==s.zeros(42,3)
free=[PAIRS.index((0,4)),PAIRS.index((1,3)),PAIRS.index((1,1))]
T=s.zeros(3,45);T[0,free[0]]=s.Rational(1,2);T[1,free[1]]=s.Rational(1,2);T[2,free[2]]=1
assert T*B.T==s.eye(3)
D=s.eye(45)-B.T*T
weights=[]
for k in range(45):
 w,params=A.T.gauss_jordan_solve(D.row(k).T)
 assert not params.rows
 weights.append(list(w))
text='''import S11D3Invariants.Rotation

/-! Three explicit admissible rotations give necessary constraints only.
Every retained equation is derived from full SO invariance. Rational linear
combinations are checked by Lean, including the exhaustive representation. -/
namespace S11D3Invariants
noncomputable section
set_option maxRecDepth 2048
set_option maxHeartbeats 3200000

'''
params=' '.join('a'+str(i) for i in range(9))
vector=vec(['a'+str(i) for i in range(9)])
expanded=' + '.join(f'c {k} * a{i} * a{j}' for k,(i,j) in enumerate(PAIRS))
text+=f'theorem polynomial_vec (c : Coefficients) ({params} : ℝ) :\n    polynomial c {vector} = {expanded} := by rfl\n\n'
text+='theorem invariant_polynomial {Q : Quad} (hQ : SOInvariant Q) (c : Coefficients)\n    (hc : ∀ G, Q G = polynomial c (coordinates G)) :\n    ∀ G, Q G = (c '+str(free[0])+'/2) * G.trace ^ 2 +\n      (c '+str(free[1])+'/2) * (G*G).trace + c '+str(free[2])+' * (G*G.transpose).trace := by\n'
for e,idx in enumerate(indices):
 ri,v=samples[idx]
 row=A.row(e)
 expr=' + '.join(f'({rat(a)}) * c {k}' for k,a in enumerate(row) if a)
 text+=f'  have e{e} : {expr} = 0 := by\n'
 w=list(rotations[ri]*s.Matrix(3,3,v)*rotations[ri].T)
 image='(Matrix.of '+mat([rat(a) for a in w])+')';input_matrix='(Matrix.of '+mat(v)+')'
 text+=f'    have image : conjugate ({rotation_names[ri]}) ({input_matrix}) = {image} := by\n'
 text+='      ext i j\n      rw [conjugate_apply]\n      fin_cases i <;> fin_cases j <;>\n'
 text+=f'        norm_num [Fin.sum_univ_three, {"rotationXY" if ri!=1 else "rotationYZ"}, Matrix.cons_val_two]\n'
 text+=f'    have raw := hQ ({rotation_names[ri]}) ({proper_names[ri]} (by norm_num)) ({input_matrix})\n'
 text+='    rw [image, hc, hc] at raw\n'
 text+=f'    change polynomial c {vec([rat(a) for a in w])} = polynomial c {vec(v)} at raw\n'
 text+='    rw [polynomial_vec, polynomial_vec] at raw\n    norm_num at raw\n    linarith only [raw]\n'

for k,w in enumerate(weights):
 target=' + '.join(f'{rat(B[r,k])} * '+([f'(c {free[0]}/2)',f'(c {free[1]}/2)',f'c {free[2]}'][r]) for r in range(3) if B[r,k]) or '0'
 if k in free:continue
 text+=f'  have h{k} : c {k} = {target} := by\n'
 combo=' + '.join(f'({rat(a)}) * e{i}' for i,a in enumerate(w) if a)
 text+=f'    linear_combination {combo}\n'
text+=f'  intro G\n  rw [hc]\n  change polynomial c {coords} = _\n  rw [polynomial_vec]\n  simp only ['+', '.join(f'h{i}' for i in range(45) if i not in free)+']\n'
text+='  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]\n  ring\n\nend\nend S11D3Invariants\n'
write('Constraints',text)
print('Generated exhaustive 45-coefficient representation and 42 necessary constraints; free coordinates',free)
