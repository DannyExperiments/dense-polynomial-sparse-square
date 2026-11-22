import Std

namespace SparseSquare

/-- The open middle third of the residue interval modulo p. -/
def middleThird (p a : Nat) : Prop := p < 3 * a ∧ 3 * a < 2 * p

/-- Adding two elements of the forbidden interval never lands back in it. -/
theorem middleThird_sumfree {p a b : Nat} (ha : middleThird p a)
    (hb : middleThird p b) : ¬ middleThird p ((a + b) % p) := by
  unfold middleThird at *
  by_cases h : a + b < p
  · rw [Nat.mod_eq_of_lt h]
    omega
  · have hp : p ≤ a + b := by omega
    have hrem : a + b - p < p := by omega
    rw [Nat.mod_eq_sub_mod hp, Nat.mod_eq_of_lt hrem]
    omega

/-- Evaluate an integer polynomial, storing the constant coefficient first. -/
def evaluate (z : Int) : List Int → Int
  | [] => 0
  | c :: cs => c + z * evaluate z cs

/-- At a large integer argument, a polynomial with small signed integer
coefficients can vanish only when all its coefficients vanish. -/
theorem coefficients_zero_of_evaluate_zero {z : Int} (hz : z ≠ 0)
    (cs : List Int) (bound : ∀ c ∈ cs, c.natAbs < z.natAbs)
    (heval : evaluate z cs = 0) : ∀ c ∈ cs, c = 0 := by
  induction cs with
  | nil => simp
  | cons a cs ih =>
    have h : a + z * evaluate z cs = 0 := heval
    have hdvd : z ∣ a := by
      refine ⟨-evaluate z cs, ?_⟩
      rw [Int.mul_neg]
      omega
    have ha := Int.eq_zero_of_dvd_of_natAbs_lt_natAbs hdvd (bound a (by simp))
    have hprod : z * evaluate z cs = 0 := by omega
    have htail := (Int.mul_eq_zero.mp hprod).resolve_left hz
    have hall := ih (fun c hc => bound c (by simp [hc])) htail
    intro c hc
    rcases List.mem_cons.mp hc with hca | hcs
    · exact hca.trans ha
    · exact hall c hcs

/-- The coefficient of the formal scaling parameter at power j in the
coefficient of x^i of E(x)V(zx). Terms with j > i contribute zero. -/
def diagonalTerm (e v : Nat → Int) (i j : Nat) : Int :=
  if j ≤ i then e (i - j) * v j else 0

/-- A dense filler V of degree t times an endpoint-nonzero E of degree d ≤ t
has no identically zero coefficient when V is formally scaled. -/
theorem diagonal_nonzero_term {d t i : Nat} (e v : Nat → Int)
    (hdt : d ≤ t) (hi : i ≤ t + d) (h₀ : e 0 ≠ 0) (hd : e d ≠ 0)
    (hv : ∀ j ≤ t, v j ≠ 0) :
    ∃ j, j ≤ t ∧ diagonalTerm e v i j ≠ 0 := by
  by_cases hit : i ≤ t
  · exact ⟨i, hit, by simpa [diagonalTerm] using Int.mul_ne_zero h₀ (hv i hit)⟩
  · have hj : i - d ≤ t := by omega
    have heq : i - (i - d) = d := by omega
    exact ⟨i - d, hj, by simpa [diagonalTerm, heq] using Int.mul_ne_zero hd (hv _ hj)⟩

/-- Ordinary coefficient convolution for E(x)V(zx), with V cut off at t. -/
def scaledProductCoeff (z : Int) (e v : Nat → Int) (t i : Nat) : Int :=
  evaluate z ((List.range (t + 1)).map (diagonalTerm e v i))

/-- Every coefficient of the scaled product through degree t+d stays nonzero
at any integer scaling parameter exceeding the absolute values of the terms. -/
theorem scaledProductCoeff_ne_zero {z : Int} {d t : Nat} (e v : Nat → Int)
    (hz : z ≠ 0) (hdt : d ≤ t) (h₀ : e 0 ≠ 0) (hd : e d ≠ 0)
    (hv : ∀ j ≤ t, v j ≠ 0) (i : Nat) (hi : i ≤ t + d)
    (bound : ∀ j ≤ t, (diagonalTerm e v i j).natAbs < z.natAbs) :
    scaledProductCoeff z e v t i ≠ 0 := by
  obtain ⟨j, hj, hterm⟩ := diagonal_nonzero_term e v hdt hi h₀ hd hv
  intro heval
  have hall := coefficients_zero_of_evaluate_zero hz _ (by
    intro c hc
    obtain ⟨k, hk, rfl⟩ := List.mem_map.mp hc
    exact bound k (by have := List.mem_range.mp hk; omega)) heval
  exact hterm (hall _ (List.mem_map.mpr ⟨j, List.mem_range.mpr (by omega), rfl⟩))

end SparseSquare
