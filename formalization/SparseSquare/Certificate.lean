import SparseSquare.Data
import Std

namespace SparseSquare.Certificate
open SparseSquare.Data

def localBound : Nat := (18 * periods.length + 1) * period * 7 ^ periods.length /
  (10 ^ periods.length)

def repairedBound : Nat :=
  5 * localBound + 18 * auxiliaryBase * 10 ^ (6 * periods.length) +
    3 * (8 * auxiliaryBase + 2) * 9 ^ periods.length

def endDegree : Nat :=
  period + 4 * (periods.map (fun p => period / p)).foldl Nat.add 0 +
    auxiliaryBase * (auxiliaryBase + 1)

set_option maxRecDepth 20000
set_option maxHeartbeats 2000000

theorem period_count : periods.length = 128 := by decide +kernel

theorem period_product : periods.foldl Nat.mul 1 = period := by decide +kernel

theorem auxiliary_base_ceiling :
    (auxiliaryBase - 1) ^ 2 < (14 * periods.length + 3) * period ∧
    (14 * periods.length + 3) * period ≤ auxiliaryBase ^ 2 := by decide +kernel

theorem repetition_count : repetitions = period / 400 := by decide +kernel

theorem repair_weight_bound : 400 * repairedBound < period := by decide +kernel

theorem endpoint_separation :
    9 * periods.length * period + 2 ≤ endDegree ∧ endDegree < repetitions * period := by
  decide +kernel

theorem recorded_degree : degree = repetitions * period + endDegree := by decide +kernel

theorem recorded_weight_upper :
    weightUpper = 5 * (2 * repetitions - 1) + 3 * repairedBound := by decide +kernel

theorem degree_size : 10 ^ 5797 ≤ degree ∧ degree < 10 ^ 5798 := by decide +kernel

theorem strict_ratio : (1000 * weightUpper) ^ 2 < 501 ^ 2 * degree := by decide +kernel

theorem original_threshold : weightUpper ^ 2 < degree := by decide +kernel

/-- Apply the numerical certificates to any proved weight bound for the witness. -/
theorem transfer_to_weight (w : Nat) (hw : w ≤ weightUpper) :
    (1000 * w) ^ 2 < 501 ^ 2 * degree ∧ w ^ 2 < degree := by
  constructor
  · exact Nat.lt_of_le_of_lt (Nat.pow_le_pow_left (Nat.mul_le_mul_left 1000 hw) 2) strict_ratio
  · exact Nat.lt_of_le_of_lt (Nat.pow_le_pow_left hw 2) original_threshold

end SparseSquare.Certificate
