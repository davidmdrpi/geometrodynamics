import Mathlib.MeasureTheory.Measure.Lebesgue.VolumeOfBalls
import Mathlib.Tactic

/-! A two-port fixed-point boundary model. The setting varies with the same
history and the same preparation measure throughout. -/
noncomputable section
open Set MeasureTheory
namespace BAM.ClosedHistory

structure History where
  right : ℝ
  left : ℝ

/-- Two port visits with rank-reversing transport on both legs. -/
def returnMap (h : History) : History := ⟨1 - h.left, 1 - h.right⟩

def Right (h : History) (p : ℝ) : Prop := h.right < p
def Left (h : History) (p : ℝ) : Prop := h.left ≤ 1 - p
def Operational (h : History) (p : ℝ) : Prop :=
  (¬ (Right h p ∧ Left h p)) ∧ (Right h p ∨ Left h p)

theorem fixed_iff_budget (h : History) :
    returnMap h = h ↔ h.right + h.left = 1 := by
  cases h with
  | mk u v =>
    change (History.mk (1-v) (1-u) = History.mk u v) ↔ u+v=1
    constructor
    · intro he
      have e := congrArg History.right he
      change 1-v=u at e
      linarith
    · intro he
      have a : 1-v=u := by linarith
      have b : 1-u=v := by linarith
      rw [a,b]

theorem fixed_operational (h : History) (hf : returnMap h = h) (p : ℝ) :
    Operational h p := by
  have he := (fixed_iff_budget h).mp hf
  unfold Operational Right Left
  constructor
  · rintro ⟨hr, hl⟩
    linarith
  · by_cases hr : h.right < p
    · exact Or.inl hr
    · exact Or.inr (by linarith)

/-- Between the two switch points the ports either both fire or both miss.
The half-open boundary implements the same deterministic tie rule at every p. -/
theorem failing_settings (h : History) :
    {p : ℝ | ¬ Operational h p} =
      Ioc (min h.right (1-h.left)) (max h.right (1-h.left)) := by
  ext p
  have e : h.left ≤ 1-p ↔ p ≤ 1-h.left := by constructor <;> intro a <;> linarith
  change (¬ ((¬ (h.right < p ∧ h.left ≤ 1-p)) ∧ (h.right < p ∨ h.left ≤ 1-p))) ↔
      min h.right (1-h.left) < p ∧ p ≤ max h.right (1-h.left)
  rw [e]
  rcases le_total h.right (1-h.left) with h₁ | h₂
  · rw [min_eq_left h₁, max_eq_right h₁]
    by_cases a : h.right < p <;> by_cases b : p ≤ 1-h.left <;> simp_all <;> linarith
  · rw [min_eq_right h₂, max_eq_left h₂]
    by_cases a : h.right < p <;> by_cases b : p ≤ 1-h.left <;> simp_all <;> linarith

theorem failure_length (h : History) :
    volume.real {p : ℝ | ¬ Operational h p} = |h.right + h.left - 1| := by
  rw [failing_settings, Real.volume_real_Ioc_of_le (min_le_max)]
  rcases le_total h.right (1-h.left) with a | b
  · rw [max_eq_right a, min_eq_left a, abs_of_nonpos (by linarith : h.right+h.left-1 ≤ 0)]
    ring
  · rw [max_eq_left b, min_eq_right b, abs_of_nonneg (by linarith : 0 ≤ h.right+h.left-1)]
    ring

/-- Physical settings span [0,1]. Necessity uses the midpoint between
the two switch points, so it does not rely on settings outside that range. -/
theorem all_settings_iff_fixed (h : History) (hu : h.right ∈ Icc (0:ℝ) 1)
    (hv : h.left ∈ Icc (0:ℝ) 1) :
    (∀ p ∈ Icc (0:ℝ) 1, Operational h p) ↔ returnMap h = h := by
  constructor
  · intro hall
    apply (fixed_iff_budget h).mpr
    rcases hu with ⟨hu0,hu1⟩
    rcases hv with ⟨hv0,hv1⟩
    let m := (h.right + (1-h.left))/2
    have hm : m ∈ Icc (0:ℝ) 1 := by constructor <;> dsimp [m] <;> linarith
    obtain ⟨hd,hm'⟩ := hall m hm
    change ¬ (h.right < m ∧ h.left ≤ 1-m) at hd
    change h.right < m ∨ h.left ≤ 1-m at hm'
    dsimp [m] at hd hm'
    rcases lt_trichotomy (h.right+h.left) 1 with lo | eqn | hi
    · exact False.elim (hd ⟨by linarith, by linarith⟩)
    · exact eqn
    · rcases hm' with a|a <;> linarith
  · intro hf p _
    exact fixed_operational h hf p

abbrev Plane := EuclideanSpace ℝ (Fin 2)
def disk : Set Plane := Metric.ball 0 1
def history (z : Plane) : History := ⟨‖z‖^2, 1-‖z‖^2⟩

theorem disk_history_fixed (z : Plane) : returnMap (history z) = history z := by
  apply (fixed_iff_budget _).mpr
  change ‖z‖^2+(1-‖z‖^2)=1
  ring

def rightBasin (θ : ℝ) := {z ∈ disk | Right (history z) (Real.cos θ ^ 2)}
def leftBasin (θ : ℝ) := {z ∈ disk | Left (history z) (Real.cos θ ^ 2)}
def fraction (s : Set Plane) : ℝ := volume.real s / volume.real disk

theorem right_disk (θ : ℝ) :
    rightBasin θ = Metric.ball 0 |Real.cos θ| := by
  have hc : |Real.cos θ| ≤ 1 := abs_le.mpr ⟨Real.neg_one_le_cos _, Real.cos_le_one _⟩
  ext z
  change (z ∈ Metric.ball 0 1 ∧ ‖z‖^2 < Real.cos θ ^ 2) ↔ z ∈ Metric.ball 0 |Real.cos θ|
  simp only [Metric.mem_ball, dist_zero_right]
  have e : ‖z‖^2 < Real.cos θ ^ 2 ↔ ‖z‖ < |Real.cos θ| := by
    simpa only [sq_abs] using (sq_lt_sq₀ (norm_nonneg z) (abs_nonneg (Real.cos θ)))
  rw [e]
  exact ⟨And.right, fun hr => ⟨hr.trans_le hc,hr⟩⟩

theorem left_complement (θ : ℝ) : leftBasin θ = disk \ rightBasin θ := by
  ext z
  change (_ ∧ (1-‖z‖^2 ≤ 1-Real.cos θ ^ 2)) ↔ (_ ∧ ¬ (_ ∧ ‖z‖^2 < Real.cos θ ^ 2))
  constructor
  · rintro ⟨hz,hl⟩
    exact ⟨hz, fun hr => (not_lt_of_ge (by linarith)) hr.2⟩
  · rintro ⟨hz,hr⟩
    have e : ¬ ‖z‖^2 < Real.cos θ ^ 2 := fun h => hr ⟨hz,h⟩
    exact ⟨hz, by linarith⟩

private theorem disk_area (r : ℝ) (hr : 0 ≤ r) :
    volume.real (Metric.ball (0 : Plane) r) = r^2*Real.pi := by
  simp [measureReal_def, EuclideanSpace.volume_ball_fin_two, ENNReal.toReal_mul,
    ENNReal.toReal_pow, ENNReal.toReal_ofReal hr, ENNReal.toReal_ofReal Real.pi_pos.le]

theorem right_fraction (θ : ℝ) : fraction (rightBasin θ) = Real.cos θ ^ 2 := by
  rw [fraction, right_disk, disk_area _ (abs_nonneg _)]
  unfold disk
  rw [disk_area _ zero_le_one, sq_abs, one_pow, one_mul]
  exact mul_div_cancel_right₀ _ Real.pi_ne_zero

theorem left_fraction (θ : ℝ) : fraction (leftBasin θ) = Real.sin θ ^ 2 := by
  have sub : rightBasin θ ⊆ disk := fun _ hz => hz.1
  have meas : MeasurableSet (rightBasin θ) := by rw [right_disk]; exact measurableSet_ball
  have fin : volume disk ≠ ⊤ := by unfold disk; exact measure_ball_lt_top.ne
  have area : volume.real disk = Real.pi := by
    unfold disk
    rw [disk_area _ zero_le_one, one_pow, one_mul]
  have right := right_fraction θ
  unfold fraction at right ⊢
  rw [area, div_eq_iff Real.pi_ne_zero] at right
  rw [left_complement, measureReal_diff sub meas, area, right, div_eq_iff Real.pi_ne_zero]
  nlinarith [Real.sin_sq_add_cos_sq θ, congrArg (fun t : ℝ => t*Real.pi) (Real.sin_sq_add_cos_sq θ)]

/-- One setting-independent preparation with zero doubles, zero misses,
and the unconditional, normalized area of each readout basin. -/
theorem closed_history_witness (θ : ℝ) :
    (∀ z ∈ disk, returnMap (history z) = history z ∧ Operational (history z) (Real.cos θ ^ 2)) ∧
    fraction (rightBasin θ) = Real.cos θ ^ 2 ∧ fraction (leftBasin θ) = Real.sin θ ^ 2 :=
  ⟨fun z _ => ⟨disk_history_fixed z, fixed_operational _ (disk_history_fixed z) _⟩,
   right_fraction θ, left_fraction θ⟩

#print axioms closed_history_witness
#print axioms all_settings_iff_fixed
#print axioms failure_length
end BAM.ClosedHistory
