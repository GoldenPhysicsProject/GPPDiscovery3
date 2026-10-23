import GppVerify.RiemannHypothesis.VariationalSchurScalar
import Mathlib.Tactic

/-!
Candidate extension for Claude's review. Not locally Lean-compiled.
An arbitrary predictor needs no Galerkin equation. No arithmetic residual bound is asserted.
-/
namespace GppResidualCertificateCandidate
variable {V : Type*} [AddCommGroup V] [Module ℝ V]

theorem translate (a : LinearMap.BilinForm ℝ V)
    (hs : ∀ u v, a u v = a v u) (ℓ : V →ₗ[ℝ] ℝ) (x v : V) :
    2 * ℓ (x + v) - a (x + v) (x + v) =
    (2 * ℓ x - a x x) + (2 * (ℓ v - a x v) - a v v) := by
  simp only [map_add, LinearMap.add_apply]
  rw [hs v x]
  ring

theorem bound_from_residual (a : LinearMap.BilinForm ℝ V)
    (hs : ∀ u v, a u v = a v u) (ℓ : V →ₗ[ℝ] ℝ) (x : V) (e : ℝ)
    (hr : ∀ v, 2 * (ℓ v - a x v) - a v v ≤ e) :
    ∀ v, 2 * ℓ v - a v v ≤ (2 * ℓ x - a x x) + e := by
  intro v
  have hid := translate a hs ℓ x (v - x)
  have hv : x + (v - x) = v := by abel
  rw [hv] at hid
  linarith [hr (v - x)]

theorem bound_from_squared_residual (a : LinearMap.BilinForm ℝ V)
    (hs : ∀ u v, a u v = a v u) (ha : ∀ v, 0 ≤ a v v)
    (ℓ : V →ₗ[ℝ] ℝ) (x : V) (e : ℝ) (he : 0 ≤ e)
    (hr : ∀ v, (ℓ v - a x v)^2 ≤ e * a v v) :
    ∀ v, 2 * ℓ v - a v v ≤ (2 * ℓ x - a x x) + e := by
  apply bound_from_residual a hs ℓ x e
  have hb := (GppVariationalSchurScalar.bound_iff a ha (ℓ - a x) e).2
    ⟨he, by simpa using hr⟩
  simpa using hb

end GppResidualCertificateCandidate
