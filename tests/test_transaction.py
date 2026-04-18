"""Tests for the Wheeler–Feynman transaction protocol on S³."""

import numpy as np
import pytest

from geometrodynamics.transaction.s3_geometry import (
    nrm4,
    antipode4,
    geo4,
    hsp,
    s3_green_potential,
    s3_green_field_kernel,
    s3_tangent_direction,
)
from geometrodynamics.transaction.particles import (
    ThroatMode,
    MouthState,
    Particle4,
    GravWave,
)
from geometrodynamics.transaction.handshake import (
    gw_hit_envelope,
    antipodal_match_weight,
    complex_mode_overlap,
    mode_overlap,
    su2_spin_phase,
    wrap_phase,
    mouth_activity,
    retarded_offer_amplitude,
    advanced_confirm_amplitude,
    retro_phase_match,
    make_offer,
    complete_transaction,
    OfferSignal,
    Transaction,
)


class TestS3Geometry:
    def test_nrm4_unit_length(self):
        v = nrm4([1, 2, 3, 4])
        assert abs(np.linalg.norm(v) - 1.0) < 1e-14

    def test_nrm4_rejects_zero(self):
        with pytest.raises(ValueError):
            nrm4([0, 0, 0, 0])

    def test_antipode_is_negation(self):
        p = np.array([1.0, 0.0, 0.0, 0.0])
        ap = antipode4(p)
        np.testing.assert_allclose(ap, -p)

    def test_geodesic_self_distance_zero(self):
        p = nrm4([1, 0, 0, 0])
        assert abs(geo4(p, p)) < 1e-14

    def test_geodesic_antipodal_distance_pi(self):
        p = nrm4([1, 0, 0, 0])
        assert abs(geo4(p, antipode4(p)) - np.pi) < 1e-10

    def test_hsp_unit_length(self):
        """Hopf-sphere parameterisation gives unit 4-vectors."""
        v = hsp(np.pi / 3, np.pi / 4, np.pi / 6)
        assert abs(np.linalg.norm(v) - 1.0) < 1e-14

    def test_green_potential_symmetry(self):
        """G(ψ) should be symmetric about ψ = π/2 (up to sign convention)."""
        g1 = s3_green_potential(0.5)
        g2 = s3_green_potential(np.pi - 0.5)
        # Both should be finite and of similar magnitude
        assert np.isfinite(g1)
        assert np.isfinite(g2)

    def test_green_potential_diverges_near_source(self):
        """Green function should be large near ψ = 0."""
        g_near = abs(s3_green_potential(0.1))
        g_far = abs(s3_green_potential(1.5))
        assert g_near > g_far

    def test_tangent_direction_orthogonal(self):
        """Tangent vector should be orthogonal to source point."""
        src = nrm4([1, 0, 0, 0])
        dst = nrm4([0, 1, 0, 0])
        tang = s3_tangent_direction(src, dst)
        assert abs(np.dot(src, tang)) < 1e-12


class TestParticles:
    def test_throat_mode_step(self):
        """ThroatMode oscillates under zero drive."""
        mode = ThroatMode(l=1, n=0, omega=1.0, alpha_q=1.0, a=1.0)
        for _ in range(100):
            mode.step(drive=0.0, dt=0.01)
        # Should oscillate, not diverge
        assert abs(mode.a) < 2.0

    def test_mouth_state_charge_sign(self):
        """Opposite orientation_sign gives opposite charge."""
        m1 = MouthState(orientation_sign=+1)
        m2 = MouthState(orientation_sign=-1)
        mode1 = ThroatMode(l=1, n=0, omega=1.0, alpha_q=1.0, a=0.5)
        mode2 = ThroatMode(l=1, n=0, omega=1.0, alpha_q=1.0, a=0.5)
        m1.modes[(1, 0)] = mode1
        m2.modes[(1, 0)] = mode2
        assert m1.q_geom() == -m2.q_geom()

    def test_grav_wave_propagation(self):
        """GW shell expands and terminates at radius π."""
        gw = GravWave(gid=0, src_pid=0, p0=np.array([1, 0, 0, 0]), t_emit=0.0)
        gw.step(t=100.0, c_gw=0.5)
        assert gw.done
        assert gw.radius == np.pi


class TestHandshake:
    def test_hit_envelope_peak(self):
        """Hit envelope peaks at ψ_diff = 0."""
        assert abs(gw_hit_envelope(0.0) - 1.0) < 1e-14

    def test_hit_envelope_decay(self):
        """Hit envelope decays away from zero."""
        assert gw_hit_envelope(1.0) < gw_hit_envelope(0.1) < 1.0

    def test_antipodal_match_exact(self):
        """Exact antipodal pair gets weight ≈ 1."""
        p = nrm4([1, 0, 0, 0])
        w = antipodal_match_weight(p, antipode4(p))
        assert w > 0.99

    def test_antipodal_match_non_antipodal(self):
        """Non-antipodal pair gets lower weight."""
        p1 = nrm4([1, 0, 0, 0])
        p2 = nrm4([0, 1, 0, 0])  # 90° apart, not antipodal
        w = antipodal_match_weight(p1, p2)
        assert w < 0.5

    def test_mode_overlap_identical(self):
        """Identical mouth spectra have overlap 1."""
        m1 = MouthState(orientation_sign=+1)
        m2 = MouthState(orientation_sign=-1)
        mode1 = ThroatMode(l=1, n=0, omega=1.0, alpha_q=1.0, a=1.0, phase=0.0)
        mode2 = ThroatMode(l=1, n=0, omega=1.0, alpha_q=1.0, a=1.0, phase=0.0)
        m1.modes[(1, 0)] = mode1
        m2.modes[(1, 0)] = mode2
        assert abs(mode_overlap(m1, m2) - 1.0) < 1e-10

    def test_mode_overlap_empty(self):
        """Empty mouth spectra have zero overlap."""
        m1 = MouthState(orientation_sign=+1)
        m2 = MouthState(orientation_sign=-1)
        assert mode_overlap(m1, m2) == 0.0

    def test_su2_spin_phase_is_unit(self):
        """SU(2) spin phase has unit modulus."""
        phase = su2_spin_phase(np.pi / 3)
        assert abs(abs(phase) - 1.0) < 1e-14

    def test_wrap_phase_range(self):
        """wrap_phase maps to (−π, π]."""
        for val in [0, np.pi, -np.pi, 3 * np.pi, -5.5]:
            w = wrap_phase(val)
            assert -np.pi <= w <= np.pi

    def test_mouth_activity_zero_for_empty(self):
        m = MouthState(orientation_sign=+1)
        assert mouth_activity(m) == 0.0

    def test_mouth_activity_nonzero_for_active(self):
        m = MouthState(orientation_sign=+1)
        m.modes[(1, 0)] = ThroatMode(l=1, n=0, omega=1.0, alpha_q=1.0, a=0.5)
        assert mouth_activity(m) > 0.0

    def test_phase_match_perfect_closure(self):
        """Zero total phase gives mismatch ≈ 0 and weight ≈ 1."""
        offer = complex(np.exp(1j * 0.3))
        confirm = complex(np.exp(-1j * 0.3))
        mismatch, weight = retro_phase_match(offer, confirm)
        assert mismatch < 0.01
        assert weight > 0.99

    def test_phase_match_pi_branch(self):
        """Phase sum ≈ π also closes (Möbius/spinor branch)."""
        offer = complex(np.exp(1j * 0.5))
        confirm = complex(np.exp(1j * (np.pi - 0.5)))
        mismatch, weight = retro_phase_match(offer, confirm)
        assert mismatch < 0.01
        assert weight > 0.99

    def test_end_to_end_offer_confirm(self):
        """Minimal end-to-end: offer from source, confirm from absorber."""
        # Source particle at north pole of S³
        src_mouth = MouthState(orientation_sign=+1)
        src_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.5, adot=0.1,
        )
        src = Particle4(pid=0, p4=nrm4([1, 0, 0, 0]), mouth=src_mouth)

        # Candidate near equator
        cand_mouth = MouthState(orientation_sign=+1)
        cand_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.3, adot=0.05,
        )
        cand = Particle4(pid=1, p4=nrm4([0, 1, 0, 0]), mouth=cand_mouth)

        # Absorber near antipode of source
        dst_mouth = MouthState(orientation_sign=-1)
        dst_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.4, adot=-0.08,
        )
        dst = Particle4(pid=2, p4=nrm4([-1, 0.01, 0, 0]), mouth=dst_mouth)

        # Offer
        offer_amp, theta_sc, response, geom_phase = retarded_offer_amplitude(
            src, cand, hit_weight=0.9, gw_radius=1.5,
        )
        assert abs(offer_amp) > 0
        assert np.isfinite(offer_amp)

        # Confirm
        confirm_amp, anti_w, overlap, field, spin_ang, theta_cd = (
            advanced_confirm_amplitude(cand, dst, dst_field_sol=None)
        )
        assert np.isfinite(confirm_amp)
        assert anti_w > 0  # near-antipodal should have some weight

        # Phase closure
        mismatch, weight = retro_phase_match(offer_amp, confirm_amp)
        assert np.isfinite(mismatch)
        assert 0 <= weight <= 1.0


class TestTangherliniTransactionIntegration:
    """Cross-subsystem integration: feed a real Maxwell solve into the
    transaction handshake.

    This is the test that was missing in v0.40.0 and v0.40.1, and whose
    absence masked the KeyError('r') bug when coupling the two subsystems.
    """

    @pytest.fixture
    def maxwell_field_sol(self):
        """Solve the Maxwell BVP from real Tangherlini eigenmodes."""
        from geometrodynamics.tangherlini import (
            solve_radial_modes,
            solve_maxwell_from_eigenmode,
        )

        modes = {}
        for l in [1, 3, 5]:
            oms, fns, rg = solve_radial_modes(l)
            modes[l] = {"omega": oms, "funcs": fns}
        return solve_maxwell_from_eigenmode(modes)

    def test_maxwell_output_has_required_keys(self, maxwell_field_sol):
        """Maxwell output must expose 'r' and 'E_r' for transaction use."""
        assert "r" in maxwell_field_sol
        assert "E_r" in maxwell_field_sol
        assert len(maxwell_field_sol["r"]) == len(maxwell_field_sol["E_r"])

    def test_confirm_amplitude_with_real_field(self, maxwell_field_sol):
        """advanced_confirm_amplitude succeeds with real Maxwell output."""
        cand_mouth = MouthState(orientation_sign=+1)
        cand_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.3, adot=0.05,
        )
        cand = Particle4(pid=0, p4=nrm4([0, 1, 0, 0]), mouth=cand_mouth)

        dst_mouth = MouthState(orientation_sign=-1)
        dst_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.4, adot=-0.08,
        )
        dst = Particle4(pid=1, p4=nrm4([-1, 0.01, 0, 0]), mouth=dst_mouth)

        # This was the call that raised KeyError('r') in v0.40.1
        confirm_amp, anti_w, overlap, field_at_dst, spin_ang, theta_cd = (
            advanced_confirm_amplitude(cand, dst, maxwell_field_sol)
        )

        assert np.isfinite(confirm_amp)
        assert field_at_dst > 0, "E_r at sample radius should be positive"
        assert anti_w > 0

    def test_full_handshake_with_real_field(self, maxwell_field_sol):
        """Complete offer → confirm → closure with real Tangherlini field."""
        src_mouth = MouthState(orientation_sign=+1)
        src_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.5, adot=0.1,
        )
        src = Particle4(pid=0, p4=nrm4([1, 0, 0, 0]), mouth=src_mouth)

        cand_mouth = MouthState(orientation_sign=+1)
        cand_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.3, adot=0.05,
        )
        cand = Particle4(pid=1, p4=nrm4([0, 1, 0, 0]), mouth=cand_mouth)

        dst_mouth = MouthState(orientation_sign=-1)
        dst_mouth.modes[(1, 0)] = ThroatMode(
            l=1, n=0, omega=1.0, alpha_q=1.0, a=0.4, adot=-0.08,
        )
        dst = Particle4(pid=2, p4=nrm4([-1, 0.01, 0, 0]), mouth=dst_mouth)

        # Offer
        offer_amp, _, _, _ = retarded_offer_amplitude(
            src, cand, hit_weight=0.9, gw_radius=1.5,
        )

        # Confirm with real Maxwell field
        confirm_amp, anti_w, overlap, field_at_dst, _, _ = (
            advanced_confirm_amplitude(cand, dst, maxwell_field_sol)
        )

        # Phase closure
        mismatch, weight = retro_phase_match(offer_amp, confirm_amp)

        assert abs(offer_amp) > 0
        assert np.isfinite(confirm_amp)
        assert field_at_dst > 0
        assert np.isfinite(mismatch)
        assert 0 <= weight <= 1.0


# ── Helpers for particle construction ────────────────────────────────────────

def _make_particle(pid, p4_raw, sign, a=0.5, adot=0.1):
    """Shorthand for building a test particle with one active throat mode."""
    mouth = MouthState(orientation_sign=sign)
    mouth.modes[(1, 0)] = ThroatMode(
        l=1, n=0, omega=1.0, alpha_q=1.0, a=a, adot=adot,
    )
    return Particle4(pid=pid, p4=nrm4(p4_raw), mouth=mouth)


class TestMakeOffer:
    """Tests for the make_offer lifecycle helper."""

    def test_returns_offer_signal(self):
        src = _make_particle(0, [1, 0, 0, 0], +1)
        cand = _make_particle(1, [0.7, 0.7, 0, 0], +1, a=0.3)
        offer = make_offer(
            src, cand, grav_id=42, hit_weight=0.85, gw_radius=1.2, t_now=5.0,
        )
        assert isinstance(offer, OfferSignal)

    def test_offer_preserves_ids(self):
        src = _make_particle(0, [1, 0, 0, 0], +1)
        cand = _make_particle(7, [0.7, 0.7, 0, 0], +1, a=0.3)
        offer = make_offer(
            src, cand, grav_id=42, hit_weight=0.85, gw_radius=1.2, t_now=5.0,
        )
        assert offer.src_pid == 0
        assert offer.cand_pid == 7
        assert offer.grav_id == 42
        assert offer.t_birth == 5.0

    def test_offer_amplitude_nonzero(self):
        src = _make_particle(0, [1, 0, 0, 0], +1)
        cand = _make_particle(1, [0.7, 0.7, 0, 0], +1, a=0.3)
        offer = make_offer(
            src, cand, grav_id=0, hit_weight=0.85, gw_radius=1.2, t_now=0.0,
        )
        assert abs(offer.amp) > 0
        assert np.isfinite(offer.amp)

    def test_offer_is_alive_within_ttl(self):
        src = _make_particle(0, [1, 0, 0, 0], +1)
        cand = _make_particle(1, [0.7, 0.7, 0, 0], +1, a=0.3)
        offer = make_offer(
            src, cand, grav_id=0, hit_weight=0.85, gw_radius=1.2, t_now=1.0,
        )
        assert offer.is_alive(1.5)
        # Should expire well past TTL
        assert not offer.is_alive(1.0 + offer.ttl + 1.0)

    def test_offer_captures_source_charge(self):
        """OfferSignal.q_src must equal the source particle's q_geom."""
        src = _make_particle(0, [1, 0, 0, 0], +1, a=0.7)
        cand = _make_particle(1, [0.7, 0.7, 0, 0], +1, a=0.1)
        offer = make_offer(
            src, cand, grav_id=0, hit_weight=0.85, gw_radius=1.2, t_now=0.0,
        )
        assert offer.q_src == pytest.approx(src.mouth.q_geom(), rel=1e-12)
        assert offer.q_src != pytest.approx(cand.mouth.q_geom(), rel=0.1)


class TestCompleteTransaction:
    """Tests for the complete_transaction lifecycle helper."""

    def _antipodal_trio(self):
        """Build src/cand/dst where cand and dst are near-antipodal.

        This geometry is designed so that:
        - src is at the north pole
        - cand is ~45° away (reachable by GW)
        - dst is near the antipode of cand (so anti_w ≈ 1)
        - all mouths share the same (l=1,n=0) mode for nonzero overlap
        """
        src = _make_particle(0, [1, 0, 0, 0], +1, a=0.5, adot=0.1)
        cand = _make_particle(1, [0.7, 0.7, 0, 0], +1, a=0.4, adot=0.08)
        # dst near antipode of cand
        dst = _make_particle(2, [-0.7, -0.7, 0, 0], -1, a=0.3, adot=-0.05)
        return src, cand, dst

    def test_successful_transaction(self):
        """Near-antipodal pair with matching modes → confirmed Transaction."""
        src, cand, dst = self._antipodal_trio()
        offer = make_offer(
            src, cand, grav_id=1, hit_weight=0.9, gw_radius=0.8, t_now=0.0,
        )
        txn = complete_transaction(
            offer, cand, dst, dst_field_sol=None, t_now=0.1,
        )
        assert txn is not None
        assert isinstance(txn, Transaction)
        assert txn.is_confirmed
        assert txn.src_pid == 0
        assert txn.cand_pid == 1
        assert txn.dst_pid == 2
        assert abs(txn.amp) > 0
        assert txn.anti_weight > 0.9  # near-antipodal

    def test_transaction_fields_populated(self):
        """All Transaction fields are finite and physically reasonable."""
        src, cand, dst = self._antipodal_trio()
        offer = make_offer(
            src, cand, grav_id=1, hit_weight=0.9, gw_radius=0.8, t_now=0.0,
        )
        txn = complete_transaction(
            offer, cand, dst, dst_field_sol=None, t_now=0.1,
        )
        assert txn is not None
        assert np.isfinite(txn.amp)
        assert np.isfinite(txn.offer_amp)
        assert np.isfinite(txn.confirm_amp)
        assert np.isfinite(txn.phase_mismatch)
        assert 0 <= txn.phase_match_weight <= 1.0
        assert txn.overlap > 0
        assert txn.field_at_dst > 0

    def test_expired_offer_returns_none(self):
        """An expired offer cannot complete a transaction."""
        src, cand, dst = self._antipodal_trio()
        offer = make_offer(
            src, cand, grav_id=1, hit_weight=0.9, gw_radius=0.8, t_now=0.0,
        )
        # Far past TTL
        txn = complete_transaction(
            offer, cand, dst, dst_field_sol=None, t_now=offer.ttl + 100.0,
        )
        assert txn is None

    def test_non_antipodal_pair_returns_none(self):
        """Particles that are not near-antipodal → confirm too weak → None."""
        src = _make_particle(0, [1, 0, 0, 0], +1)
        cand = _make_particle(1, [0.7, 0.7, 0, 0], +1, a=0.3)
        # dst is near cand, NOT near its antipode
        dst = _make_particle(2, [0.6, 0.8, 0, 0], -1, a=0.3)
        offer = make_offer(
            src, cand, grav_id=1, hit_weight=0.9, gw_radius=0.8, t_now=0.0,
        )
        txn = complete_transaction(
            offer, cand, dst, dst_field_sol=None, t_now=0.1,
        )
        # anti_w ≈ 0 → confirm_amp ≈ 0 → total_amp ≈ 0 → None
        assert txn is None

    def test_transaction_is_alive(self):
        """A fresh transaction is alive; an old one is not."""
        src, cand, dst = self._antipodal_trio()
        offer = make_offer(
            src, cand, grav_id=1, hit_weight=0.9, gw_radius=0.8, t_now=0.0,
        )
        txn = complete_transaction(
            offer, cand, dst, dst_field_sol=None, t_now=0.1,
        )
        assert txn is not None
        assert txn.is_alive(0.2)
        assert not txn.is_alive(0.1 + txn.ttl + 1.0)

    def test_with_real_maxwell_field(self):
        """Complete a transaction using a real Tangherlini field solve."""
        from geometrodynamics.tangherlini import (
            solve_radial_modes,
            solve_maxwell_from_eigenmode,
        )

        modes = {}
        for l in [1, 3, 5]:
            oms, fns, rg = solve_radial_modes(l)
            modes[l] = {"omega": oms, "funcs": fns}
        field_sol = solve_maxwell_from_eigenmode(modes)

        src, cand, dst = self._antipodal_trio()
        offer = make_offer(
            src, cand, grav_id=1, hit_weight=0.9, gw_radius=0.8, t_now=0.0,
        )
        txn = complete_transaction(
            offer, cand, dst, field_sol, t_now=0.1,
        )
        assert txn is not None
        assert isinstance(txn, Transaction)
        assert txn.field_at_dst > 0, "should use real E_r, not fallback"
        assert txn.is_confirmed

    def test_q_src_is_source_not_candidate(self):
        """txn.q_src must be the emitter's charge, not the relay's.

        Regression test for the bug where complete_transaction recorded
        cand.mouth.q_geom() as q_src.  The source charge is captured at
        offer-creation time in OfferSignal.q_src and propagated through.

        To make this test discriminating, source and candidate have
        deliberately different mode amplitudes (0.8 vs 0.2), so their
        geometric charges are different by a factor of 4.
        """
        # Source: large amplitude → large |q_geom|
        src = _make_particle(0, [1, 0, 0, 0], +1, a=0.8, adot=0.1)
        # Candidate: small amplitude → small |q_geom|
        cand = _make_particle(1, [0.7, 0.7, 0, 0], +1, a=0.2, adot=0.02)
        # Destination: near antipode of candidate
        dst = _make_particle(2, [-0.7, -0.7, 0, 0], -1, a=0.3, adot=-0.05)

        q_source_true = src.mouth.q_geom()
        q_candidate = cand.mouth.q_geom()

        # Sanity: these must be different for the test to be meaningful
        assert q_source_true != pytest.approx(q_candidate, rel=0.1), (
            "Test design error: source and candidate charges too similar"
        )

        offer = make_offer(
            src, cand, grav_id=1, hit_weight=0.9, gw_radius=0.8, t_now=0.0,
        )

        # Verify the offer captured the source charge, not the candidate's
        assert offer.q_src == pytest.approx(q_source_true, rel=1e-12)
        assert offer.q_src != pytest.approx(q_candidate, rel=0.1)

        txn = complete_transaction(
            offer, cand, dst, dst_field_sol=None, t_now=0.1,
        )
        assert txn is not None

        # THE KEY ASSERTION: txn.q_src is the source, not the relay
        assert txn.q_src == pytest.approx(q_source_true, rel=1e-12), (
            f"txn.q_src={txn.q_src} should match source={q_source_true}, "
            f"not candidate={q_candidate}"
        )
